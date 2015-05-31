# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import zip
from builtins import object
import numpy as np
import gzip
import string
import ftplib

from fastcache import clru_cache
from tempfile import NamedTemporaryFile
from collections import defaultdict
import six

import Bio
from Bio.PDB import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from . import errors
from . import helper_functions as hf


A_DICT = {
    'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU',
    'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS',
    'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP',
    'Y':'TYR', 'V':'VAL', 'U':'SEC', 'O':'PYL',
    'B':'ASX', 'Z':'GLX', 'J':'XLE', 'X':'XAA', '*':'TER'
}
AAA_DICT = dict([(value,key) for key,value in list(A_DICT.items())])
AAA_DICT['UNK'] = 'X'
AAA_DICT['MSE'] = 'M'
AAA_DICT['CSD'] = 'C'
# Phosphorylated residues
#AAA_DICT['SEP'] = 'S' # PHOSPHOSERINE
#AAA_DICT['TPO'] = 'T' # PHOSPHOTHREONINE
#AAA_DICT['SEP'] = 'Y' # O-PHOSPHOTYROSINE

# Methylated lysines
AAA_DICT['MLZ'] = 'K'
AAA_DICT['MLY'] = 'K'
AAA_DICT['M3L'] = 'K'
methylated_lysines = ['MLZ', 'MLY', 'M3L']
amino_acids = list(AAA_DICT.keys())
lysine_atoms = ['N', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'C', 'O']



#%%################################################################################################
### Functions for downloading and parsing pdb files

class MMCIFParserMod(MMCIFParser):
    def __init__(self, tmp_path):
        self.tmp_path = tmp_path

    def get_structure(self, structure_id, gzip_fh):
        """
        Altered ``get_structure`` method which accepts gzip file-handles as input.
        """
        with NamedTemporaryFile(mode='w', dir=self.tmp_path) as temp_fh:
            temp_fh.writelines(gzip_fh.readlines())
            temp_fh.flush()
            temp_fh.seek(0)
            return super(MMCIFParserMod, self).get_structure(structure_id, temp_fh.name)



def download_pdb(pdb_id, pdb_path_suffix, tmp_path='/tmp/'):
    """
    Download the specified pdb file from the internet

    Parameters
    ----------
    pdb_id : str
        Structure ID of the pdb file
    pdb_path_suffix : str
        Path leading from the divided pdb dir to the pdb file
    tmp_path : str, optional
        Folder to store the downloaded file
    """
    ftp_hostname = 'ftp.wwpdb.org'
    ftp_pdb_path = '/pub/pdb/data/structures/divided/pdb/'
    temp_filename = tmp_path + pdb_path_suffix.split('/')[-1]

    # Log into server
    ftp = ftplib.FTP()
    ftp.connect(ftp_hostname)
    try:
        ftp.login()
        ftp.retrbinary("RETR {}".format(ftp_pdb_path + pdb_path_suffix), open(temp_filename, 'wb').write)
    except:
        raise
    finally:
        ftp.quit()

    return temp_filename



def _get_pdb_path_suffix(pdb_id, pdb_type='ent'):
    pdb_path_suffix = ''
    if pdb_type == 'ent':
        # Original PDB structure
        prefix = 'pdb'
        suffix = '.ent.gz'

    elif pdb_type == 'cif':
        # mmCIF pdb structure
        prefix = ''
        suffix = '.cif.gz'
        pdb_path_suffix += '../mmCIF/'

    elif pdb_type == 'pdb':
        # The first biological unit
        prefix = ''
        suffix = '.pdb1.gz'
        pdb_path_suffix += '../../../biounit/coordinates/divided/'

    pdb_path_suffix += pdb_id[1:3].lower() + '/' + prefix + pdb_id.lower() + suffix
    return pdb_path_suffix


@clru_cache(maxsize=128, typed=False)
def get_pdb(pdb_id, pdb_path, tmp_path='/tmp/', pdb_type='ent', use_external=True):
    """
    Parse a pdb file with biopythons PDBParser() and return the structure.

    Parameters
    ----------
    pdb_code : str
        Four letter code of the PDB file
    pdb_path : str
        Biopython pdb structure
    tmp_path : str, optional, default='/tmp/'
        Path to the folder for storing temporary files
    pdb_type : 'ent'/'pdb'/'cif', optional, default='ent'
        The extension of the pdb to use

    Raises
    ------
    PDBNotFoundError
        If the pdb file could not be retrieved from the local (and remote) databases
    """
    pdb_path_suffix = _get_pdb_path_suffix(pdb_id, pdb_type)

    if pdb_type in {'ent', 'pdb'}:
        parser = PDBParser(QUIET=True)
    elif pdb_type in {'cif'}:
        parser = MMCIFParserMod(tmp_path=tmp_path)
    else:
        raise Exception('Unsupported ``pdb_type`` {}'.format(pdb_type))

    try:
        structure = parser.get_structure(pdb_id, gzip.open(pdb_path + pdb_path_suffix, 'rt'))
    except IOError:
        error_message = (
            'PDB not found! (pdb_id: {}, pdb_path: {}, pdb_type: {})'
            .format(pdb_id, pdb_path, pdb_type)
        )
        if not use_external:
            raise errors.PDBNotFoundError(error_message)
        print('Retrieving pdb from the wwpdb ftp server...')
        temp_filename = download_pdb(pdb_id, pdb_path_suffix, tmp_path)
        structure = parser.get_structure(pdb_id, gzip.open(temp_filename, 'rt'))

    return structure


#%%
def euclidean_distance(a, b):
    """Calculate the Euclidean distance between two lists or tuples of arbitrary length.
    """
    return np.sqrt(sum((a - b)**2 for a, b in zip(a, b)))


def calculate_distance(atom_1, atom_2, cutoff=None):
    """Calculate the distance between two points in 3D space.

    Parameters
    ----------
    cutoff : float, optional
        The maximum distance allowable between two points.
    """

    if ((type(atom_1) == type(atom_2) == list) or
        (type(atom_1) == type(atom_2) == tuple)):
            a = atom_1
            b = atom_2
    elif hasattr(atom_1, 'coord') and hasattr(atom_2, 'coord'):
            a = atom_1.coord
            b = atom_2.coord
    else:
        raise Exception('Unsupported format {} {}'.format(type(atom_1), type(atom_2)))

    assert(len(a) == 3 and len(b) == 3)
    if cutoff is None or all(abs(p - q) <= cutoff for p, q in zip(a, b)):
        return euclidean_distance(a, b)


#%%
def get_chain_seqres_sequence(chain, aa_only=False):
    """Get the amino acid sequence for the construct coding for the given chain.

    Extracts a sequence from a PDB file. Usefull when interested in the
    sequence that was used for crystallization and not the ATOM sequence.

    Parameters
    ----------
    aa_only : bool
        If aa_only is set to `False`, selenomethionines will be included in the sequence.
        See: http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html.
    """
    sequence = Seq('', IUPAC.protein)
    for pb in PPBuilder().build_peptides(chain, aa_only=aa_only):
        sequence += sequence + pb.get_sequence()
    return sequence



def get_chain_sequence_and_numbering(chain, domain_def_tuple=None, include_hetatms=False):
    """Get the amino acid sequence and a list of residue ids for the given chain.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        The chain for which to get the amino acid sequence and numbering.
    """
    if domain_def_tuple is not None:
        start_resid, end_resid = domain_def_tuple

    chain_numbering = []
    chain_numbering_extended = []
    chain_sequence = []
    inside_domain = False
    for res in chain:
        #
        resid = str(res.id[1]) + res.id[2].strip()

        if domain_def_tuple is None or resid == start_resid:
            inside_domain = True

        if inside_domain and (include_hetatms or res.resname in amino_acids):
            chain_numbering.append(res.id[1])
            chain_numbering_extended.append(resid)
            chain_sequence.append(AAA_DICT[res.resname])

        if domain_def_tuple is not None and resid == end_resid:
            inside_domain = False

    chain_sequence = ''.join(chain_sequence)
    return chain_sequence, chain_numbering_extended



def convert_position_to_resid(chain, positions, domain_def_tuple=None):
    """
    """
    __, chain_numbering = get_chain_sequence_and_numbering(chain, domain_def_tuple)
    return [chain_numbering[p-1] for p in positions]



def convert_resid_to_position(chain, resids, domain_def_tuple=None):
    """
    """
    __, chain_numbering = get_chain_sequence_and_numbering(chain, domain_def_tuple)
    return [chain_numbering.index(resid) for resid in resids]



def get_structure_sequences(file_or_structure, seqres_sequence=False):
    """
    Convenience function returining a dictionary of sequences for a given file
    or a Biopython Structure, Model or Chain.
    """
    if isinstance(file_or_structure, six.string_types):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('ID', file_or_structure)
        model = structure[0]
    elif isinstance(file_or_structure, Bio.PDB.Structure.Structure):
        model = file_or_structure[0]
    elif isinstance(file_or_structure, Bio.PDB.Model.Model):
        model = file_or_structure
    elif isinstance(file_or_structure, Bio.PDB.Chain.Chain):
        model = [file_or_structure]
    else:
        raise Exception(
            'Unexpected type {} for input ``file_or_structure`` {}!'
            .format(file_or_structure, type(file_or_structure)))

    chain_sequences = defaultdict(list)
    for chain in model:
        if seqres_sequence:
            chain_sequence = get_chain_seqres_sequence(chain)
        else:
            chain_sequence, __ = get_chain_sequence_and_numbering(chain)
        chain_sequences[chain.id] = chain_sequence
    return chain_sequences



#%%
def convert_aa(aa):
    """Convert amino acids from three letter code to one letter code or vice versa

    .. note:: Deprecated!
    
       Use ``''.join(AAA_DICT[aaa] for aaa in aa)`` and ``''.join(A_DICT[a] for a in aa)``.
    """
    if len(aa) == 3:
        try:
            return AAA_DICT[aa.upper()]
        except KeyError:
            print('Not a valid amino acid')
            return
    if len(aa) == 1:
        try:
            return A_DICT[aa.upper()]
        except KeyError:
            print('Not a valid amino acid')
            return
    print('Not a valid amino acid')


def convert_resnum_alphanumeric_to_numeric(resnum):
    """
    Convert residue numbering that has letters (i.e. 1A, 1B, 1C...) to
    residue numbering without letters (i.e. 1, 2, 3...).

    .. note:: Deprecated!
    
        Use ``get_chain_sequence_and_numbering()``.
    """
    idx_increment = 0
    while string.ascii_letters.find(resnum[-1]) != -1:
        idx_increment += string.ascii_letters.find(resnum[-1])
        resnum = resnum[:-1]
    resnum = int(resnum) + idx_increment
    return resnum




#%%################################################################################################

class SelectChains(Select):
    """Only accept the specified chains when saving.
    """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)



class PDBTemplate(object):
    """
    Attributes
    ----------
    domain_boundaries : list of lists of lists
        Elements in the outer list correspond to domains in each chain of the
        pdb. Elements of the inner list contain the start and end of each
        fragment of each domain. For example, if there is only one chain
        with pdb domain boundaries 1-10:20-45, this would correspond to
        domain_boundaries [[[1,10],[20,45]]].
    
    """
    def __init__(self, pdb_path, pdb_id, chain_ids, domain_defs, output_path, tmp_path, logger):
        """
        Parameters
        ----------
        
        """
        self.pdb_path = pdb_path
        self.pdb_id = pdb_id
        self.output_path = output_path
        self.structure = get_pdb(self.pdb_id, self.pdb_path, tmp_path).copy()
        self.chain_ids = (
            chain_ids if chain_ids else [chain.id for chain in self.structure[0].child_list]
        )
        self.logger = logger
        self.r_cutoff = 6 # remove hetatms more than x A away from the main chain(s)

        self.domain_boundaries = []
        for domain_def in domain_defs:
            self.domain_boundaries.append(hf.decode_domain_def(domain_def, merge=False, return_string=True))

        self.unique_id = (
            'pdb_id: {}, pdb_chain: {}, pdb_domain_def: {}'.format(
            self.pdb_id, self.chain_ids, self.domain_boundaries) )


    def extract(self):
        """Extract the wanted chains out of the PDB file. Removes water atoms
        and selects the domain regions (i.e. selects only those parts of the
        domain that are within the domain boundaries specified).
        """
        self.logger.debug('Extracting {}...'.format(self.unique_id))

        model = self.structure[0] # assuming that model 0 is always the desired one
        new_structure = Bio.PDB.Structure.Structure(self.pdb_id)
        new_model = Bio.PDB.Model.Model(0)
        hetatm_chain = Bio.PDB.Chain.Chain('Z')

        # Loop over every chain and every residue and make sure that everything is ok
        for chain_idx, chain_id in enumerate(self.chain_ids):
            chain = model[chain_id]
            res_idx = 0

            chain_numbering, domain_start_idxs, domain_end_idxs = (
                self._get_domain_def_idxs_for_chain(chain, chain_idx)
            )
            self.logger.debug(
                'domain_def: {}, domain_start_idxs: {}, domain_end_idxs: {}'.format(
                self.domain_boundaries[chain_idx], domain_start_idxs, domain_end_idxs))

            while res_idx < len(chain):
                res = chain.child_list[res_idx]
                original_res_id = res.id

                # Move water to the hetatm chain
                if res.id[0] == 'W':
                    self._move_hetatm_to_hetatm_chain(chain, hetatm_chain, res, echo=False)
                    continue

#                # Move heteroatoms to the hetatm chain
#                if res.id[0] != ' ':
#                    self._move_hetatm_to_hetatm_chain(chain, hetatm_chain, res, echo=True)
#                    continue
                
                # Now treating all unusual amino acids as hetatms
                # Convert methylated lysines to regular lysines
                if res.resname in methylated_lysines:
                    self._correct_methylated_lysines(res)

                # Move hetatms to the hetatm chain
                if res.resname not in amino_acids:
                    self._move_hetatm_to_hetatm_chain(chain, hetatm_chain, res)
                    continue

                # Cut each chain to domain boundaries
                residue_is_outside_domain = (
                    self._residue_outside_domain(
                        chain, chain_numbering, domain_start_idxs, domain_end_idxs, res)
                )
                if residue_is_outside_domain:
                    chain.detach_child(original_res_id)
                    continue

                res_idx += 1
            new_model.add(chain)

        # Make sure that the new model is not empty
        if not list(new_model.get_atoms()):
            raise errors.PDBEmptySequenceError(self.unique_id)

        # Remove hetatms if they are > 6A away from the chains of interest
        self._remove_distant_hatatms(new_model, hetatm_chain)

        if hetatm_chain:
            self.logger.debug('Adding hetatm chain of length {}'.format(len(hetatm_chain)))
            new_model.add(hetatm_chain)

        # If the hetatm chain is not empty, add it to the model
        new_structure.add(new_model)
        self.structure = new_structure

        self.logger.debug('PDB {} extracted successfully\n'.format(self.pdb_id))


    def get_chain_sequence_and_numbering(self, chain_id, *args, **varargs):
        """Call ``get_chain_sequence_and_numbering`` using chain with id ``chain_id``
        """
        chain = self.structure[0][chain_id]
        return get_chain_sequence_and_numbering(chain, *args, **varargs)


    def get_chain_seqres_sequence(self, chain_id, *args, **varargs):
        """Call ``get_chain_seqres_sequence`` using chain with id ``chain_id``
        """
        chain = self.structure[0][chain_id]
        return get_chain_seqres_sequence(chain, *args, **varargs)


    def _get_domain_def_idxs_for_chain(self, chain, chain_idx):
        if not self.domain_boundaries or not self.domain_boundaries[chain_idx]:
            return None, None

        __, chain_numbering = get_chain_sequence_and_numbering(chain)
        try:
            domain_start_idxs, domain_end_idxs = [
                tuple(chain_numbering.index(resid) for resid in resids)
                for resids in zip(*self.domain_boundaries[chain_idx])]
        except Exception as e:
            print(str(e))
            raise errors.PDBDomainDefsError(self.unique_id)

        return chain_numbering, domain_start_idxs, domain_end_idxs


    def _correct_methylated_lysines(self, res):
        new_resname = 'LYS'
        new_resid = (' ', res.id[1], res.id[2])
        self.logger.debug(
            'Renaming residue {} {} to {} {}'
            .format(res.resname, res.id, new_resname, new_resid))
        res.resname = new_resname
        res.id = new_resid
        atom_idx = 0
        while atom_idx < len(res):
            atom_id = res.child_list[atom_idx].id
            if atom_id not in lysine_atoms:
                self.logger.debug('Removing atom {} from residue {} {}.'.format(atom_id, res.resname, res.id))
                res.detach_child(atom_id)
            else:
                atom_idx += 1


    def _move_hetatm_to_hetatm_chain(self, chain, hetatm_chain, res, echo=False):
        self.logger.debug('Moving hetatm residue {} {} to the hetatm chain'.format(res.resname, res.id))
        chain.detach_child(res.id)
        hetatm_res = res
        hetatm_res.id = (hetatm_res.id[0], len(hetatm_chain)+1, hetatm_res.id[2],)
        hetatm_chain.add(hetatm_res)


    def _residue_outside_domain(self, chain, chain_numbering, domain_start_idxs, domain_end_idxs, res):
        """Returns `True` if residue ``res`` is outside the domain
        """
        if domain_start_idxs is None or domain_end_idxs is None:
            return False

        resid = str(res.id[1]) + res.id[2].strip()
        resid_idx = chain_numbering.index(resid)
        for domain_start_idx, domain_end_idx in zip(domain_start_idxs, domain_end_idxs):
            if resid_idx >= domain_start_idx and resid_idx <= domain_end_idx:
                # Residue is inside the domain
                return False

        # Residue is outside the domain
        return True


    def _remove_distant_hatatms(self, new_model, hetatm_chain):
        """Detach hetatms that are more than ``self.r_cutoff`` away from the main chain(s)
        """
        try:
            from Bio.PDB import NeighborSearch
        except ImportError as e:
            self.logger.error('Importing Biopython NeighborSearch failed with an error: {}'.format(e))
            raise Exception('Alternative not yet implemented!')
        ns = NeighborSearch(list(new_model.get_atoms()))
        hetatm_chain.id = [c for c in reversed(hf.uppercase) if c not in self.chain_ids][0]
        res_idx = 0
        while res_idx < len(hetatm_chain):
            res_1 = hetatm_chain.child_list[res_idx]
            in_contact = False
            for atom_1 in res_1:
                interacting_residues = ns.search(atom_1.get_coord(), self.r_cutoff, 'R')
                if interacting_residues:
                    # self.logger.debug(res_1.id)
                    # self.logger.debug(interacting_residues)
                    in_contact = True
            if in_contact:
                res_idx += 1
                continue
            # self.logger.debug('Detaching child: {}'.format(res_1.id))
            hetatm_chain.detach_child(res_1.id)


    def _unset_disordered_flags(self):
        """Change atom and residue ``disordered`` flag to `False`.
        Otherwise, Biopython may crash when saving the PDB structure.
        """
        self.logger.debug('Setting all residues and atoms marked as disorded to non-disorded')
        for m in self.structure:
            for c in m:
                for r in c:
                    if r.is_disordered() or r.disordered:
                        self.logger.debug(
                            'Changing disordered_flag on residue {} from {} to 0'
                            .format(r, r.disordered))
                        r.disordered = 0
                    for a in r:
                        if a.is_disordered() or a.disordered_flag:
                            self.logger.debug(
                                'Changing disordered_flag on atom {} from {} to 0'
                                .format(a, a.disordered_flag))
                            a.disordered_flag = 0


    def save_structure(self, remove_disordered=False):
        if remove_disordered:
            self._unset_disordered_flags()

        io = PDBIO()
        io.set_structure(self.structure)
        try:
            outFile = self.output_path + self.pdb_id + ''.join(self.chain_ids) + '.pdb'
            io.save(outFile)
            if len(self.chain_ids) > 1:
                for chain_id in self.chain_ids:
                    outFile = self.output_path + self.pdb_id + chain_id + '.pdb'
                    io.save(outFile, select=SelectChains(chain_id))
        except AttributeError as e:
            if remove_disordered:
                raise(e)
            self.save_structure(remove_disordered=True)


    def save_sequences(self):
        self.chain_numbering_extended_dict = {}
        self.chain_sequence_dict = {}
        for chain_id in self.chain_ids:
            chain_sequence, chain_numbering_extended = self.get_chain_sequence_and_numbering(chain_id)
            self.chain_numbering_extended_dict[chain_id] = chain_numbering_extended
            self.chain_sequence_dict[chain_id] = chain_sequence
            with open(self.output_path + self.pdb_id + chain_id + '.seq.txt', 'w') as f:
                f.write('>' + self.pdb_id + chain_id + '\n')
                f.write(chain_sequence + '\n')
                f.write('\n')




#%%
if __name__ == '__main__':
    pdbPath = '/home/niklas/pdb_database/structures/divided/pdb/'
    pdbCode = '1VOK'
    chains = ['A', 'B']
    domainBoundaries = [[23, 115], [29, 198]]
    output_path = './'

    import logging
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
#    handler = logging.FileHandler(tmp_path + 'templates.log', mode='w', delay=True)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)


