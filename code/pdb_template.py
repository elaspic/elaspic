# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 18:58:50 2012

@author: niklas
"""
import numpy as np
import gzip
import subprocess
import string

import Bio
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import errors

from collections import defaultdict

A_DICT = {
    'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', 
    'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 
    'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 
    'Y':'TYR', 'V':'VAL', 'U':'SEC', 'O':'PYL', 
    'B':'ASX', 'Z':'GLX', 'J':'XLE', 'X':'XAA', '*':'TER'
}
AAA_DICT = dict([(value,key) for key,value in A_DICT.items()])
AAA_DICT['UNK'] = 'X'
AAA_DICT['MSE'] = 'M'
amino_acids = AAA_DICT.keys()
methylated_lysines = ['MLZ', 'MLY', 'M3L']
lysine_atoms = ['N', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'C', 'O']


def get_pdb(pdb_code, pdb_path, tmp_path='/tmp/'):
    """
    parse a pdb file with biopythons PDBParser() and return the structure

    input: pdb_code  type String     four letter code of the PDB file

    return: Biopython pdb structure
    """
    if False:
        pdb_path = pdb_path + '../../../biounit/coordinates/divided/'
        parser = PDBParser(QUIET=True)
        prefix = ''
        suffix = '.pdb1.gz'
    if False:
        pdb_path = pdb_path + '../mmCIF/'
        parser = MMCIFParser()
        prefix = ''
        suffix = '.cif.gz'
    if True:
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        prefix = 'pdb'
        suffix = '.ent.gz'
    pdb_file = pdb_path + pdb_code[1:3].lower() + '/' + prefix + pdb_code.lower() + suffix
    try:
        with gzip.open(pdb_file, 'r') as ifh, \
        open(tmp_path + pdb_code + suffix.replace('.gz',''), 'w') as ofh:
            ofh.write(ifh.read())
    except IOError as e:
        raise e
        raise errors.NoPDBFound(pdb_file)
#    # Convert methylated lysines to lysines directly in the pdb file
#    with open(tmp_path + pdb_code + suffix.replace('.gz',''), 'r') as ifh:
#        pdb_file = ifh.readlines()
#    line_idx = 0
#    while line_idx < len(pdb_file):
#        increment = True
#        for resname in methylated_lysines:
#            if ' '+resname+' ' not in pdb_file[line_idx]:
#                continue # check the next residue
#            elif (pdb_file[line_idx][:6] == 'HETATM' 
#                    and pdb_file[line_idx][17:20] == resname
#                    and pdb_file[line_idx][13:16] in ['N  ', 'CA ', 'CB ', 'CG ', 'CD ', 'CE ', 'NZ ', 'C  ', 'O  ']):
#                pdb_file[line_idx] = pdb_file[line_idx].replace('HETATM', 'ATOM  ').replace(resname, 'LYS')
#            elif pdb_file[line_idx][:6] == 'SEQRES':
#                pdb_file[line_idx] = pdb_file[line_idx].replace(resname, 'LYS')
#            else:
#                print 'Removing line: {}'.format(pdb_file[line_idx])
#                del pdb_file[line_idx]
#                increment = False
#                break
#        if increment:
#            line_idx += 1
#    with open(tmp_path + pdb_code + suffix.replace('.gz','.corrected'), 'w') as ofh:
#        ofh.writelines(pdb_file)
#    result = parser.get_structure('ID', tmp_path + pdb_code + suffix.replace('.gz','.corrected'))
    result = parser.get_structure('ID', tmp_path + pdb_code + suffix.replace('.gz',''))
    return result


def convert_aa(aa):
    """
    Convert amino acids from three letter code to one letter code or vice versa
    """
    if len(aa) == 3:
        try:
            return AAA_DICT[aa.upper()]
        except KeyError:
            print  'Not a valid amino acid'
            return
    if len(aa) == 1:
        try:
            return A_DICT[aa.upper()]
        except KeyError:
            print  'Not a valid amino acid'
            return
    print 'Not a valid amino acid'


def convert_position_to_resid(model, chain_id, positions):
    """ 
    Maps the mutation sequence position of the pdb to pdb numbering

    input
    model       class 'Bio.PDB.Model.Model'
    mutation    type 'str'                      ; B_Q61L

    return:
    chainNumbering[position-1]      type 'int'
    """

    chain = model[chain_id]
    chainNumbering = get_chain_numbering(chain, return_extended=True)
    return [chainNumbering[p-1] for p in positions]


def convert_resid_to_position(model, pdb_chain, resids, pdb_domain_start, pdb_domain_end):
    """ 
    Maps the mutation sequence position of the pdb to pdb numbering

    input
    model       class 'Bio.PDB.Model.Model'
    mutation    type 'str'                      ; B_Q61L

    return:
    chainNumbering[position-1]      type 'int'
    """

    chain = model[pdb_chain]
    chain_numbering = get_chain_numbering(chain, return_extended=True)
#    try:
#        chain_numbering = chain_numbering[chain_numbering.index(pdb_domain_start):chain_numbering.index(pdb_domain_end)+1]
#    except ValueError:
#        raise errors.PDBDomainDefsError(
#        'pdb domain start %s or pdb domain end %s not found in model %s and chain %s, with chain numbering: %s' % \
#        (pdb_domain_start, pdb_domain_end, model.id, chain.id, ','.join(chainNumbering)))
    return [chain_numbering.index(resid) for resid in resids if resid in chain_numbering]


def get_chain_numbering(chain, return_sequence=False, return_extended=False, include_hetatms=False):
    """
    Returns a list with the numbering of the chains
    """
    chain_numbering = []
    chain_numbering_extended = []
    chain_sequence = []
    for res in chain:
        if include_hetatms or res.resname in amino_acids:
            chain_numbering.append(res.id[1])
            chain_numbering_extended.append(str(res.id[1]) + res.id[2].strip())
            chain_sequence.append(AAA_DICT[res.resname])

    chain_sequence = ''.join(chain_sequence)
#        chainNumbering = [residue.id[1] for residue in chain if is_aa(residue, standard=True)]
    if return_sequence and return_extended:
        return chain_numbering_extended, chain_sequence
    elif return_sequence:
        return chain_numbering, chain_sequence
    elif return_extended:
        return chain_numbering_extended
    else:
        return chain_numbering


def get_seqres_sequence(chain):
    """
    Extracts a sequence from a PDB file.
    Usefull when interested in the sequence that was used for crystallization
    and not the ATOM sequence
    """
    # setting aa_only to False Selenomethionines are reported in the
    # sequence as well
    # see: http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html
    sequence = Seq('', IUPAC.protein)
    for pb in PPBuilder().build_peptides(chain, aa_only=False):
#        for pb in PPBuilder().build_peptides(chain, aa_only=True):
        tmp = sequence + pb.get_sequence()
        sequence = tmp
    return sequence


def get_chain_sequences(file_or_structure, seqres_sequence=False):
    """
    """
    if isinstance(file_or_structure, basestring):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('ID', file_or_structure)
        model = structure[0]
    elif isinstance(file_or_structure, Bio.PDB.Structure.Structure):
        model = file_or_structure[0]
    elif isinstance(file_or_structure, Bio.PDB.Model.Model):
        model = file_or_structure
    elif isinstance(file_or_structure, Bio.PDB.Chain.Chain):
        model = [file_or_structure]
        
    chain_sequences = defaultdict(list)
    for chain in model:
        if seqres_sequence:
            chain_sequence = get_seqres_sequence(chain)
        else:
            __, chain_sequence = get_chain_numbering(chain, return_sequence=True)
        chain_sequences[chain.id] = chain_sequence
    return chain_sequences

    
def convert_resnum_alphanumeric_to_numeric(resnum):
    idx_increment = 0
    while string.letters.find(resnum[-1]) != -1:
        idx_increment += string.letters.find(resnum[-1])
        resnum = resnum[:-1]
    resnum = int(resnum) + idx_increment
    return resnum


class PDBTemplate():

    def __init__(self, pdbPath, pdb_id, chain_ids, domain_boundaries, outputPath, tmp_path, log):
        """
        Fetches the pdb file from the local database, extracts the desired chains
        and exports the amino acid sequence of the chains

        input:  pdbPath:    type string     path of the locally installed pdb database
                pdbCode:    type string     four letter code of the pdb structure
                chains:     type list       list containing the chain names
                                            as strings, i.e. chains = ['A', 'B']
                domainBoundaries
                outputPath  type string     path where the extracted pdb structure
                                            and sequence are saved
        """
        self.pdbPath = pdbPath
        self.pdb_id = pdb_id
        self.outputPath = outputPath
        self.structure = get_pdb(self.pdb_id, self.pdbPath, tmp_path)
        if chain_ids:
            self.chain_ids = chain_ids
        else: # If [], extract all chains
            self.chain_ids = [chain.id for chain in self.structure[0].child_list]
        self.domain_boundaries = domain_boundaries # If [], extract the entire domain
        self.log = log


    def extract(self):
        """
        Extract the wanted chains out of the PDB file
        Removes water atoms and selects the domain regions
        """
        io = PDBIO()
        model = self.structure[0] # assuming that model 0 is always the desired one
        new_structure = Bio.PDB.Structure.Structure('ID')
        new_model = Bio.PDB.Model.Model(0)
        chain_for_hetatms = Bio.PDB.Chain.Chain('Z')
        for i, chain_id in enumerate(self.chain_ids):
            chain = model[chain_id]
            BEGIN = False
            END = False
            res_idx = 0
            while res_idx < len(chain):
                res = chain.child_list[res_idx]
                old_res_id = res.id
                # Remove water
                if res.id[0] == 'W':
#                    self.log.debug('Removing water molecule {}'.format(res.id))
                    chain.detach_child(res.id)
                    continue
                # Convert methylated lysines to regular lysines
                if res.resname in methylated_lysines:
                    new_resname = 'LYS'
                    new_resid = (' ', res.id[1], res.id[2])
                    self.log.debug(
                        'Renaming residue {} {} to {} {}'
                        .format(res.resname, res.id, new_resname, new_resid))
                    res.resname = new_resname
                    res.id = new_resid
                    atom_idx = 0
                    while atom_idx < len(res):
                        atom_id = res.child_list[atom_idx].id
                        if atom_id not in lysine_atoms:
                            self.log.debug('Removing atom {} from residue {} {}.'.format(atom_id, res.resname, res.id))
                            res.detach_child(atom_id)
                        else:
                            atom_idx += 1
                # Move hetatms to the hetatm chain
                if res.resname not in amino_acids:
                    self.log.debug('Moving hetatm residue {} {} to the hetatm chain'.format(res.resname, res.id))
                    chain.detach_child(res.id)
                    hetatm_res = res
                    hetatm_res.id = (hetatm_res.id[0], len(chain_for_hetatms)+1, hetatm_res.id[2],)
                    chain_for_hetatms.add(hetatm_res)
                    continue
                # Cut each chain to domain boundaries
                if self.domain_boundaries and self.domain_boundaries[i]:                   
                    resnum_extended = convert_resnum_alphanumeric_to_numeric(str(res.id[1]) + res.id[2].strip())
                    domain_start_extended = convert_resnum_alphanumeric_to_numeric(self.domain_boundaries[i][0])
                    domain_end_extended = convert_resnum_alphanumeric_to_numeric(self.domain_boundaries[i][1]) + 1
                    if resnum_extended >= domain_start_extended:
                        BEGIN = True # From now on, keep the residues
                    if resnum_extended >= domain_end_extended:
                        END = True # From now on, detach them again
                    if not BEGIN or END: # We are before the beginning or past the end
                        chain.detach_child(old_res_id)
                        continue
                res_idx += 1
            new_model.add(chain)
        # Remove hetatms if they are > 10A away from the chains of interest
        chain_for_hetatms.id = [c for c in reversed(string.uppercase) if c not in self.chain_ids][0]
        for res_1 in chain_for_hetatms:
            in_contact = False
            for chain_2 in new_model:
                for res_2 in chain_2:
                    for atom_1 in res_1:
                        for atom_2 in res_2:
                            a = atom_1.coord
                            b = atom_2.coord
                            r = np.sqrt(sum( (a - b)**2 for a, b in zip(a, b)))
                            if r < 10:
                                in_contact = True
            if not in_contact:
                chain_for_hetatms.detach_child(res_1.id)
        if chain_for_hetatms:
            self.log.debug('Adding hetatm chain of length {}'.format(len(chain_for_hetatms)))
            new_model.add(chain_for_hetatms)
        # If the hetatm chain is not empty, add it to the model
        new_structure.add(new_model)
        self.structure = new_structure
        # Save the structure to a pdb file
        io.set_structure(self.structure)
        outFile = self.outputPath + self.pdb_id + ''.join(self.chain_ids) + '.pdb'
        io.save(outFile)
        # Save the sequence of each chain to a fasta file
        for chain_id in self.chain_ids:
            chain_numbering_extended, chain_sequence = self.get_chain_numbering(chain_id, return_sequence=True, return_extended=True)        
            with open(self.outputPath + self.pdb_id + chain_id + '.seq.txt', 'w') as f:
                f.write('>' + self.pdb_id + chain_id + '\n')
                f.write(chain_sequence + '\n')
                f.write('\n')


    def get_chain_numbering(self, chain_id, return_sequence=False, return_extended=False, include_hetatms=False):
        """
        """
        chain = self.structure[0][chain_id]
        return get_chain_numbering(chain, return_sequence=return_sequence, return_extended=return_extended, include_hetatms=include_hetatms)


    def get_seqres_sequence(self, chain_id):
        """
        Extracts a sequence from a PDB file.
        Usefull when interested in the sequence that was used for crystallization
        and not the ATOM sequence
        """
        if not chain_id.strip():
            raise Exception('chain_id must not be empty!')
        if len(chain_id) > 1:
            raise Exception('chain_id must contain only a single chain!')
        if isinstance(chain_id, list):
            chain_id = chain_id[0]
        chain = self.pdbStructure[0][chain_id]
        sequence = get_seqres_sequence(chain)
        return sequence


if __name__ == '__main__':
    pdbPath = '/home/niklas/pdb_database/structures/divided/pdb/'
    pdbCode = '1VOK'
    chains = ['A', 'B']
    domainBoundaries = [[23, 115], [29, 198]]
    outputPath = './'

    import logging
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
#    handler = logging.FileHandler(tmp_path + 'templates.log', mode='w', delay=True)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    log = logger

    p = PDBTemplate(pdbPath, pdbCode, chains, domainBoundaries, outputPath, log)
    print p.extract()
