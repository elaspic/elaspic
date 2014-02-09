# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:29:13 2013

@author: niklas
"""

import subprocess
import gzip 
from math import sqrt

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import DSSP

import class_error as error



class AnalyzeStructure():
    """
    Runs the program pops to calculate the interface size of the complexes
    This is done by calculating the surface of the complex and the seperated parts.
    The interface is then given by the substracting
    """
    
    def __init__(self, tmpPath, pdb_file, chains, log, domain_defs):

        self.working_path = tmpPath # modeller_path, foldx_path
        
        self.log = log
        self.pdb_file = pdb_file
        self.chain_ids = chains
        self.domain_defs = domain_defs
        
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('ID', self.working_path + self.pdb_file)
        
        self.__split_pdb_into_chains()
        
        self.amino_acids = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN',
               'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 
               'MET', 'PHE', 'TYR', 'TRP']
        
    ###########################################################################
    
    def get_sasa(self):
        
        # Use pops to calculate sasa score for each chain
        for chain_id in self.chain_ids:
            termination, rc, e = self.__run_pops_aa(chain_id)
            if termination != 'Clean termination':
                self.log.error('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
                self.log.error(e)
                raise error.PopsError(e, self.working_path + self.pdb_file, self.chain_ids)
            else:
                self.log.warning('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
                self.log.warning(e)
        
        # Read the sasa scores from text files into a dictionary
        sasa_score = {}
        for chain_id in self.chain_ids:
            sasa_score[chain_id] = self.__read_pops_aa(self.working_path + chain_id + '.out')
        
        # Confirm that the dictionary has the right chains and output
        assert set(sasa_score.keys()) == set(self.chain_ids)
        return sasa_score


    
    def get_interface_area(self):
        
        termination, rc, e = self.__run_pops_area(self.working_path + self.pdb_file)
        if rc != 0:
            if termination != 'Clean termination':
                self.log.error('Pops error for pdb: %s:' % self.pdb_file)
                self.log.error(e)
                return '0', '0', '0'
            else:
                self.log.warning('Pops warning for pdb: %s:' % self.pdb_file)
                self.log.warning(e)
        result = self.__read_pops_area(self.working_path + self.pdb_file.replace('pdb', 'out'))
        
        # Distinguish the surface area by hydrophobic, hydrophilic, and total
        for item in result:
            if item[0] == 'hydrophobic:':
                hydrophobic = float(item[1])
            elif item[0] == 'hydrophilic:':
                hydrophilic = float(item[1])
            elif item[0] == 'total:':
                total = float(item[1])
        sasa_complex = hydrophobic, hydrophilic, total
        
        # calculate SASA for chain, i.e. part one of the complex:
        termination, rc, e = self.__run_pops_area(self.working_path + self.chain_ids[0] + '.pdb')
        if rc != 0:
            if termination != 'Clean termination':
                self.log.error('Error in pops for pdb: %s:' % self.pdb_file)
                return '0', '0', '0'
            else:
                self.log.warning('Warning in pops for pdb: %s:' % self.pdb_file)
        result = self.__read_pops_area(self.working_path + self.chain_ids[0] + '.out')
        
        # Distinguish the surface area by hydrophobic, hydrophilic, and total
        for item in result:
            if item[0] == 'hydrophobic:':
                hydrophobic = float(item[1])
            elif item[0] == 'hydrophilic:':
                hydrophilic = float(item[1])
            elif item[0] == 'total:':
                total = float(item[1])
        sasa_chain = hydrophobic, hydrophilic, total
        
        # calculate SASA for oppositeChain, i.e. the second part of the complex:
        termination, rc, e = self.__run_pops_area(self.working_path + self.chain_ids[1] + '.pdb')
        if rc != 0:
            if termination != 'Clean termination':
                self.log.error('Error in pops for pdb: %s:' % self.pdb_file)
                self.log.error(e)
                return '0', '0', '0'
            else:
                self.log.error('Warning in pops for pdb: %s:' % self.pdb_file)
                self.log.error(e)
        result = self.__read_pops_area(self.working_path + self.chain_ids[1] + '.out')
        
        for item in result:
            if item[0] == 'hydrophobic:':
                hydrophobic = float(item[1])
            elif item[0] == 'hydrophilic:':
                hydrophilic = float(item[1])
            elif item[0] == 'total:':
                total = float(item[1])
        sasa_oppositeChain = hydrophobic, hydrophilic, total
        
        sasa = [ 0, 0, 0 ]
        # hydrophobic
        sasa[0] = (sasa_chain[0] + sasa_oppositeChain[0] - sasa_complex[0]) / 2.0
        # hydrophilic
        sasa[1] = (sasa_chain[1] + sasa_oppositeChain[1] - sasa_complex[1]) / 2.0
        # total
        sasa[2] = (sasa_chain[2] + sasa_oppositeChain[2] - sasa_complex[2]) / 2.0
        
        return str(sasa[0]), str(sasa[1]), str(sasa[2])



    def __split_pdb_into_chains(self):
        
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        io = PDBIO()
#        if len(self.chain_ids) == 1:
#            # If there's only one pdb, it might not have it's chain numbered
#            structure = parser.get_structure('ID', self.working_path + self.pdb_file)
#            model = self.structure[0]
#            for child in model.get_list():
#                
#            io.set_structure(structure)
#            outFile = self.working_path + self.chain_ids[0] + '.pdb'
#            io.save(outFile)
#        else:
        for chain_id in self.chain_ids:
            # save chain, i.e. part one of the complex:
            self.log.debug(self.working_path + self.pdb_file)
            structure = parser.get_structure('ID', self.working_path + self.pdb_file)
            model = structure[0]
            for child in model.get_list():
                self.log.debug('child id:' + child.id)
                if child.id == '' or child.id == ' ':
                    child.id = chain_id
                    for c in child:
                        c.id = (c.id[0], c.id[1]+100, c.id[2],)
                if child.id != chain_id:
                    model.detach_child(child.id)
            io.set_structure(structure)
            outFile = self.working_path + chain_id + '.pdb'
            io.save(outFile)
            

        
    def __run_pops_aa(self, chain_id):
        system_command = (self.working_path + 'pops --pdb ' + self.working_path + chain_id + '.pdb' + 
                        ' --noHeaderOut --noTotalOut --atomOut --popsOut ' + self.working_path + chain_id + '.out')
        childProcess = subprocess.Popen(system_command, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True)
        result, e = childProcess.communicate()
        rc = childProcess.returncode
        
        # the returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        self.log.debug('result: %.30s' % result.replace('\n','; '))
        output = [ line for line in result.split('\n') if line != '' ]
        return output[-1], rc, e

            

    def __read_pops_aa(self, filename):
        """
        Read pops sasa results atom by atom, ignoring all main chain atoms except for Ca
        """
        # The new way
        ignore = ['N', 'C', 'O']
        per_residue_sasa_scores = []
        current_residue_number = None
        with open(filename, 'r') as fh:
            for line in fh:
                row = line.split()
                if len(row) != 11:
                    continue
                atom_number, atom_name, residue_name, chain, residue_number, sasa, __, __, __, __, sa = line.split()
                atom_number, residue_number, sasa, sa = int(atom_number), int(residue_number), float(sasa), float(sa)
                if atom_name in ignore:
                    continue
                if current_residue_number != residue_number:
                    if current_residue_number:
                        per_residue_sasa_scores.append(total_sasa/total_sa)
                    current_residue_number = residue_number
                    total_sasa = 0
                    total_sa = 0
                total_sasa += sasa
                total_sa += sa
            per_residue_sasa_scores.append(total_sasa/total_sa)
        return per_residue_sasa_scores



    def __run_pops_area(self, full_filename):
        system_command = self.working_path + 'pops --pdb ' + full_filename + ' --chainOut --popsOut ' + full_filename.replace('pdb', 'out')
        childProcess = subprocess.Popen(system_command, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True, 
                                        )
        result, e = childProcess.communicate()
        rc = childProcess.returncode
        self.log.debug('result: %.30s' % result.replace('\n','; '))
        # the returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        output = [ line for line in result.split('\n') if line != '' ]
        return output[-1], rc, e



    def __read_pops_area(self, filename):
        # The old way
        keep = ['hydrophobic:', 'hydrophilic:', 'total:']
        with open(filename, 'r') as pops:
            result = [ x.split(' ') for x in pops.readlines() if x != '' and x.split(' ')[0] in keep ]
        return [ [ x.strip() for x in item if x != '' ] for item in result ]


###############################################################################

    def get_interacting_aa(self):
        
        model = self.structure[0]
        chains = [ chain for chain in model]
        
        interacting_aa = {}
        for chain_1 in chains:
            interacting_aa[chain_1.id] = set()
            for idx, residue_1 in enumerate(chain_1):
                
                for chain_2 in [c for c in chains if c != chain_1]:
                    for residue_2 in chain_2:
                        
                        if residue_1.resname not in self.amino_acids \
                        or residue_2.resname not in self.amino_acids:
                            continue
                
                        for atom_1 in residue_1:
                            for atom_2 in residue_2:
                                r = self.calculate_distance(atom_1, atom_2)
                                if r <= 5.0:
                                    interacting_aa[chain_1.id].add(idx+1) # pdb domain defs and indices always start from 1
            # Change set to list
            interacting_aa[chain_1.id] = list(interacting_aa[chain_1.id])
        
        self.log.debug('interacting_aa_keys:')
        self.log.debug(interacting_aa.keys())
        self.log.debug('chain ids: %s' % self.chain_ids)
        return interacting_aa
        
     
    def calculate_distance(self, atom1, atom2):
        """
        returns the distance of two points in three dimensional space
        
        input: atom instance of biopython: class 'Bio.PDB.Atom.Atom
        
        return: type 'float'
        """
        a = atom1.coord
        b = atom2.coord
        assert(len(a) == 3 and len(b) == 3)
        return sqrt(sum( (a - b)**2 for a, b in zip(a, b)))



###############################################################################
   
    def get_dssp(self, residue, position, dssp_bin_name='dssp-2.0.4-linux-amd64'):
        """
        Make sure position is a Biophython position entry like
        <Residue ASP het=  resseq=15 icode= >, i.e. (' ', 15, ' ')
        """
        if len(residue) == 1:
            residue = self.convert_aa(residue)
        
        model = self.structure[0]
        dssp = DSSP(model, self.working_path + self.pdb_file, self.working_path + dssp_bin_name)
        
        result = dssp[self.chain, (' ', int(position), ' ')]
        
        # make sure that one is looking at the correct amino acid
        assert( result[0].resname == residue.upper() )
        
        # result[1] == secondary structure information
        # result[2] == solvent accessibilty
        return result[1], result[2]
    
    
    def convert_aa(self, aa):
        """
        convert amino acids from three letter code to one letter code or vice versa
        """
        A_DICT = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', \
                  'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', \
                  'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', \
                  'Y':'TYR', 'V':'VAL', 'U':'SEC', 'O':'PYL', \
                  'B':'ASX', 'Z':'GLX', 'J':'XLE', 'X':'XAA', '*':'TER'}
        
        AAA_DICT = dict([(value,key) for key,value in A_DICT.items()])
        
        if len(aa) == 3:
            try:
                return AAA_DICT[aa.upper()]
            except KeyError:
                self.log.error('Not a valid amino acid')
                return
        if len(aa) == 1:
            try:
                return A_DICT[aa.upper()]
            except KeyError:
                self.log.error('Not a valid amino acid')
                return
        self.log.error('Not a valid amino acid')
        
        
        

###############################################################################        
# Old code used to deal with single-point mutations
# We now calculate dssp and interface amino acids for the entire model
    def __check_structure(self, pdbCode, chainID, mutation):
        """ checks if the mutation falls into the interface, i.e. is in contact with
        another chain
        'mutation' has to be of the form A_T70H, mutation in chain A, from Tyr at
        position 70 to His
        NOTE: takes the mutation as numbered ins sequence! The conversion is done
        within this function!
        
        input
        pdbCode     type 'str'
        chainID     type 'str'
        mutation    type 'str'      ; B_Q61L
        
        return:
        contacts    type 'dict'     ; {'C': False, 'B': True}
                                      key:   chainID                type 'str'
                                      value: contact to chainID     type boolean
        """
        structure = self.getPDB(pdbCode, self.pdbPath)
        model = structure[0]
        
        chains   = [ chain for chain in model]
        chainIDs = [ chain.id for chain in model]
    
        position = self.convert_mutation_position(model, mutation) # convert the position numbering
        contacts = { chainID: False for chainID in chainIDs if not chainID == mutation[0] }
        
        for i in range(len(chains)):
            if chains[i].id == mutation[0]:
                chain = chains[i]
                # use list expansion to select only the 'opposing chains'
                oppositeChains = [ x for x in chains if x != chains[i] ]
        
        # If the residues do not match, issue a warning.
        # To obtain a better model one could restrict to templates that have
        # the same amino acid as the uniprot sequence at the position of the mutation.
#        if chain[position].resname != self.convert_aa(fromAA):
#            print 'Residue missmatch while checking the structure!'
#            print 'pdbCode', pdbCode
#            print 'mutation', mutation
#            print chain[position].resname, self.convert_aa(fromAA)
#            print 'position', position

       
        for oppositeChain in oppositeChains:
            # check each residue
            for residue in oppositeChain:
                # for each residue each atom of the mutated residue has to be checked
                for atom1 in chain[position]: # chain[position] is the residue that should be mutated
                    # and each atom
                    for atom2 in residue:
                        r = self.distance(atom1, atom2)
                        if r <= 5.0:
                            contacts[oppositeChain.id] = True

        return contacts



    def convert_mutation_position(self, model, mutation):
        """ maps the mutation sequence position of the pdb to pdb numbering
        
        input
        model       class 'Bio.PDB.Model.Model'
        mutation    type 'str'                      ; B_Q61L
        
        return:
        chainNumbering[position-1]      type 'int'
        """
        chain = model[mutation[0]]
        position = int(mutation[3:-1])
        
        chainNumbering = self.getChainNumberingNOHETATMS(chain)

        return chainNumbering[position-1]



    def getChainNumberingNOHETATMS(self, chain):
        """
        returns a list with the numbering of the chains
        
        input:
        chain               class 'Bio.PDB.Chain.Chain'
        
        return:
        chainNumbering      type 'list' of 'int'
        """
        chainNumbering = list()
        for residue in chain:
            if residue.resname in self.amino_acids and residue.id[0] == ' ':
                chainNumbering.append(residue.id[1])

#        chainNumbering = [residue.id[1] for residue in chain if is_aa(residue, standard=True)]
        return chainNumbering    
 

   
    def getPDB(self, pdbCode, pdbPath):
        """
        parse a pdb file with biopythons PDBParser() and return the structure
        
        input: pdbCode  type String     four letter code of the PDB file
        
        return: Biopython pdb structure
        
        input:
        pdbCode     type 'str'
        pdbPath     type 'str'
        
        return:
        result      type class 'Bio.PDB.Structure.Structure'
        """
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        pdbFile = pdbPath + pdbCode[1:3].lower() + '/pdb' + pdbCode.lower() + '.ent.gz'
        pdbFileUncompressed = gzip.open(pdbFile, 'r')
        result = parser.get_structure('ID', pdbFileUncompressed)

        return result

    
###############################################################################
    
if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    analyze_structure = AnalyzeStructure('/tmp/cosmic_driver_proteins/Consumer-1/modeller/', 'Q9UKG1.BL00030001.pdb', ['A',], logger)
    sasa_score = analyze_structure.get_sasa()
#    interacting_aa = analyze_structure.get_interacting_aa()
#    interface_area = analyze_structure.get_interface_area()
    
    print sasa_score
#    print interacting_aa
#    print interface_area
    
    
    
    
    