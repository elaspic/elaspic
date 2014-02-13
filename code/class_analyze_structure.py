# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:29:13 2013

@author: niklas
"""
import os
import subprocess
import gzip 
from math import sqrt

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import DSSP

import logging
import class_error as error



class AnalyzeStructure(object):
    """
    Runs the program pops to calculate the interface size of the complexes
    This is done by calculating the surface of the complex and the seperated parts.
    The interface is then given by the substracting
    """
    
    def __init__(self, data_path, working_path, pdb_file, chains, domain_defs, logger):

        self.data_path = data_path # modeller_path, foldx_path        
        self.working_path = working_path 
        
        self.pdb_file = pdb_file
        self.chain_ids = chains
        self.domain_defs = domain_defs
        
        #
        if logger:
            self.log = logger
        else:
            logger = logging.getLogger(__name__)
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            self.log = logger
        
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('ID', self.data_path + self.pdb_file)
        
        self.__split_pdb_into_chains()
        
        self.amino_acids = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN',
               'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 
               'MET', 'PHE', 'TYR', 'TRP']
               
    
    def __split_pdb_into_chains(self):
        
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        io = PDBIO()
        for chain_id in self.chain_ids:
            # save chain, i.e. part one of the complex:
            self.log.debug(self.data_path + self.pdb_file)
            structure = parser.get_structure('ID', self.data_path + self.pdb_file)
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
        
        # All chains together
        structure = parser.get_structure('ID', self.data_path + self.pdb_file)
        model = structure[0]        
        io.set_structure(structure)
        outFile = self.working_path +  self.pdb_file
        io.save(outFile)
    
    
    ###########################################################################
    
    def get_sasa(self, program_to_use='naccess'):
        
        if program_to_use == 'naccess':
            run_sasa_atom = self._run_naccess_atom
        elif program_to_use == 'pops':
            run_sasa_atom = self._run_pops_atom
        else:
            raise Exception('Unknown program specified!')
            
        sasa_score_splitchains = {}
        for chain_id in self.chain_ids:
            sasa_score_splitchains.update(run_sasa_atom(chain_id + '.pdb'))
        sasa_score_allchains = run_sasa_atom(self.pdb_file)
        return [sasa_score_splitchains, sasa_score_allchains]


    def _run_naccess_atom(self, filename):
        # run naccess
        current_path = os.getcwd()
        os.chdir(self.working_path)
        system_command = ('./naccess ' + filename)
        self.log.debug('naccess system command: %s' % system_command)
        childProcess = subprocess.Popen(system_command, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True)
        result, e = childProcess.communicate()
        rc = childProcess.returncode
        
        # collect results
        sasa_scores = {}
        with open(self.working_path + filename.split('.')[0] + '.rsa') as fh:
            for line in fh:
                row = line.split()
                if row[0] != 'RES':
                    continue
                (line_id, res, chain, num, all_abs, all_rel, 
                 sidechain_abs, sidechain_rel, mainchain_abs, mainchain_rel, 
                 nonpolar_abs, nonpolar_rel, polar_abs, polar_rel) = row
                sasa_scores.setdefault(chain, []).append(sidechain_rel) # percent sasa on sidechain
                
        os.chdir(current_path)
        return sasa_scores
               
            
    def _run_pops_atom(self, chain_id):
        # Use pops to calculate sasa score for the given chain
        termination, rc, e = self.__run_pops_atom(chain_id)
        if termination != 'Clean termination':
            self.log.error('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
            self.log.error(e)
            raise error.PopsError(e, self.data_path + self.pdb_file, self.chain_ids)
        else:
            self.log.warning('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
            self.log.warning(e)
        
        # Read the sasa scores from a text file
        sasa_scores = self.__read_pops_atom(chain_id)
        return sasa_scores
        
        
    def __run_pops_atom(self, chain_id):
        system_command = (self.data_path + 'pops --pdb ' + self.working_path + chain_id + '.pdb' + 
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
            

    def __read_pops_atom(self, chain_id):
        """
        Read pops sasa results atom by atom, ignoring all main chain atoms except for Ca
        """
        # The new way
        ignore = ['N', 'C', 'O']
        per_residue_sasa_scores = []
        current_residue_number = None
        with open(self.working_path + chain_id + '.out', 'r') as fh:
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


###############################################################################

    def get_dssp(self):
        """
        """
        
        current_path = os.getcwd()
        os.chdir(self.working_path)
        system_command = ('./dssp -i ' + self.pdb_file + ' -o ' + 'dssp_results.txt')
        self.log.debug('dssp system command: %s' % system_command)
        childProcess = subprocess.Popen(system_command, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True)
        result, e = childProcess.communicate()
        rc = childProcess.returncode
        self.log.debug('dssp return code: %i' % rc)
        self.log.debug('dssp result: %s' % result)
        self.log.debug('dssp result: %s' % e)
        
        # collect results
        dssp_ss = {}
        dssp_acc = {}
        start = False
        with open(self.working_path + 'dssp_results.txt') as fh:
            for l in fh:
                row = l.split()
                if not row or len(row) < 2:
                    continue
                if row[1] == "RESIDUE":
                    # Start parsing from here
                    start = True
                    continue
                if not start:
                    continue
                if l[9] == ' ':
                    # Skip -- missing residue
                    continue
                resseq, icode, chainid, aa, ss = int(l[5:10]), l[10], l[11], l[13], l[16]
                if ss == ' ':
                    ss = '-'
                try:
                    acc = int(l[34:38])
                    phi = float(l[103:109])
                    psi = float(l[109:115])
                except ValueError, exc:
                    # DSSP output breaks its own format when there are >9999
                    # residues, since only 4 digits are allocated to the seq num
                    # field.  See 3kic chain T res 321, 1vsy chain T res 6077.
                    # Here, look for whitespace to figure out the number of extra
                    # digits, and shift parsing the rest of the line by that amount.
                    if l[34] != ' ':
                        shift = l[34:].find(' ')
                        acc = int((l[34+shift:38+shift]))
                        phi = float(l[103+shift:109+shift])
                        psi = float(l[109+shift:115+shift])
                    else:
                        raise ValueError(exc)
                dssp_ss.setdefault(chainid, []).append(ss) # percent sasa on sidechain
                dssp_acc.setdefault(chainid, []).append(acc)
        for key in dssp_ss.keys():
            dssp_ss[key] = ''.join(dssp_ss[key])
        os.chdir(current_path)
        return dssp_ss, dssp_acc


###############################################################################

    def get_interchain_distances(self, chain_mutation_structure=None):
        """
        """
        
        model = self.structure[0]
        chains = [ chain for chain in model ]
        
        if chain_mutation_structure:
            position = self.convert_mutation_position(model, chain_mutation_structure) # convert the position numbering
        
        shortest_interchain_distances = {}
        for chain_1 in chains:
            if chain_mutation_structure and chain_1.id != chain_mutation_structure[0]:
                continue # skip chains that we are not interested in
            shortest_interchain_distances[chain_1.id] = list()
            for idx, residue_1 in enumerate(chain_1):
                if chain_mutation_structure and residue_1.id[1] != position:
                    continue # skip all residues that we are not interested in
                min_r = None
                for chain_2 in [c for c in chains if c != chain_1]:
                    for residue_2 in chain_2:
                        if residue_1.resname not in self.amino_acids \
                        or residue_2.resname not in self.amino_acids:
                            continue
                        
                        for atom_1 in residue_1:
                            for atom_2 in residue_2:
                                r = self.calculate_distance(atom_1, atom_2)
                                if not min_r or min_r > r:
                                    min_r = r
                shortest_interchain_distances[chain_1.id].append(min_r)
        
        self.log.debug('interacting_aa_keys:')
        self.log.debug(shortest_interchain_distances.keys())
        self.log.debug('chain ids: %s' % self.chain_ids)
        return shortest_interchain_distances
        
     
    def calculate_distance(self, atom1, atom2):
        """
        returns the distance of two points in three dimensional space
        input: atom instance of biopython: class 'Bio.PDB.Atom.Atom
        return: type 'float'
        """
        
        a = atom1.coord
        b = atom2.coord
        assert(len(a) == 3 and len(b) == 3)
        return sqrt(sum( (a - b)**2 for a, b in zip(a, b) ))


    def convert_mutation_position(self, model, chain_mutation_structure):
        """ maps the mutation sequence position of the pdb to pdb numbering
        
        input
        model       class 'Bio.PDB.Model.Model'
        mutation    type 'str'                      ; B_Q61L
        
        return:
        chainNumbering[position-1]      type 'int'
        """
        
        chain = model[chain_mutation_structure[0]]
        position = int(chain_mutation_structure[3:-1])
        
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
        
      
###############################################################################
      
    def get_interface_area(self):
        
        termination, rc, e = self.__run_pops_area(self.data_path + self.pdb_file)
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
        result = self.__read_pops_area(self.data_path + self.chain_ids[0] + '.out')
        
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
        result = self.__read_pops_area(self.data_path + self.chain_ids[1] + '.out')
        
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


    def __run_pops_area(self, full_filename):
        system_command = (self.working_path + ' pops --chainOut'
            ' --pdb ' + full_filename + 
            ' --popsOut ' + self.working_path + full_filename.split('/')[-1].replace('pdb', 'out'))
        childProcess = subprocess.Popen(system_command, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True, 
                                        )
        result, e = childProcess.communicate()
        rc = childProcess.returncode
        
        # the returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        output = [ line for line in result.split('\n') if line != '' ]
        
        self.log.debug('result: %s' % result.replace('\n','; '))
        self.log.debug('output: %s' % output.replace('\n','; '))
        return output[-1], rc, e


    def __read_pops_area(self, filename):
        # The old way
        keep = ['hydrophobic:', 'hydrophilic:', 'total:']
        with open(filename, 'r') as pops:
            result = [ x.split(' ') for x in pops.readlines() if x != '' and x.split(' ')[0] in keep ]
        return [ [ x.strip() for x in item if x != '' ] for item in result ]
        

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
        structure = get_PDB(pdbCode, self.pdbPath)
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
    
    ###########################################################################
    

    

    
if __name__ == '__main__':
    import logging
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    
    pdb_file = 'Q9Y6K1_Q9UBC3.BL00030001.pdb'
    data_path = '/home/kimlab1/database_data/elaspic/human/Q9Y/6K/Q9Y6K1/PWWP*291-374/PWWP*224-307/Q9UBC3/'
    working_path = '/tmp/elaspic/NUORFs/analyze_structure/'
    chain = ['A','B']
    
    analyze_structure = AnalyzeStructure(data_path, working_path, pdb_file, chain, None, logger)
    dssp_score = analyze_structure.get_dssp()
    sasa_score = analyze_structure.get_sasa()
    interchain_distances = analyze_structure.get_interchain_distances('A_Q10N')
        
    print dssp_score[0]['A'][10]
    print sasa_score[0]['A'][10]
    print interchain_distances['A'][0]
    
    
    
    
    