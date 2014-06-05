# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 15:17:07 2013

@author: niklas
"""
from Bio.PDB.PDBParser import PDBParser
from math import sqrt
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from difflib import SequenceMatcher
import pdb_template

#from Bio.SeqRecord import SeqRecord
#from Bio import AlignIO
#from Bio import SeqIO
#import subprocess
#import os


#import shutil
#import subprocess
#import os



        

class pysiChem():
    
    def __init__(self, vdW, d, unique, log):
        self.vdW_distance = float(vdW)
        self.contact_distance = float(d)
        self.unique = unique
        self.log = log
    
    
    def __call__(self, pdbFile, chains, chainIDs_complex, mutation_position, mutation):
        """
        return the atomic contact vector (as list with four items)
        item 0: electrostatic of equal sign
        item 1: electrostatic of opposite sign 
        item 2: hydrogen bond
        item 3: van der Waals interactions
        """
        return self.getVector(pdbFile, chains, chainIDs_complex, mutation_position, mutation)
        
        
    def __getPDB(self, pdbFile):
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        return parser.get_structure('ID', pdbFile)
    
    

    

    def __distance(self, atom1, atom2):
        """
        returns the distance of two points in three dimensional space
        
        input: atom instance of biopython
        """
        a = atom1.coord
        b = atom2.coord
        assert(len(a) == 3 and len(b) == 3)
        return sqrt(sum( (a - b)**2 for a, b in zip(a, b)))
        
        
    def __what_am_I(self, residue, atom):
        """
        Checks what type of atom it is (i.e. charged, polar, carbon)
        
        In order to see what "interaction" type two atoms are forming, check the
        individual label which every atom in every residue has (see pdb file
        convention for an explanation of the labels).
        With this label, one can determine which atom of the residue one is looking
        at, and hence, one can determine which "interaction" two atoms are forming.
        
        """
        # This is based on the naming convention for the atoms in crystalography
        # Each atom has a unique label (ask Joan he knows more)
        # label of positvely charged atoms
        charged_plus  = ['NH1', 'NH2', 'NZ']
        # label of negatively charged atoms
        charged_minus = ['OD1', 'OD2', 'OE1', 'OE2']
        # label of polar atoms
        polar         = ['OG', 'OG1', 'OD1', 'OD2', 'ND1', 'OE1', 'NE', 'NE1', \
                         'NE2', 'ND1', 'ND2', 'SG', 'OH', 'O', 'N']
        
        if residue.upper() in ['ARG', 'R', 'LYS', 'K']:
            if atom.name in charged_plus:
                return 'charged_plus'
            if atom.name in polar:
                return 'polar'
            if atom.name[0] == 'C' or atom.name == 'SD':
                return 'carbon'
            return 'ignore'
        
        elif residue.upper() in ['ASP', 'D', 'GLU', 'E']:
            if atom.name in charged_minus:
                return 'charged_minus'
            if atom.name in polar:
                return 'polar'
            if atom.name[0] == 'C' or atom.name == 'SD':
                return 'carbon'
            return 'ignore'
            
        else:
            if atom.name in polar:
                return 'polar'
            if atom.name[0] == 'C' or atom.name == 'SD':
                return 'carbon'
            return 'ignore'

        
    def find_position(self, chain, mutation_position, mutation_seq):
        """
        aligns the short sequence snipped to the chain sequence and
        determines the position
        """
        # extract the sequence of the chain
        sequence = Seq('', IUPAC.protein)
        for pb in PPBuilder().build_peptides(chain, aa_only=False):
            tmp = sequence + pb.get_sequence()
            sequence = tmp
        
        # convert the biopython objects into strings
        chain_seq       = str(sequence)
        mutation_seq    = str(mutation_seq)

        # Match the sequences:
        # Gives you the position where the shorter string is
        # found within the longer string
        s = SequenceMatcher(None, chain_seq, mutation_seq)
        match = s.get_matching_blocks()
#        if len(match) != 2:
#            print 'match not equals two!'
#            print 'matching:'
#            print 'mutation_seq', mutation_seq
#            print 'chain_seq', chain_seq
#            print 'match', match
#            print 'len(match)', len(match)
        
#        assert( len(match) == 2 ) # if the mutation_seq is to short several hits might be found
        
        return match[0][0] + int(mutation_position)
    


    def getVector(self, pdbFile, chainID, chainIDs_complex, mutation_position, mutation):
        """
        Return the atomic contact vector, that is, counting how many interactions
        between charged, polar or "carbon" residues there are. The "carbon"
        interactions give you information about the Van der Waals packing of
        the residues. Comparing the wildtype vs. the mutant values is used in
        the machine learning algorithm.
        
        'mutation' is of the form: 'A16' where A is the chain identifier and 16
        the residue number (in pdb numbering) of the mutation 
        chainIDs is a list of strings with the chain identifiers to be used
        if more than two chains are given, the chains not containing the mutation
        are considered as "opposing" chain
        """
        mainChain = ['CA', 'C', 'N', 'O']
        structure = self.__getPDB(pdbFile)
        model = structure[0]
        
        chains = [ c for c in model ]
        chain          = model[chainID]
        oppositeChains = [ c for c in chains if c.id not in chainIDs_complex ]
        
        residue_counter = 0
        for residue in chain:
            if residue.resname in pdb_template.amino_acids and residue.id[0] == ' ':
                residue_counter += 1
                if residue_counter == mutation_position:
                    if not (pdb_template.AAA_DICT[residue.resname] == mutation.upper()):
                        self.log.error(residue)
                        self.log.error(pdb_template.AAA_DICT[residue.resname])
                        self.log.error(mutation.upper())
                        self.log.error(mutation_position)
                        raise Exception
                    mutated_residue = residue
                    break
                
#            chainNumbering.append(residue.id[1])
#            chain_sequence.append(AAA_DICT[residue.resname])
    
    # DELETE THE FOLLOWING LINES IF EVERYTHING WORKS
#        # get the sequence position of the mutation in the modelled structure
##        mutPosition = self.find_position(chain, mutation_position, mutation)
#        mutPosition = mutation_position #AS
#        # get the pdb numbering of the modelled structure
#        chainNumbering, chain_sequence = analyze_structure.get_chain_numbering(chain, True)
#        # the mutated residue is determined by getting the numbering of the pdb
#        # file at the corresponding sequence position
#        mutatedResidue = chain[chainNumbering[mutPosition-1]]
#
#        #######################################################################
#        
#        self.log.debug('mutatedResidue: %s' % mutatedResidue)
#        self.log.debug('mutatedPosition in chain: %s' % chain[chainNumbering[mutPosition-1]])
#        self.log.debug('chain sequence:') analyze_structure.AAA_DICT[chain[chainNumbering[mutPosition-1]].resname]
#        self.log.debug(chain_sequence)
#        assert(analyze_structure.AAA_DICT[mutatedResidue.resname] == mutation)

        #######################################################################
        
        mutatedResidue = mutated_residue
        
        mutatedAtoms = list()
        for atom in mutatedResidue:
            if atom.name not in mainChain:
                mutatedAtoms.append(atom)
        
        # atomicContactVector = electrostatic_equal, electrostatic_opposite, hydrogen bond, vdW
        atomicContactVector = [ 0, 0, 0, 0 ]
        atomicContactVector_carbon_count = set()
        # if also the charged and polar thingys should be counted nonredundantly uncomment
#        atomicContactVector_charge_equal_count = set()
#        atomicContactVector_charge_opposite_count = set()
#        atomicContactVector_hbond_count = set()
        count = 0
        count_charge_equal = 0
        count_charge_opposite = 0
        count_hbond_count = 0
        ## go through each chain
        for oppositeChain in oppositeChains:
            # check each residue
            for residue in oppositeChain:
                # for each residue each atom of the mutated residue has to be checked
                for mutatedAtom in mutatedAtoms:
                    mutatedAtomType = self.__what_am_I(mutatedResidue.resname, mutatedAtom)
                    # and each atom
                    for atom in residue:
                        r = self.__distance(mutatedAtom, atom)
                        if r <= self.vdW_distance:
                            count +=1
                            I_am = self.__what_am_I(residue.resname, atom)
                            if I_am == 'ignore':
                                continue
                            # to avoid duplicate counts, a set of the atom coordinates
                            # is used to keep track of the number interactions
                            unique = str(atom.coord[0]) + str(atom.coord[1]) + str(atom.coord[2])
                            if mutatedAtomType == 'carbon' and I_am == 'carbon':
                                # The Van der Waals packing should be determined
                                # nonredundant. Thus, the atomic coordinates are
                                # used to keep track of which interactions where
                                # already counted. Can easily be implemented for
                                # the other types by uncommenting the lines below.
                                atomicContactVector_carbon_count.add(unique)
                            
                            if r <= self.contact_distance:
                                unique_tmp = unique + str(mutatedAtom.coord[0]) + str(mutatedAtom.coord[1]) + str(mutatedAtom.coord[2])
                                unique = unique_tmp
                                if mutatedAtomType == 'charged_plus' and I_am == 'charged_plus':
                                    count_charge_equal += 1
#                                    atomicContactVector_charge_equal_count.add(unique)
                                if mutatedAtomType == 'charged_minus' and I_am == 'charged_plus':
                                    count_charge_opposite += 1
#                                    atomicContactVector_charge_opposite_count.add(unique)
                                if mutatedAtomType == 'charged_plus' and I_am == 'charged_minus':
                                    count_charge_opposite += 1
#                                    atomicContactVector_charge_opposite_count.add(unique)
                                if mutatedAtomType == 'charged_minus' and I_am == 'charged_minus':
                                    count_charge_equal += 1
#                                    atomicContactVector_charge_equal_count.add(unique)
                                    
                                if mutatedAtomType == 'charged' and I_am == 'polar':
                                    count_hbond_count += 1
#                                    atomicContactVector_hbond_count.add(unique)
                                if mutatedAtomType == 'polar' and I_am == 'charged':
                                    count_hbond_count += 1
#                                    atomicContactVector_hbond_count.add(unique)
                                if mutatedAtomType == 'polar' and I_am == 'polar':
                                    count_hbond_count += 1
#                                    atomicContactVector_hbond_count.add(unique)

#        atomicContactVector[0] = len(atomicContactVector_charge_equal_count)
#        atomicContactVector[1] = len(atomicContactVector_charge_opposite_count)
#        atomicContactVector[2] = len(atomicContactVector_hbond_count)
        
        atomicContactVector[0] = count_charge_equal
        atomicContactVector[1] = count_charge_opposite
        atomicContactVector[2] = count_hbond_count
        atomicContactVector[3] = len(atomicContactVector_carbon_count)
        
        
        atomicContactVector_ownChain = [ 0, 0, 0, 0 ]
        atomicContactVector_carbon_count_ownChain = set()
        count_ownChain = 0
        count_charge_equal_ownChain = 0
        count_charge_opposite_ownChain = 0
        count_hbond_count_ownChain = 0
        # calculate the vector for the own chain
        for residue in chain:
            # if the don't count the mutated residue
            if residue == mutatedResidue:
                continue
            # for each residue each atom of the mutated residue has to be checked
            for mutatedAtom in mutatedAtoms:
                mutatedAtomType = self.__what_am_I(mutatedResidue.resname, mutatedAtom)
                # and each atom
                for atom in residue:
                    r = self.__distance(mutatedAtom, atom)
                    if r <= self.vdW_distance:
                        count_ownChain +=1
                        I_am = self.__what_am_I(residue.resname, atom)
                        if I_am == 'ignore':
                            continue
                        # to avoid duplicate counts, a set of the atom coordinates
                        # is used to keep track of the number interactions
                        unique = str(atom.coord[0]) + str(atom.coord[1]) + str(atom.coord[2])
                        if mutatedAtomType == 'carbon' and I_am == 'carbon':
                            atomicContactVector_carbon_count_ownChain.add(unique)
                        
                        if r <= self.contact_distance:
                            unique_tmp = unique + str(mutatedAtom.coord[0]) + str(mutatedAtom.coord[1]) + str(mutatedAtom.coord[2])
                            unique = unique_tmp
                            if mutatedAtomType == 'charged_plus' and I_am == 'charged_plus':
                                count_charge_equal_ownChain += 1
                            if mutatedAtomType == 'charged_minus' and I_am == 'charged_plus':
                                count_charge_opposite_ownChain += 1
                            if mutatedAtomType == 'charged_plus' and I_am == 'charged_minus':
                                count_charge_opposite_ownChain += 1
                            if mutatedAtomType == 'charged_minus' and I_am == 'charged_minus':
                                count_charge_equal_ownChain += 1
                                
                            if mutatedAtomType == 'charged' and I_am == 'polar':
                                count_hbond_count_ownChain += 1
                            if mutatedAtomType == 'polar' and I_am == 'charged':
                                count_hbond_count_ownChain += 1
                            if mutatedAtomType == 'polar' and I_am == 'polar':
                                count_hbond_count_ownChain += 1

#        atomicContactVector[0] = len(atomicContactVector_charge_equal_count)
#        atomicContactVector[1] = len(atomicContactVector_charge_opposite_count)
#        atomicContactVector[2] = len(atomicContactVector_hbond_count)
        
        atomicContactVector_ownChain[0] = count_charge_equal_ownChain
        atomicContactVector_ownChain[1] = count_charge_opposite_ownChain
        atomicContactVector_ownChain[2] = count_hbond_count_ownChain
        atomicContactVector_ownChain[3] = len(atomicContactVector_carbon_count_ownChain)
        
        return atomicContactVector, atomicContactVector_ownChain
        

if __name__ == '__main__':
    test = pysiChem(6,5, '/tmp/test/')
    
    print test('WT_RepairPDB_1DFJ_IE_0_Y434Achain2.BL00010003_1.pdb', 'C', 'QLVLYDIYWSE')
    
    print test('RepairPDB_1DFJ_IE_0_Y434Achain2.BL00010003_1.pdb', 'C', 'QLVLADIYWSE')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
