# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 18:58:50 2012

@author: niklas
"""


import gzip
import subprocess

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import class_error



class pdbTemplate():
    
    def __init__(self, pdbPath, pdbCode, chains, domainBoundaries, outputPath):
        """
        fetches the pdb file from the local database, extracts the desired chains
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
        self.pdbCode = pdbCode
        self.outputPath = outputPath
        
        self.chains = chains
        self.domainBoundaries = domainBoundaries

        childProcess = subprocess.Popen('echo $TMPDIR', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, e = childProcess.communicate()
        self.tmpPath = result.strip() + '/tmp/'
        
        self.pdbStructure = self.__getPDB(self.pdbCode)
        
       
    def extract(self, returnHETATM=False):
        """
        extract and save the pdb data and the sequences of the desired chains
        
        input:  returnHETATM    type boolean    if True returns dict with elements
                                                of self.chains as keys and the
                                                corresponding number of HETATM
                                                residues as values
        """
        # first, remove the chains from the structure instance that are not
        # needed. HETATM have to be included as '.' in the modeller .pir file
        # HETATMs can not only occur at the end of a chain but also somewhere
        # in the middle. By recording where they are, one can later determine
        # if a '.' has to be included in the sequence
        # Modeller needs the chains in the order in which they appear in the pdb
        # file. If the ordering in the pdb differs from the given target, one
        # has to switch them.
        HETATMpositions, chains_pdb_order = self.__extractChains(self.pdbStructure, self.chains, self.domainBoundaries)

        # HETFlag is a dictionary with key as chainID and the value True/False
        # depending wether HETATMs are found for chainID.
        HETFlag = dict()
        for chainID in chains_pdb_order:
            chainNumbering = self.getChainNumberingNOHETATMS(chainID)
            if chainNumbering == []:
                raise class_error.PDBChainError(self.pdbCode, self.chains)
                
            seq = self.extractSequence(chainID)
            
            HETFlag[chainID] = False
            if HETATMpositions[chainID] != []:
                for item in HETATMpositions[chainID]:
                    if item not in range(min(chainNumbering),max(chainNumbering)+1):
                        HETFlag[chainID] = True
            
            # write the sequence to a file
            with open(self.outputPath + str(self.pdbCode) + str(chainID) + '.seq.txt', 'w') as f:
                f.write('>' + str(self.pdbCode) + str(chainID) + '\n')
                f.write(str(seq) + '\n')
                f.write('\n')
        
        return HETATMpositions, HETFlag, chains_pdb_order
        
    
    def getChainNumberingNOHETATMS(self, chainID):
        """
        extracts the interface of a given chain
        
        input:  chainID         type string
        
        return  chainNumbering  type list   containing the atom numbering from the pdb file
        """
        model = self.pdbStructure[0]
        chain = model[chainID]
        
        chainNumbering = list()
        amino_acids = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', \
                       'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', \
                       'MET', 'PHE', 'TYR', 'TRP']
                       
        for residue in chain:
            if residue.resname in amino_acids and residue.id[0] == ' ':
                chainNumbering.append(residue.id[1])
        
        return chainNumbering
    
        
    def getChainNumbering(self, chainID):
        """
        extracts the interface of a given chain
        
        input:  chainID         type string
        
        return  chainNumbering  type list   containing the atom numbering from the pdb file or residue.id
        """
        model = self.pdbStructure[0]
        chain = model[chainID]
        
        chainNumbering = list()

        for residue in chain:
            chainNumbering.append(residue.id[1])
        
        return chainNumbering

            
    def __getPDB(self, pdbCode):
        """
        parse a pdb file with biopythons PDBParser() and return the structure
        
        input: pdbCode  type String     four letter code of the PDB file
        
        return: Biopython pdb structure
        """
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        pdbFile = self.pdbPath + pdbCode[1:3].lower() + '/pdb' + pdbCode.lower() + '.ent.gz'
        try:
            pdbFileUncompressed = gzip.open(pdbFile, 'r')
        except IOError:
            raise class_error.NoPDBFound('pdb' + pdbCode.lower() + '.ent.gz')
            
        return parser.get_structure('ID', pdbFileUncompressed)
        
    def __extractChains(self, pdbStructure, chains, domainBoundaries):
        """
        extract the wanted chains out of the PDB file
        removes water atoms and selects the domain regions
        
        input:  pdbStructure    type: Biopython pdb structure
                chains          type: list  containing strings with chain names
                domainBoundarys type: list  containing a tuple of beginning and ending for
                                            each chain in chains
        """
        io = PDBIO()
        chains_pdb_order = list()
        HETATMpositions = dict()

        model = pdbStructure[0] # assuming that model 0 is always the desired one
                                # should be true for most cases
        
        # detach the chain from the model if the ID is not found in the list
        # and add the chain id to a list to record the order
        for child in model.get_list():
            if child.id not in chains:
                model.detach_child(child.id)
            else:
                chains_pdb_order.append(child.id)
        

        # remove the water molecules and select the domain part
        i = -1
        for chainID in chains:
            i += 1
            chain = model[chainID]
            residues = list()
#            residues_numbering = list()
                
            # add the IDs of the form (' ', 123, ' ') for a list
            # the first entry specifies the type, the second is the position,
            # the third the icode (see Biopython manual)
            for residue in chain:
                residues.append(residue.id)
#                if residue.id[0] == ' ':
#                    residues_numbering.append(residue.id[1])
            
            HETATMpositions[chainID] = []
            BEGIN = False
            END = False
            for residue in residues:
                ## add special treatment for some PDB files...
                ## not elegant but as a first try
                # this last residue is not connected to the main chain and not
                # complete. It is not considered in the sequence and is thus
                # removed in the PDB file for modeller
                if self.pdbCode == '1H9D' and residue[1] == 178:
                    chain.detach_child(residue)
                    continue
                # in chain A residue 1010-1016 are disconected and not needed for
                # the structure. But they cause problems
                # in chain B redisue 1004 is incomplete
                if self.pdbCode == '2JIU' and residue[1] in [1004, 1010, 1011, 1012, 1013, 1014, 1015, 1016]:
                    chain.detach_child(residue)
                    continue
                
                #remove water
                if residue[0] == 'W':
                    chain.detach_child(residue)
                # select the domain only and remove the rest
                else:
                    # record the position of the HETATM
                    if residue[0] != ' ':
                        # H_MSE is Selenomethionine and is treated as M by modeller
                        if residue[0] == 'H_MSE':
                            pass
                        else:
                            HETATMpositions[chainID].append(residue[1])
                    
                    # I had at least one case where a non-standard residue was
                    # not listed as HETATM in the pdb file. Thus a check if
                    # the residue is a standart residue and if not, consider
                    # it as HETATM
                    amino_acids = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', \
                                   'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', \
                                   'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', \
                                   'TYR', 'TRP']
                    try:                    
                        if chain[residue[1]].resname not in amino_acids:
                            HETATMpositions[chainID].append(residue[1])
                    except KeyError:
                        pass

                    if domainBoundaries[i] != ['Null', 'Null']:
                        # to select the domain boundaries one has to start counting
                        # when the first number appears and has to stop after the
                        # second. PDB files are not necessarily numbered strictly 
                        # ascending!
                        #
                        # residue[0] == ' ' is for ignoring all HETATMs
                        # this way they are kept in the pdb file
                                                
                        # Selenomethionine is consideres a normal residue and
                        # should be removed if outside of the domain boundaries
                        if residue[0] == ' ' and residue[1] >= int(domainBoundaries[i][0]):
                            BEGIN = True # from now on, keep the residues
                        if residue[0] == 'H_MSE' and residue[1] >= int(domainBoundaries[i][0]):
                            BEGIN = True # from now on, keep the residues
                        
                        if residue[0] == ' ' and residue[1] >= int(domainBoundaries[i][1]) + 1:
                            END = True # now detach them again
                        if residue[0] == 'H_MSE' and residue[1] >= int(domainBoundaries[i][1]) + 1:
                            END = True # now detach them again
                        
                        if BEGIN == False: # we are before the beginning
                                           # i.e. remove the residue
                            if residue[0] == ' ' or residue[0] == 'H_MSE':
                                chain.detach_child(residue)
                        elif END: # we are after the end
                                  # i.e. remove the residue
                            if residue[0] == ' ' or residue[0] == 'H_MSE':
                                chain.detach_child(residue)
                        # else: we are inbetween the domain boundaries
                        #       i.e. keep the residues

        # save the structure to a pdb file
        io.set_structure(pdbStructure)
        outFile = self.outputPath + self.pdbCode + ''.join(chains_pdb_order) + '.pdb'
        io.save(outFile)
        
        return HETATMpositions, chains_pdb_order
    
    
                        
    def extractSequence(self, chainID):
        """
        Extracts a sequence from a PDB file.
        Usefull when interested in the sequence that was used for crystallization
        and not the SEQRES sequence
        
        Note that chain must be a single item, not a list!
        
        input:  pdbStructure    type: Biopython PDB structure
                chainID         type: string
        """
        # sanity checks
        if chainID.strip() == '':
            return ''
        if isinstance(chainID, list):
            if len(chainID) == 1:
                chainID = chainID[0]
            else:
                print 'extractSequence takes a string with chainID! Not a list'
                return ''
        
        model = self.pdbStructure[0]
        chain = model[chainID]

        # setting aa_only to False Selenomethionines are reported in the
        # sequence as well
        # see: http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html
        sequence = Seq('', IUPAC.protein)
        for pb in PPBuilder().build_peptides(chain, aa_only=False):
#        for pb in PPBuilder().build_peptides(chain, aa_only=True):
            tmp = sequence + pb.get_sequence()
            sequence = tmp
        return sequence







        
        
if __name__ == '__main__':
    pdbPath = '/home/niklas/pdb_database/structures/divided/pdb/'
    pdbCode = '1VOK'
    chains = ['A', 'B']
    domainBoundaries = [[23, 115], [29, 198]]
    outputPath = './'
    p = pdbTemplate(pdbPath, pdbCode, chains, domainBoundaries, outputPath)
    print p.extract()
