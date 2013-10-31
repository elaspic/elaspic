# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:29:13 2013

@author: niklas
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:35:48 2013

@author: niklas
"""

import subprocess
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import DSSP

class interfaceSize():
    """
    Runs the program pops to calculate the interface size of the complexes
    This is done by calculating the surface of the complex and the seperated parts.
    The interface is then given by the substracting
    """
    def __init__(self, pdbPath, tmpPath):
        self.pdbPath = pdbPath
        self.tmpPath = tmpPath
    
    def __get_structures(self, pdbFile, chains):
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        io = PDBIO()

        # save chain, i.e. part one of the complex:
        structure = parser.get_structure('ID', pdbFile)
        model = structure[0]
        for child in model.get_list():
            if child.id not in chains:
                model.detach_child(child.id)
        io.set_structure(structure)
        outFile = self.tmpPath + 'chain.pdb'
        io.save(outFile)
        
        # save opposite chain, i.e. the other part of the complex:
        structure = parser.get_structure('ID', pdbFile)
        model = structure[0]
        for child in model.get_list():
            if child.id in chains:
                model.detach_child(child.id)
        io.set_structure(structure)
        outFile = self.tmpPath + 'oppositeChain.pdb'
        io.save(outFile)
        return
        
    def __run_pops(self, pdb):
        system_command = self.tmpPath + 'pops --pdb ' + pdb + ' --chainOut --popsOut ' + self.tmpPath + 'pops.out'
        childProcess = subprocess.Popen(system_command, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True, 
                                        )
        result, error = childProcess.communicate()
        rc = childProcess.returncode
        
        # the returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        output = [ line for line in result.split('\n') if line != '' ]
        return output[-1], rc, error

    def __read_pops(self, f):
        keep = ['hydrophobic:', 'hydrophilic:', 'total:']
        with open(f, 'r') as pops:
            result = [ x.split(' ') for x in pops.readlines() if x != '' and x .split(' ')[0] in keep ]
        return [ [ x.strip() for x in item if x != '' ] for item in result ]
            
            
    def __call__(self, pdb, chains):
        self.__get_structures(pdb, chains)
        
        # calculate SASA for the full complex:
        termination, rc, error = self.__run_pops(pdb)
        if rc != 0:
            if termination != 'Clean termination':
                print 'Error in pops for pdb', pdb
                print error
                return '0', '0', '0'
            else:
                print 'Warning in pops for pdb', pdb
                print error
        result = self.__read_pops(self.tmpPath + 'pops.out')
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
        termination, rc, error = self.__run_pops(self.tmpPath + 'chain.pdb')
        if rc != 0:
            if termination != 'Clean termination':
                print 'Error in pops for pdb', pdb
                print error
                return '0', '0', '0'
            else:
                print 'Warning in pops for pdb', pdb
                print error
        result = self.__read_pops(self.tmpPath + 'pops.out')
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
        termination, rc, error = self.__run_pops(self.tmpPath + 'oppositeChain.pdb')
        if rc != 0:
            if termination != 'Clean termination':
                print 'Error in pops for pdb', pdb
                print error
                return '0', '0', '0'
            else:
                print 'Warning in pops for pdb', pdb
                print error
        result = self.__read_pops(self.tmpPath + 'pops.out')
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
    

class getDSSP():
    
    def __init__(self, program_name='dssp-2.0.4-linux-amd64'):
        
        self.program_name = program_name
    
    
    def __call__(self, pdbFile, chain, residue, position):
        """
        Make sure position is a Biophython position entry like
        <Residue ASP het=  resseq=15 icode= >, i.e. (' ', 15, ' ')
        """
        if len(residue) == 1:
            residue = self.convert_aa(residue)
        
        model = self.get_structures(pdbFile)
        dssp = DSSP(model, pdbFile, self.program_name)
        
        result = dssp[chain, (' ', int(position), ' ')]

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
                print  'Not a valid amino acid'
                return
        if len(aa) == 1:
            try:
                return A_DICT[aa.upper()]
            except KeyError:
                print  'Not a valid amino acid'
                return
        print 'Not a valid amino acid'
    
    
    def get_structures(self, pdbFile):
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.

        # save chain, i.e. part one of the complex:
        structure = parser.get_structure('ID', pdbFile)
        model = structure[0]
        
        return model
    
    
    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    
    pops = interfaceSize('tmp/', 'tmp/')
    interface_size = pops('tmp/Mut_R_RepairPDB_1LFDAB_1-HETATM.pdb', ['A', ])
    interface_size1 = pops('tmp/Mut_R_RepairPDB_1LFDAB_1.pdb', ['A', ])
    
    print interface_size
    print interface_size1
    
    
    