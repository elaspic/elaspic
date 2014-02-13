# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 09:49:26 2013

@author: niklas
"""

import subprocess
import shutil
from class_error import FoldXError

class foldX():
    
    def __init__(self, tmpPath, pdbFile, chainID, unique, buildModel_runs, foldX_WATER):

        # In case the pipeline gets extended to handle more than two chains,
        # the chainID becomes relevant for the energy calculation. Otherwise
        # it is OK to simply take the first chain
        self.chainID = chainID
#        self.chainID = 'A'
        
        self.PATH = tmpPath + 'FoldX/'
        self.runFile = self.PATH + 'runfile_FoldX.txt'
        self.executablePath = self.PATH
        
        self.pdbFile = pdbFile.split('/')[-1]
        
        self.buildModel_runs = buildModel_runs
        self.water = foldX_WATER
        try:
            subprocess.check_call('cp ' + pdbFile + ' ' + self.pdbFile, shell=True)
        except:
            pass
        
        
    def run(self, whatToRun, mutCodes=[]):
        """
        Select which action should be performed by FoldX by setting 'whatToRun'
        Possible values are: 'AnalyseComplex', 'Stability', 'RepairPDB', 'BuildModel'
        See the FoldX manual for an explanation on what they do
        c.f. (http://foldx.crg.es/manual3.jsp)
        """
        if whatToRun == 'AnalyseComplex':
            self.__writeRunfileAnalyseComplex(self.pdbFile, self.chainID)
            shutil.copy(self.runFile, self.PATH + 'run-analyseComplex.txt')
        if whatToRun == 'Stability':
            self.__writeRunfileStability(self.pdbFile)
            shutil.copy(self.runFile, self.PATH + 'run-stability.txt')
        if whatToRun == 'RepairPDB':
            self.__writeRunfileRepairPDB(self.pdbFile)
            shutil.copy(self.runFile, self.PATH + 'run-repair.txt')
        if whatToRun == 'BuildModel':
            self.__writeRunfileBuildModel(self.pdbFile, mutCodes)
            shutil.copy(self.runFile, self.PATH + 'run-build.txt')
            
        
        rc, error = self.__executeRunfile(self.runFile, whatToRun)
        if rc == 0:
            if whatToRun == 'AnalyseComplex':
                return self.__readResult(self.PATH + 'Interaction_AnalyseComplex_resultFile.txt', self.pdbFile, whatToRun)
            if whatToRun == 'Stability':
                return self.__readResult(self.PATH + 'Stability.txt', self.pdbFile, whatToRun)
            if whatToRun == 'BuildModel':
                if self.buildModel_runs == str(1):
                    # see the FoldX manual for the naming of the generated structures
                    mutants = [self.PATH + self.pdbFile[:-4] + '_1.pdb', ]
                    wiltype = [self.PATH + 'WT_' + self.pdbFile[:-4] + '_1.pdb', ]
                    return wiltype, mutants
                else:
                    mutants = [ self.PATH + self.pdbFile[:-4] + '_1_' + str(x) + '.pdb' for x in range(0,int(self.buildModel_runs)) ]
                    wiltype = [ self.PATH + 'WT_' + self.pdbFile[:-4] + '_1_' + str(x) + '.pdb' for x in range(0,int(self.buildModel_runs)) ]
                    return wiltype, mutants
            if whatToRun == 'RepairPDB':
                return self.PATH + 'RepairPDB_' + self.pdbFile
            return
        else:
            print error
            print rc
            print 'FoldX error!!!', whatToRun
            print 'self.pdbFile', self.pdbFile
            print 'error', error
            raise FoldXError(error)
#            return 1, error
    
    
    
    def __writeRunfileAnalyseComplex(self, pdbFile, chainID):
        #[C] = [K] - 273.15, i.e. 298K=25C
        foldX_runfile = """\
<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>""" + pdbFile + """;
<BATCH>#;
<COMMANDS>FOLDX_commandfile;
<AnalyseComplex>AnalyseComplex_resultFile.txt,""" + chainID + """;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>""" + self.water + """;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>false;
<pdb_hydrogens>false;
<END>#;
<JOBEND>#;
<ENDFILE>#;"""
             
        with open(self.runFile, 'w') as f:
            f.write(foldX_runfile)
    
    def __writeRunfileStability(self, pdbFile):
        #[C] = [K] - 273.15, i.e. 298K=25C
        foldX_runfile = """\
<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>""" + pdbFile + """;
<BATCH>#;
<COMMANDS>FOLDX_commandfile;
<Stability>Stability.txt;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>""" + self.water + """;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>false;
<pdb_hydrogens>false;
<END>#;
<JOBEND>#;
<ENDFILE>#;"""
             
        with open(self.runFile, 'w') as f:
            f.write(foldX_runfile)
    
    def __writeRunfileRepairPDB(self, pdbFile):
        #[C] = [K] - 273.15, i.e. 298K=25C
        foldX_runfile = """\
<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>""" + pdbFile + """;
<BATCH>#;
<COMMANDS>FOLDX_commandfile;
<RepairPDB>#;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>""" + self.water + """;
<metal>-CRYSTAL;
<VdWDesign>2;
<pdb_waters>true;
<OutPDB>true;
<pdb_hydrogens>false;
<END>#;
<JOBEND>#;
<ENDFILE>#;"""
             
        with open(self.runFile, 'w') as f:
            f.write(foldX_runfile)


    def __writeRunfileBuildModel(self, pdbFile, mutCodes):
        #[C] = [K] - 273.15, i.e. 298K=25C
        foldX_runfile = """\
<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>""" + pdbFile + """;
<BATCH>#;
<COMMANDS>FOLDX_commandfile;
<BuildModel>BuildModel,mutant_file.txt;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<numberOfRuns>""" + self.buildModel_runs + """;
<water>""" + self.water + """;
<metal>-CRYSTAL;
<VdWDesign>2;
<pdb_waters>false;
<OutPDB>true;
<pdb_hydrogens>false;
<END>#;
<JOBEND>#;
<ENDFILE>#;"""
             
        with open(self.runFile, 'w') as f:
            f.write(foldX_runfile)

        with open(self.PATH + 'mutant_file.txt', 'w') as f:
            mutation = ''
            FIRST = 0
            for mut in mutCodes:
                mutation = mutation + FIRST*',' + mut
                FIRST = 1
                f.write(mut + '\n')
            


    def __executeRunfile(self, runfile, whatToRun):

        system_command = './FoldX.linux64 -runfile ' + runfile
        cmd = 'cd ' + self.PATH + ' && ' + system_command
        
        childProcess = subprocess.Popen(cmd, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True, 
                                        )
        result, error = childProcess.communicate()
        
        return childProcess.returncode, error

      
    def __readResult(self, outFile, pdb, whatToRead):
        with open(outFile, 'r') as f:        
            lines = f.readlines()
            line = lines[-1].split('\t')
            if whatToRead == 'AnalyseComplex':
                intraclashesEnergy1 = line[3]
                intraclashesEnergy2 = line[4]
                interactionEnergy = line[5]
                Backbone_Hbond = line[6]
                Sidechain_Hbond = line[7]
                Van_der_Waals = line[8]
                Electrostatics = line[9]
                Solvation_Polar = line[10]
                Solvation_Hydrophobic = line[11]
                Van_der_Waals_clashes = line[12]
                entropy_sidechain = line[13]
                entropy_mainchain = line[14]
                sloop_entropy = line[15]
                mloop_entropy = line[16]
                cis_bond = line[17]
                torsional_clash = line[18]
                backbone_clash = line[19]
                helix_dipole = line[20]
                water_bridge = line[21]
                disulfide = line[22]
                electrostatic_kon = line[23]
                partial_covalent_bonds = line[24]
                energy_Ionisation = line[25]
                Entropy_Complex = line[26]
                Number_of_Residues = line[27].strip()
                FoldX_vector = intraclashesEnergy1, intraclashesEnergy2, interactionEnergy,\
                               Backbone_Hbond, Sidechain_Hbond, Van_der_Waals, Electrostatics,\
                               Solvation_Polar, Solvation_Hydrophobic, Van_der_Waals_clashes,\
                               entropy_sidechain, entropy_mainchain, sloop_entropy,\
                               mloop_entropy, cis_bond, torsional_clash, backbone_clash,\
                               helix_dipole, water_bridge, disulfide, electrostatic_kon,\
                               partial_covalent_bonds, energy_Ionisation,\
                               Entropy_Complex, Number_of_Residues
                return FoldX_vector
            if whatToRead == 'Stability':
                totalEnergy = line[1]
                Backbone_Hbond = line[2]
                Sidechain_Hbond = line[3]
                Van_der_Waals = line[4]
                Electrostatics = line[5]
                Solvation_Polar = line[6]
                Solvation_Hydrophobic = line[7]
                Van_der_Waals_clashes = line[8]
                entropy_sidechain = line[9]
                entropy_mainchain = line[10]
                sloop_entropy = line[11]
                mloop_entropy = line[12]
                cis_bond = line[13]
                torsional_clash = line[14]
                backbone_clash = line[15]
                helix_dipole = line[16]
                water_bridge = line[17]
                disulfide = line[18]
                electrostatic_kon = line[19]
                partial_covalent_bonds = line[20]
                energy_Ionisation = line[21]
                Entropy_Complex = line[22]
                Number_of_Residues = line[23].strip()
                FoldX_vector = totalEnergy,\
                               Backbone_Hbond, Sidechain_Hbond, Van_der_Waals, Electrostatics,\
                               Solvation_Polar, Solvation_Hydrophobic, Van_der_Waals_clashes,\
                               entropy_sidechain, entropy_mainchain, sloop_entropy,\
                               mloop_entropy, cis_bond, torsional_clash, backbone_clash,\
                               helix_dipole, water_bridge, disulfide, electrostatic_kon,\
                               partial_covalent_bonds, energy_Ionisation,\
                               Entropy_Complex, Number_of_Residues
                return FoldX_vector
            if whatToRead == 'BuildModel':
                totalEnergyDifference = line[1]
                return totalEnergyDifference  

        



if __name__ == '__main__':
    import os
    os.chdir('/tmp/FoldX/')
    fX = foldX('/tmp', '1LFD.pdb', 'A', 'unique', '1', '-IGNORE')
    print '--------'
#    print fX.run('RepairPDB')
    print '--------'
    print fX.run('AnalyseComplex')
    print '--------'
#    print fX.run('Stability')
    

#    print fX.readResult('AnalyseComplex_1ABO.txt', '1ABO.pdb')






