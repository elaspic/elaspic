# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 09:49:26 2013

@author: niklas
"""

import shutil
import errors
import helper_functions as hf

class foldX():

    def __init__(self, tmp_path, pdb_file, chain_id, buildModel_runs, foldX_WATER, log):
        """
        """
        # In case the pipeline gets extended to handle more than two chains,
        # the chainID becomes relevant for the energy calculation. Otherwise
        # it is OK to simply take the first chain
        self.chain_id = chain_id
#        self.chain_id = 'A'
        self.foldx_path = tmp_path + 'FoldX/'
        self.foldx_runfile = self.foldx_path + 'runfile_FoldX.txt'
        self.pdb_filename = pdb_file.split('/')[-1]
        self.buildModel_runs = buildModel_runs
        self.water = foldX_WATER
        self.log = log


    def run(self, whatToRun, mutCodes=[]):
        """
        Select which action should be performed by FoldX by setting 'whatToRun'
        Possible values are: 'AnalyseComplex', 'Stability', 'RepairPDB', 'BuildModel'
        See the FoldX manual for an explanation on what they do
        c.f. (http://foldx.crg.es/manual3.jsp)
        """
        self.log.debug('Running FoldX {}'.format(whatToRun))
        self.__write_runfile(self.pdb_filename, self.chain_id, whatToRun, mutCodes)
        self.__run_runfile()
        if whatToRun == 'AnalyseComplex':
            results = self.__readResult(self.foldx_path + 'Interaction_AnalyseComplex_resultFile.txt', self.pdb_filename, whatToRun)
        elif whatToRun == 'Stability':
            results = self.__readResult(self.foldx_path + 'Stability.txt', self.pdb_filename, whatToRun)
        elif whatToRun == 'RepairPDB':
            results = self.foldx_path + 'RepairPDB_' + self.pdb_filename
        elif whatToRun == 'BuildModel':
            # see the FoldX manual for the naming of the generated structures
            if self.buildModel_runs == '1':
                mutants = [self.foldx_path + self.pdb_filename[:-4] + '_1.pdb', ]
                wiltype = [self.foldx_path + 'WT_' + self.pdb_filename[:-4] + '_1.pdb', ]
                results = [wiltype, mutants]
            else:
                mutants = [ self.foldx_path + self.pdb_filename[:-4] + '_1_' + str(x) + '.pdb' for x in range(0,int(self.buildModel_runs)) ]
                wiltype = [ self.foldx_path + 'WT_' + self.pdb_filename[:-4] + '_1_' + str(x) + '.pdb' for x in range(0,int(self.buildModel_runs)) ]
                results = [wiltype, mutants]
        return results


    def __write_runfile(self, pdbFile, chainID, whatToRun, mutCodes):

        if whatToRun == 'AnalyseComplex':
            copy_filename = 'run-analyseComplex.txt'
            command_line = '<AnalyseComplex>AnalyseComplex_resultFile.txt,{chainID};'\
                .format(chainID=chainID)
            output_pdb = 'false'
        elif whatToRun == 'Stability':
            copy_filename = 'run-stability.txt'
            command_line = '<Stability>Stability.txt;'
            output_pdb = 'false'
        elif whatToRun == 'RepairPDB':
            copy_filename = 'run-repair.txt'
            command_line = '<RepairPDB>#;'
            output_pdb = 'true'
        elif whatToRun == 'BuildModel':
            copy_filename = 'run-build.txt'
            # file_with_mutations = 'mutant_file.txt'
            file_with_mutations = 'individual_list.txt'
            with open(self.foldx_path + file_with_mutations, 'w') as fh:
                fh.writelines(','.join(mutCodes) + ';\n')
            command_line = '<BuildModel>BuildModel,{file_with_mutations};'\
                .format(file_with_mutations=file_with_mutations)
            output_pdb = 'true'

        foldX_runfile = """\
            <TITLE>FOLDX_runscript;
            <JOBSTART>#;
            <PDBS>{pdbFile};
            <BATCH>#;
            <COMMANDS>FOLDX_commandfile;
            {command_line}
            <END>#;
            <OPTIONS>FOLDX_optionfile;
            <Temperature>298;
            <R>#;
            <pH>7;
            <IonStrength>0.050;
            <numberOfRuns>{buildModel_runs};
            <water>{water};
            <metal>-CRYSTAL;
            <VdWDesign>2;
            <pdb_waters>false;
            <OutPDB>{output_pdb};
            <pdb_hydrogens>false;
            <END>#;
            <JOBEND>#;
            <ENDFILE>#;
            """.replace(' ', '').format(
            pdbFile=pdbFile,
            command_line=command_line,
            buildModel_runs=self.buildModel_runs,
            water=self.water,
            output_pdb=output_pdb)

        # This just makes copies of the runfiles for debugging...
        with open(self.foldx_runfile, 'w') as f:
            f.write(foldX_runfile)
        shutil.copy(self.foldx_runfile, self.foldx_path + copy_filename)


    def __run_runfile(self):
        system_command = './FoldX.linux64 -runfile ' + self.foldx_runfile
        self.log.debug('FoldX system command: {}'.format(system_command))
        childProcess = hf.run_subprocess_locally(self.foldx_path, system_command)
        result, error_message = childProcess.communicate()
        return_code = childProcess.returncode
        if return_code != 0:
            self.log.debug('FoldX result: %s' % result)
            self.log.debug('FoldX error: %s' % error_message)
            if 'Cannot allocate memory' in error_message:
                raise errors.ResourceError(error_message)


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
    import logging

    unique = '6OpXOj'
    mutCodes = ['GGEALGRLLVVYPWT\nGGEALGTLLVVYPWT']
    repairedPDB_wt = '/tmp/elaspic/6OpXOj/FoldX/RepairPDB_A0N071_P02042.BL00040001.pdb'
    chains_modeller = [u'C', u'A']

    tmp_path = '/tmp/elaspic/'
    pdb_path = '/home/kimlab1/database_data/pdb/data/data/structures/divided/pdb/'

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
#    handler = logging.FileHandler(tmp_path + 'templates.log', mode='w', delay=True)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)

    foldX_path = tmp_path + unique + 'FoldX/'
    buildModel_runs = '1'
    foldX_WATER = '-IGNORE'
    log = logger
    fX_wt = foldX(tmp_path + unique + '/', repairedPDB_wt, chains_modeller[0], unique,
                  buildModel_runs, foldX_WATER, log)
    # do the mutation with foldX
    repairedPDB_wt_list, repairedPDB_mut_list = fX_wt.run('BuildModel', mutCodes)
    log.debug('repairedPDB_wt_list: %s' % str(repairedPDB_wt_list))
    log.debug('repairedPDB_mut_list: %s' % str(repairedPDB_mut_list))




