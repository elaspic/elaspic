# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 09:49:26 2013

@author: niklas
"""
from os import environ
import shutil
import errors
import helper_functions as hf

names_rows_stability = [
    ['dg', 1], # totalEnergy
    ['backbone_hbond', 2],
    ['sidechain_hbond', 3],
    ['van_der_waals', 4],
    ['electrostatics', 5],
    ['solvation_polar', 6],
    ['solvation_hydrophobic', 7],
    ['van_der_waals_clashes', 8],
    ['entropy_sidechain', 9],
    ['entropy_mainchain', 10],
    ['sloop_entropy', 11],
    ['mloop_entropy', 12],
    ['cis_bond', 13],
    ['torsional_clash', 14],
    ['backbone_clash', 15],
    ['helix_dipole', 16],
    ['water_bridge', 17],
    ['disulfide', 18],
    ['electrostatic_kon', 19],
    ['partial_covalent_bonds', 20],
    ['energy_ionisation', 21],
    ['entropy_complex', 22],
    ['number_of_residues', 23],]
names_stability_wt = (
    [name + '_wt' for name in zip(*names_rows_stability)[0][:-1]] +
    ['number_of_residues'])
names_stability_mut = (
    [name + '_mut' for name in zip(*names_rows_stability)[0][:-1]] +
    ['number_of_residues'])

names_rows_stability_complex = (
    [ ['intraclashes_energy_1', 3], ['intraclashes_energy_2', 4], ] +
    [ [x[0], x[1] + 4] for x in names_rows_stability ])
names_stability_complex_wt = (
    [name + '_wt' for name in zip(*names_rows_stability_complex)[0][:-1]] +
    ['number_of_residues'])
names_stability_complex_mut = (
    [name + '_mut' for name in zip(*names_rows_stability_complex)[0][:-1]] +
    ['number_of_residues'])


class FoldX():

    def __init__(self, tmp_path, pdb_file, chain_id, buildModel_runs, foldX_WATER, logger):
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
        self.logger = logger




    def __call__(self, whatToRun, mutCodes=[]):
        """
        Select which action should be performed by FoldX by setting 'whatToRun'
        Possible values are: 'AnalyseComplex', 'Stability', 'RepairPDB', 'BuildModel'
        See the FoldX manual for an explanation on what they do
        c.f. (http://foldx.crg.es/manual3.jsp)
        """
        self.logger.debug('Running FoldX {}'.format(whatToRun))
        self.__write_runfile(self.pdb_filename, self.chain_id, whatToRun, mutCodes)
        self.__run_runfile()
        if whatToRun == 'AnalyseComplex':
            return self.__read_result(self.foldx_path + 'Interaction_AnalyseComplex_resultFile.txt', self.pdb_filename, whatToRun)
        elif whatToRun == 'Stability':
            return self.__read_result(self.foldx_path + 'Stability.txt', self.pdb_filename, whatToRun)
        elif whatToRun == 'RepairPDB':
            return self.foldx_path + 'RepairPDB_' + self.pdb_filename
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

        foldX_runfile = (
            '<TITLE>FOLDX_runscript;\n'
            '<JOBSTART>#;\n'
            '<PDBS>{pdbFile};\n'
            '<BATCH>#;\n'
            '<COMMANDS>FOLDX_commandfile;\n'
            '{command_line}\n'
            '<END>#;\n'
            '<OPTIONS>FOLDX_optionfile;\n'
            '<Temperature>298;\n'
            '<R>#;\n'
            '<pH>7;\n'
            '<IonStrength>0.050;\n'
            '<numberOfRuns>{buildModel_runs};\n'
            '<water>{water};\n'
            '<metal>-CRYSTAL;\n'
            '<VdWDesign>2;\n'
            '<pdb_waters>false;\n'
            '<OutPDB>{output_pdb};\n'
            '<pdb_hydrogens>false;\n'
            '<END>#;\n'
            '<JOBEND>#;\n'
            '<ENDFILE>#;\n').replace(' ', '').format(
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
#        system_command = './FoldX.linux64 -runfile ' + self.foldx_runfile
        my_env = environ.copy()
        my_env['LD_PRELOAD'] = self.foldx_path + 'libfaketime.so.1'
        my_env['FAKETIME'] = "-1y"
        system_command = './foldx64Linux -runfile ' + self.foldx_runfile
        self.logger.debug('FoldX system command: {}'.format(system_command))
        childProcess = hf.run_subprocess_locally(self.foldx_path, system_command, env=my_env)
        result, error_message = childProcess.communicate()
        return_code = childProcess.returncode
        if return_code != 0:
            self.logger.debug('FoldX result: %s' % result)
            self.logger.debug('FoldX error: %s' % error_message)
            if 'Cannot allocate memory' in error_message:
                raise errors.ResourceError(error_message)


    def __read_result(self, outFile, pdb, whatToRead):
        with open(outFile, 'r') as f:
            lines = f.readlines()
            line = lines[-1].split('\t')
            if whatToRead == 'BuildModel':
                total_energy_difference = line[1]
                return total_energy_difference
            if whatToRead == 'Stability':
                stability_values = [ line[x[1]].strip() for x in names_rows_stability ]
                return stability_values
            if whatToRead == 'AnalyseComplex':
                complex_stability_values = [ line[x[1]].strip() for x in names_rows_stability_complex ]
                return complex_stability_values


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
    fX_wt = FoldX(tmp_path + unique + '/', repairedPDB_wt, chains_modeller[0], unique,
                  buildModel_runs, foldX_WATER, logger)
    # do the mutation with foldX
    repairedPDB_wt_list, repairedPDB_mut_list = fX_wt.run('BuildModel', mutCodes)
    logger.debug('repairedPDB_wt_list: %s' % str(repairedPDB_wt_list))




