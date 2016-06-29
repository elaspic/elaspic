import os
import os.path as op
import shutil
import logging

import faketime.config
from . import conf, errors, helper

logger = logging.getLogger(__name__)

names_rows_stability = [
    ['dg', 1],  # totalEnergy
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
    ['number_of_residues', 23]
]
names_stability_wt = (
    [name + '_wt' for name in list(zip(*names_rows_stability))[0][:-1]] +
    ['number_of_residues'])
names_stability_mut = (
    [name + '_mut' for name in list(zip(*names_rows_stability))[0][:-1]] +
    ['number_of_residues'])

names_rows_stability_complex = (
    [['intraclashes_energy_1', 3], ['intraclashes_energy_2', 4], ] +
    [[x[0], x[1] + 4] for x in names_rows_stability]
)
names_stability_complex_wt = (
    [name + '_wt' for name in list(zip(*names_rows_stability_complex))[0][:-1]] +
    ['number_of_residues'])
names_stability_complex_mut = (
    [name + '_mut' for name in list(zip(*names_rows_stability_complex))[0][:-1]] +
    ['number_of_residues'])


class FoldX(object):

    def __init__(self, pdb_file, chain_id, foldx_dir=None):
        """
        """
        self.pdb_filename = op.basename(pdb_file)
        self.chain_id = chain_id
        if foldx_dir is None:
            self.foldx_dir = conf.CONFIGS['foldx_dir']
        else:
            self.foldx_dir = foldx_dir
        self.foldx_runfile = op.join(self.foldx_dir, 'runfile_FoldX.txt')

    def __call__(self, whatToRun, mutCodes=[]):
        """
        Select which action should be performed by FoldX by setting `whatToRun`.

        Possible values are:

            - AnalyseComplex
            - Stability
            - RepairPDB
            - BuildModel

        See the `FoldX manual`_ for an explanation on what they do.

        .. _FoldX manual: http://foldx.crg.es/manual3.jsp
        """
        logger.debug('Running FoldX {}'.format(whatToRun))
        self.__write_runfile(self.pdb_filename, self.chain_id, whatToRun, mutCodes)
        self.__run_runfile()
        if whatToRun == 'AnalyseComplex':
            return self.__read_result(
                op.join(self.foldx_dir, 'Interaction_AnalyseComplex_resultFile.txt'), whatToRun)
        elif whatToRun == 'Stability':
            return self.__read_result(
                op.join(self.foldx_dir, 'Stability.txt'), whatToRun)
        elif whatToRun == 'RepairPDB':
            return op.join(self.foldx_dir, 'RepairPDB_' + self.pdb_filename)
        elif whatToRun == 'BuildModel':
            # see the FoldX manual for the naming of the generated structures
            if conf.CONFIGS['foldx_num_of_runs'] == 1:
                mutants = [op.join(self.foldx_dir, self.pdb_filename[:-4] + '_1.pdb'), ]
                wiltype = [op.join(self.foldx_dir, 'WT_' + self.pdb_filename[:-4] + '_1.pdb'), ]
                results = [wiltype, mutants]
            else:
                mutants = [
                    op.join(self.foldx_dir, self.pdb_filename[:-4] + '_1_' + str(x) + '.pdb')
                    for x in range(0, conf.CONFIGS['foldx_num_of_runs'])
                ]
                wiltype = [
                    op.join(
                        self.foldx_dir,
                        'WT_' + self.pdb_filename[:-4] + '_1_' + str(x) + '.pdb')
                    for x in range(0, conf.CONFIGS['foldx_num_of_runs'])
                ]
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
            with open(op.join(self.foldx_dir, file_with_mutations), 'w') as fh:
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
                buildModel_runs=conf.CONFIGS['foldx_num_of_runs'],
                water=conf.CONFIGS['foldx_water'],
                output_pdb=output_pdb)

        # This just makes copies of the runfiles for debugging...
        with open(self.foldx_runfile, 'w') as f:
            f.write(foldX_runfile)
        shutil.copy(self.foldx_runfile, op.join(self.foldx_dir, copy_filename))

    def __run_runfile(self):
        """.

        .. todo:: Add a fallback plan using libfaketime.
        """
        # system_command = './FoldX.linux64 -runfile ' + self.foldx_runfile
        system_command = "foldx -runfile '{}'".format(self.foldx_runfile)
        logger.debug("FoldX system command: '{}'".format(system_command))
        env = os.environ.copy()
        env['LD_PRELOAD'] = faketime.config.libfaketime_so_file
        env['FAKETIME'] = '2015-12-26 00:00:00'
        p = helper.run(system_command, cwd=self.foldx_dir, env=env)
        if p.stderr.strip():
            logger.debug('foldx result:\n{}'.format(p.stdout.strip()))
            logger.error('foldx error message:\n{}'.format(p.stderr.strip()))
            if 'Cannot allocate memory' in p.stderr:
                raise errors.ResourceError(p.stderr)
        if 'There was a problem' in p.stdout:
            logger.error('foldx result:\n{}'.format(p.stdout.strip()))
            if 'Specified residue not found.' in p.stdout:
                raise errors.MutationMismatchError()

    def __read_result(self, outFile, whatToRead):
        with open(outFile, 'r') as f:
            lines = f.readlines()
            line = lines[-1].split('\t')
            if whatToRead == 'BuildModel':
                total_energy_difference = line[1]
                return total_energy_difference
            if whatToRead == 'Stability':
                stability_values = [
                    line[x[1]].strip() for x in names_rows_stability
                ]
                return stability_values
            if whatToRead == 'AnalyseComplex':
                complex_stability_values = [
                    line[x[1]].strip() for x in names_rows_stability_complex
                ]
                return complex_stability_values
