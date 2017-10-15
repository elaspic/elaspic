import logging
import os
import os.path as op
import shlex
import shutil
import subprocess

import faketime.config
import pandas as pd

from elaspic import conf, errors

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
names_stability = list(next(zip(*names_rows_stability)))
names_stability_wt = [name + '_wt' for name in names_stability[:-1]] + ['number_of_residues']
names_stability_mut = [name + '_mut' for name in names_stability[:-1]] + ['number_of_residues']

names_rows_stability_complex = ([
    ['intraclashes_energy_1', 3],
    ['intraclashes_energy_2', 4],
] + [[x[0], x[1] + 4] for x in names_rows_stability])
names_stability_complex = list(next(zip(*names_rows_stability_complex)))
names_stability_complex_wt = [name + '_wt'
                              for name in names_stability_complex[:-1]] + ['number_of_residues']
names_stability_complex_mut = [name + '_mut'
                               for name in names_stability_complex[:-1]] + ['number_of_residues']


def read_build_model(output_file, wt_pdb_id, mut_pdb_id):
    df = pd.read_csv(output_file, sep='\t', skiprows=8)
    # Format dataframes
    df = df.rename(columns=str.lower)
    logger.debug(df.head())
    logger.info(df.head())
    df_wt = df.loc[df['pdb'] == wt_pdb_id, :].drop('pdb', axis=1)
    df_mut = df.loc[df['pdb'] == mut_pdb_id, :].drop('pdb', axis=1)
    assert df_wt.shape[0] == 1 and df_mut.shape[0] == 1
    # Compile results
    stability_values_wt = df_wt.iloc[0].tolist()
    stability_values_mut = df_mut.iloc[0].tolist()
    return stability_values_wt, stability_values_mut


def read_stability(output_file):
    df = pd.read_csv(output_file, sep='\t', names=['pdb'] + names_stability, index_col=False)
    # Format dataframe
    df = df.rename(columns=str.lower)
    logger.debug(df.head())
    assert df.shape[0] == 1
    result = df.drop('pdb', axis=1).iloc[0].tolist()
    return result


def read_analyse_complex(output_file):
    df = pd.read_csv(output_file, sep='\t', index_col=False, skiprows=8)
    # Format dataframe
    df = df.rename(columns=lambda s: s.lower().replace(' ', '_'))
    logger.debug(df.head())
    assert df.shape[0] == 1
    result = df.drop(pd.Index(['pdb', 'group1', 'group2']), axis=1).iloc[0].tolist()
    return result


class FoldX:

    def __init__(self, foldx_dir=None):
        self._tempdir = op.abspath(foldx_dir or conf.CONFIGS['foldx_dir'])
        logger.debug("FoldX._tempdir: %s", self._tempdir)
        self._foldx_rotabase = self._find_rotabase()
        logger.debug("FoldX._foldx_rotabase: %s", self._foldx_rotabase)

    def _find_rotabase(self):
        system_command = 'which rotabase.txt'
        process = subprocess.run(
            shlex.split(system_command),
            stdout=subprocess.PIPE,
            universal_newlines=True,
            check=True)
        return process.stdout.strip()

    def _run(self, system_command, cwd):
        env = os.environ.copy()
        env['LD_PRELOAD'] = faketime.config.libfaketime_so_file
        env['FAKETIME'] = '2015-12-26 00:00:00'
        logger.debug(system_command)
        process = subprocess.run(
            shlex.split(system_command),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            cwd=cwd,
            env=env)
        if process.stderr.strip():
            logger.debug('foldx result:\n{}'.format(process.stdout.strip()))
            logger.error('foldx error message:\n{}'.format(process.stderr.strip()))
            if 'Cannot allocate memory' in process.stderr:
                raise errors.ResourceError(process.stderr.strip())
        if 'There was a problem' in process.stdout:
            logger.error('foldx result:\n{}'.format(process.stdout.strip()))
            if 'Specified residue not found.' in process.stdout:
                raise errors.MutationMismatchError()
        if process.returncode != 0:
            logger.error("Command finished with return code %s", process.returncode)
            logger.error(process.stdout.strip())
            logger.error(process.stderr.strip())
            raise subprocess.CalledProcessError(
                returncode=process.returncode,
                cmd=system_command,
                output=process.stdout.strip(),
                stderr=process.stderr.strip())

    def _repair_pdb(self, structure_file):
        """Run FoldX ``RepairPDB``."""
        # Run FoldX
        system_command = (
            "foldx --rotabaseLocation {} --command=RepairPDB ".format(self._foldx_rotabase) +
            "--pdb='{}'".format(op.basename(structure_file)))
        self._run(system_command, op.dirname(structure_file))

        # Read results
        repaired_structure_file = shutil.move(
            op.splitext(structure_file)[0] + '_Repair.pdb',
            op.splitext(structure_file)[0] + '-foldx.pdb')
        return repaired_structure_file

    def build_model(self, pdb_file, foldx_mutation):
        """Run FoldX ``BuildModel``.

        .. note::

            For some reason, the results of ``BuildModel``
            do not include ``number_of_residues``.
        """
        pdb_file = op.abspath(pdb_file)
        try:
            pdb_file = shutil.copy(pdb_file, op.join(self._tempdir, op.basename(pdb_file)))
        except shutil.SameFileError:
            pass
        structure_file = self._repair_pdb(op.abspath(pdb_file))

        pdb_id = op.basename(op.splitext(structure_file)[0])
        cwd = op.dirname(structure_file)
        mutation_file = self._get_mutation_file(foldx_mutation, cwd)

        # Run FoldX
        system_command = (
            "foldx --rotabaseLocation {} --command=BuildModel ".format(self._foldx_rotabase) +
            "--pdb='{}' ".format(op.basename(structure_file)) +
            "--mutant-file='{}'".format(mutation_file))
        self._run(system_command, cwd)

        # Copy FoldX results
        wt_pdb_id = 'WT_{}_1.pdb'.format(pdb_id)
        mut_pdb_id = '{}_1.pdb'.format(pdb_id)
        structure_file_wt = shutil.move(
            op.join(cwd, wt_pdb_id), op.join(cwd, '{}-{}-wt.pdb'.format(pdb_id, foldx_mutation)))
        structure_file_mut = shutil.move(
            op.join(cwd, mut_pdb_id), op.join(cwd, '{}-{}-mut.pdb'.format(pdb_id, foldx_mutation)))

        # Read results
        output_file = op.join(cwd, 'Raw_{}.fxout'.format(pdb_id))
        stability_values_wt, stability_values_mut = read_build_model(output_file, wt_pdb_id,
                                                                     mut_pdb_id)
        return structure_file_wt, structure_file_mut, stability_values_wt, stability_values_mut

    def stability(self, structure_file) -> dict:
        """Run FoldX ``Stability``.

        .. deprecated:: `FoldX._build_model` already gives you the same information.
        """
        pdb_id = op.basename(op.splitext(structure_file)[0])
        cwd = op.dirname(structure_file)

        # Run FoldX
        system_command = (
            "foldx --rotabaseLocation {} --command=Stability ".format(self._foldx_rotabase) +
            "--pdb='{}'".format(op.basename(structure_file)))
        self._run(system_command, cwd=cwd)

        # Read results
        output_file = op.join(cwd, '{}_0_ST.fxout'.format(pdb_id))
        result = read_stability(output_file)
        return result

    def analyse_complex(self, structure_file, chain_ids):
        """Run FoldX ``AnalyseComplex``."""
        pdb_id = op.basename(op.splitext(structure_file)[0])
        chain_id_1, chain_id_2 = chain_ids

        # Run FoldX
        system_command = (
            "foldx --rotabaseLocation {} --command=AnalyseComplex ".format(self._foldx_rotabase) +
            "--pdb='{}' ".format(op.basename(structure_file)) +
            "--analyseComplexChains={},{}".format(chain_id_1, chain_id_2))
        self._run(system_command, cwd=op.dirname(structure_file))

        # Read results
        output_file = op.join(self._tempdir, 'Interaction_{}_AC.fxout'.format(pdb_id))
        result = read_analyse_complex(output_file)
        return result

    def _get_mutation_file(self, foldx_mutation, cwd) -> str:
        """
        Parameters
        ----------
        foldx_mutation:
            Mutation specified in the following format:
            {mutation.residue_wt}{chain_id}{residue_id}{mutation.residue_mut}
        """
        mutation_file = op.join(cwd, 'individual_list_{}.txt'.format(foldx_mutation))
        with open(mutation_file, 'wt') as fout:
            fout.write('{};\n'.format(foldx_mutation))
        return mutation_file
