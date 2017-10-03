import logging
import os
import os.path as op
import shlex
import shutil
import subprocess

import faketime.config
import pandas as pd

from elaspic import conf, errors, helper

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
    [name + '_wt' for name in list(zip(*names_rows_stability))[0][:-1]] + ['number_of_residues'])
names_stability_mut = (
    [name + '_mut' for name in list(zip(*names_rows_stability))[0][:-1]] + ['number_of_residues'])

names_rows_stability_complex = ([
    ['intraclashes_energy_1', 3],
    ['intraclashes_energy_2', 4],
] + [[x[0], x[1] + 4] for x in names_rows_stability])
names_stability_complex_wt = (
    [name + '_wt'
     for name in list(zip(*names_rows_stability_complex))[0][:-1]] + ['number_of_residues'])
names_stability_complex_mut = (
    [name + '_mut'
     for name in list(zip(*names_rows_stability_complex))[0][:-1]] + ['number_of_residues'])


class FoldX:
    """FoldX

    Examples
    --------
    >>> import os
    >>> import tempfile
    >>> from elaspic.structure_tools import download_pdb_file
    >>> tmp_dir = tempfile.mkdtemp()
    >>> pdb_file = download_pdb_file('3zml', tmp_dir)
    >>> foldx = FoldX(pdb_file, 'A', tmp_dir)
    >>> structure_file_wt, structure_file_mut, stability_values_wt, stability_values_mut = \
            foldx.build_model('QA93A')

    # >>> stability_values_wt_2 = foldx.stability(structure_file_wt)
    # >>> stability_values_mut_2 = foldx.stability(structure_file_mut)
    # >>> assert stability_values_wt == stability_values_wt_2
    # >>> assert stability_values_mut == stability_values_mut_2
    # >>> foldx.analyze_complex(structure_file_wt)
    # >>> foldx.analyze_complex(structure_file_mut)
    """

    def __init__(self, pdb_file, chain_id, foldx_dir=None):
        """Initialize class for calling FoldX."""
        pdb_file = op.abspath(pdb_file)
        self._tempdir = op.abspath(foldx_dir or conf.CONFIGS['foldx_dir'])
        self._foldx_rotabase = self._find_rotabase()
        try:
            pdb_file = shutil.copy(pdb_file, op.join(self._tempdir, op.basename(pdb_file)))
        except shutil.SameFileError:
            pass
        self.structure_file = self._repair_pdb(op.abspath(pdb_file))
        self.chain_id = chain_id

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
        system_command = (f"foldx --rotabaseLocation {self._foldx_rotabase} --command=RepairPDB "
                          f"--pdb='{op.basename(structure_file)}'")
        self._run(system_command, op.dirname(structure_file))

        # Read results
        repaired_structure_file = shutil.move(
            op.splitext(structure_file)[0] + '_Repair.pdb',
            op.splitext(structure_file)[0] + '-foldx.pdb')
        return repaired_structure_file

    def build_model(self, foldx_mutation):
        """Run FoldX ``BuildModel``."""
        pdb_id = op.basename(op.splitext(self.structure_file)[0])
        cwd = op.dirname(self.structure_file)
        mutation_file = self._get_mutation_file(foldx_mutation, cwd)

        # Run FoldX
        system_command = (f"foldx --rotabaseLocation {self._foldx_rotabase} --command=BuildModel "
                          f"--pdb='{op.basename(self.structure_file)}' "
                          f"--mutant-file='{mutation_file}'")
        self._run(system_command, cwd)

        # Copy FoldX results
        wt_pdb_id = f'WT_{pdb_id}_1.pdb'
        mut_pdb_id = f'{pdb_id}_1.pdb'
        structure_file_wt = shutil.move(
            op.join(cwd, wt_pdb_id), op.join(cwd, f'{pdb_id}-{foldx_mutation}-wt.pdb'))
        structure_file_mut = shutil.move(
            op.join(cwd, mut_pdb_id), op.join(cwd, f'{pdb_id}-{foldx_mutation}-mut.pdb'))

        # Read results
        output_file = op.join(cwd, f'Raw_{pdb_id}.fxout')
        df = pd.read_csv(output_file, sep='\t', skiprows=8)
        os.remove(output_file)

        # Format dataframes
        df.columns = df.rename(columns=str.lower)
        df_wt = df.loc[df['pdb'] == wt_pdb_id, :].drop('pdb', axis=1)
        df_mut = df.loc[df['pdb'] == mut_pdb_id, :].drop('pdb', axis=1)
        assert df_wt.shape[0] == 1 and df_mut.shape[0] == 1

        # Compile results
        stability_values_wt = df_wt.iloc[0].tolist()
        stability_values_mut = df_mut.iloc[0].tolist()
        return structure_file_wt, structure_file_mut, stability_values_wt, stability_values_mut

    def stability(self, structure_file) -> dict:
        """Run FoldX ``Stability``.

        .. deprecated:: `FoldX._build_model` already gives you the same information.
        """
        pdb_id = op.basename(op.splitext(structure_file)[0])
        cwd = op.dirname(structure_file)

        # Run FoldX
        system_command = (f"foldx --rotabaseLocation {self._foldx_rotabase} --command=Stability "
                          f"--pdb='{op.basename(structure_file)}'")
        self._run(system_command, cwd=cwd)

        # Read results
        output_file = op.join(cwd, f'{pdb_id}_0_ST.fxout')
        df = pd.read_csv(output_file, sep='\t', skiprows=8)
        os.remove(output_file)

        # Format dataframe
        df.columns = df.rename(columns=str.lower)
        df.drop('pdb', axis=1, inplace=True)
        assert df.shape[0] == 1
        return df.iloc[0].to_dict()

    def analyse_complex(self, structure_file, chain_ids):
        """Run FoldX ``AnalyseComplex``."""
        pdb_id = op.basename(op.splitext(structure_file)[0])
        chain_id_1, chain_id_2 = chain_ids

        # Run FoldX
        system_command = (
            f"foldx --rotabaseLocation {self._foldx_rotabase} --command=AnalyseComplex "
            f"--pdb='{op.basename(structure_file)}' "
            f"--analyseComplexChains={chain_id_1},{chain_id_2}")
        self._run(system_command, cwd=op.dirname(structure_file))

        # Read results
        output_file = op.join(self._tempdir, f'Interaction_{pdb_id}_AC.fxout')
        df = pd.read_csv(output_file, sep='\t', skiprows=8)
        os.remove(output_file)

        # Format dataframe
        df.columns = df.rename(columns=str.lower)
        df.drop('pdb', axis=1, inplace=True)
        assert df.shape[0] == 1
        return df.iloc[0].to_dict()

    def _get_mutation_file(self, foldx_mutation, cwd) -> str:
        """
        Parameters
        ----------
        foldx_mutation:
            Mutation specified in the following format:
            {mutation.residue_wt}{chain_id}{residue_id}{mutation.residue_mut}
        """
        mutation_file = op.join(cwd, f'individual_list_{foldx_mutation}.txt')
        with open(mutation_file, 'wt') as fout:
            fout.write(f'{foldx_mutation};\n')
        return mutation_file
