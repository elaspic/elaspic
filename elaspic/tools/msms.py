import functools
import os.path as op
from collections import namedtuple

import pandas as pd

from elaspic.tools._abc import ToolError, StructureAnalyzer
from kmtools import py_tools, structure_tools, system_tools

logger = py_tools.get_logger(__name__)


class MSMSError(ToolError):
    pass


class MSMS(StructureAnalyzer):

    _result_slots = [
        'structure_seasa_by_chain', 'structure_seasa_by_residue',
        'chain_seasa_by_chain', 'chain_seasa_by_residue',
    ]

    def build(self):
        # Calculate for the entire structure
        structure = self.structure
        structure_file = op.join(self.tempdir, structure.id + '.pdb')
        structure_tools.save_structure(structure, structure_file)
        seasa_by_chain, seasa_by_residue = self._build(structure_file)
        self.result.update({
            'structure_seasa_by_chain': seasa_by_chain,
            'structure_seasa_by_residue': seasa_by_residue,
        })
        # Calculate chain by chain
        chain_seasa_by_chain = []
        chain_seasa_by_residue = []
        for chain in structure[0]:
            chain_structure = structure[0].extract([chain.id])
            chain_structure_file = op.join(self.tempdir, structure.id + chain.id + '.pdb')
            structure_tools.save_structure(chain_structure, chain_structure_file)
            seasa_by_chain, seasa_by_residue = self._build(chain_structure_file)
            chain_seasa_by_chain.append(seasa_by_chain)
            chain_seasa_by_residue.append(seasa_by_residue)
        self.result.update({
            'chain_seasa_by_chain': pd.concat(chain_seasa_by_chain, ignore_index=True),
            'chain_seasa_by_residue': pd.concat(chain_seasa_by_residue, ignore_index=True),
        })

    @functools.lru_cache(maxsize=512)
    def analyze(self, chain_id, residue_id, aa):
        assert self.done
        if isinstance(residue_id, int):
            residue_id = (' ', residue_id, ' ')

        structure_seasa = (
            self.result['structure_seasa_by_residue'][
                (self.result['structure_seasa_by_residue']['chain_id'] == chain_id) &
                (self.result['structure_seasa_by_residue']['residue_id'].astype(int) ==
                    residue_id[1])
            ])
        assert len(structure_seasa) == 1, structure_seasa
        structure_seasa = structure_seasa.iloc[0]
        assert structure_tools.AAA_DICT.get(
            structure_seasa['res_name'], structure_seasa['res_name']) == aa

        chain_seasa = (
            self.result['chain_seasa_by_residue'][
                (self.result['chain_seasa_by_residue']['chain_id'] == chain_id) &
                (self.result['chain_seasa_by_residue']['residue_id'].astype(int) ==
                    residue_id[1])
            ])
        assert len(chain_seasa) == 1, chain_seasa
        chain_seasa = chain_seasa.iloc[0]
        assert structure_tools.AAA_DICT.get(
            chain_seasa['res_name'], chain_seasa['res_name']) == aa

        return {
            'solvent_accessibility': structure_seasa['rel_sasa'],
            'solvent_occlusion': structure_seasa['rel_sasa'] - chain_seasa['rel_sasa'],
        }

    # Helper methods
    def _build(self, structure_file):
        xyzrn_file = self._run_pdb_to_xyzrn(structure_file)
        area_file = self._run_msms(xyzrn_file)
        file_data = self._parse_area_file(area_file)
        seasa_by_chain, seasa_by_residue = self._generate_df(file_data)
        return seasa_by_chain, seasa_by_residue

    def _run_pdb_to_xyzrn(self, structure_file):
        """Generate *.xyzrn file."""
        system_command = "pdb_to_xyzrn '{}'".format(structure_file)
        p = system_tools.run(system_command, cwd=self.tempdir)
        p.check_returncode()
        xyzrn_file = op.join(self.tempdir, op.splitext(structure_file)[0] + '.xyzrn')
        with open(xyzrn_file, 'wt') as ofh:
            ofh.write(p.stdout)
        return xyzrn_file

    def _run_msms(self, xyzrn_file):
        """.
        In the future, could add an option to measure residue depth
        using Bio.PDB.ResidueDepth().residue_depth()...
        """
        # Calculate solvent accessible (SASA) and solvent excluded (SESA) surface area
        probe_radius = 1.4
        area_file = op.join(self.tempdir, op.splitext(xyzrn_file)[0] + '.area')
        system_command = (
            "msms "
            " -probe_radius {probe_radius:.1f} "
            " -surface ases "
            " -if '{input_file}' "
            " -af '{output_file}'"
            .format(
                probe_radius=probe_radius,
                input_file=xyzrn_file,
                output_file=area_file))
        p = system_tools.run(system_command, cwd=self.tempdir)
        number_of_tries = 0
        while p.returncode != 0 and number_of_tries < 5:
            logger.warning('MSMS exited with an error!')
            probe_radius -= 0.1
            logger.debug('Reducing probe radius to {}', probe_radius)
            p = system_tools.run(system_command, cwd=self.tempdir)
            number_of_tries += 1
        p.check_returncode()
        return area_file

    def _parse_area_file(self, area_file):
        """Read and parse MSMS output."""
        msms_columns = [
            'chain_id', 'residue_id', 'res_name', 'atom_id', 'atom_num',
            'abs_sesa', 'abs_sasa', 'rel_sasa']
        MSMSRow = namedtuple("MSMSRow", msms_columns)

        def absolute_to_relative(res_name, abs_sasa):
            try:
                ref_sasa = structure_tools.STANDARD_SASA[res_name]
            except KeyError:
                # logger.warning("No reference SASA value for residue {}", res_name)
                return 0.0
            else:
                return abs_sasa / ref_sasa

        def parse_row(row):
            atom_num = int(row[0]) + 1  # atom_num needs to be incremented by 1 for some reason?
            abs_sesa = float(row[1])
            abs_sasa = float(row[2])

            _extra_fields = row[3].split('_')
            atom_id = _extra_fields[0].strip()
            res_name = _extra_fields[1].strip()
            residue_id = _extra_fields[2]
            chain_id = _extra_fields[3]

            rel_sasa = absolute_to_relative(res_name, abs_sasa)

            row = MSMSRow(
                atom_num=atom_num,
                abs_sesa=abs_sesa,
                abs_sasa=abs_sasa,
                atom_id=atom_id,
                res_name=res_name,
                residue_id=residue_id,
                chain_id=chain_id,
                rel_sasa=rel_sasa,
            )
            return row

        with open(area_file, 'r') as fh:
            file_data = fh.readlines()
        file_data = [[l.strip() for l in line.split()] for line in file_data[1:]]
        file_data = [parse_row(row) for row in file_data if row]
        return file_data

    def _generate_df(self, file_data):
        df = pd.DataFrame(data=file_data)
        df_by_chain = df.groupby(['chain_id']).sum().reset_index()
        df_by_residue = df.groupby(['chain_id', 'res_name', 'residue_id']).sum().reset_index()
        return df_by_chain, df_by_residue
