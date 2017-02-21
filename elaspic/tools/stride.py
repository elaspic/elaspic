import functools
import os.path as op

import pandas as pd

from elaspic.tools._abc import StructureAnalyzer, ToolError
from kmtools import py_tools, structure_tools, system_tools

logger = py_tools.get_logger()


class StrideError(ToolError):
    pass


class Stride(StructureAnalyzer):

    _result_slots = ['secondary_structure']

    def build(self):
        if self.done:
            logger.info("Already built!")
            return
        self.structure_file = op.join(self.tempdir, self.structure.id + '.pdb')
        structure_tools.save_structure(self.structure, self.structure_file)
        self.result['secondary_structure'] = self._build()

    @functools.lru_cache(maxsize=512)
    def analyze(self, chain_id, residue_id, aa):
        assert self.done
        if isinstance(residue_id, int):
            residue_id = (' ', residue_id, ' ')

        secondary_structure = (
            self.result['secondary_structure'][
                (self.result['secondary_structure']['chain_id'] == chain_id) &
                (self.result['secondary_structure']['resnum'].astype(int) == residue_id[1])
            ])
        assert len(secondary_structure) == 1
        secondary_structure = secondary_structure.iloc[0]
        assert structure_tools.AAA_DICT.get(
            secondary_structure['amino_acid'], secondary_structure['amino_acid']) == aa
        # secondary_structure = secondary_structure_df.ss_code
        return {'ss_code': secondary_structure['ss_code']}

    # Helper
    def _build(self):
        stride_results_file = self._run_stride()
        file_data = self._parse_output(stride_results_file)
        secondary_structure = self._generate_df(file_data)
        return secondary_structure

    def _run_stride(self):
        """Run `stride` to calculate protein secondary structure."""
        stride_results_file = op.join(
            self.tempdir, op.splitext(self.structure_file)[0] + '.stride')
        system_command = 'stride {} -f{}'.format(self.structure_file, stride_results_file)
        p = system_tools.run(system_command, cwd=self.tempdir)
        p.check_returncode()
        return stride_results_file

    def _parse_output(self, stride_results_file):
        def parse_row(row):
            return [
                structure_tools.AAA_DICT[row.split()[1]],
                row.split()[2],
                row.split()[3],
                int(row.split()[4]),
                row.split()[5]
            ]
        with open(stride_results_file) as fh:
            file_data = [parse_row(row) for row in fh if row[:3] == 'ASG']
        return file_data

    def _generate_df(self, file_data):
        file_data_df = pd.DataFrame(
            file_data, columns=['amino_acid', 'chain_id', 'resnum', 'idx', 'ss_code'])
        return file_data_df
