import functools
import os.path as op

import pandas as pd

from elaspic.tools._abc import ToolError
from kmtools import py_tools, structure_tools, system_tools

logger = py_tools.get_logger()


class StrideError(ToolError):
    pass


class StrideAnalyser:

    _result_slots = ['secondary_structure']

    def __init__(self, structure_file):
        super().__init__()

    def build(self):
        self.result['secondary_structure'] = self._build()

    @functools.lru_cache(maxsize=512)
    def mutate(self, chain_id, mutation):
        assert self.done

        residue_id = int(mutation[1:-1])
        secondary_structure = (
            self.result['secondary_structure'][
                (self.result['secondary_structure']['chain_id'] == chain_id) &
                (self.result['secondary_structure']['residue_id'] == residue_id)
            ])

        _df = (
            secondary_structure_df[
                (secondary_structure_df.chain == chain_id) &
                (secondary_structure_df.resnum == mutation[1:-1])
            ]
        )
        assert len(secondary_structure_df) == 1
        secondary_structure_df = secondary_structure_df.iloc[0]
        self._validate_mutation(seasa_info['res_name'], mutation)
        secondary_structure = secondary_structure_df.ss_code


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
        p = system_tools.run(system_command, cwd=self.working_dir)
        p.check_returncode()
        return stride_results_file

    def _parse_output(self):
        ...

    def _generate_df(sef):
        ...

    def read_stride_results(stride_results_file):
        def parse_row(row):
            return [
                structure_tools.AAA_DICT[row.split()[1]],
                row.split()[2],
                row.split()[3],
                int(row.split()[4]),
                row.split()[5]
            ]

        with open(stride_results_file) as fh:
            file_data_df = pd.DataFrame(
                [parse_row(row) for row in fh if row[:3] == 'ASG'],
                columns=['amino_acid', 'chain', 'resnum', 'idx', 'ss_code'])
        return file_data_df
