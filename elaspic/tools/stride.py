from elaspic.tools._abc import ToolError


class StrideError(ToolError):
    pass


def __init__(self):
    # Secondary structure
    secondary_structure_df = self.get_secondary_structure()
    secondary_structure_df = (
        secondary_structure_df[
            (secondary_structure_df.chain == chain_id) &
            (secondary_structure_df.resnum == mutation[1:-1])
        ]
    )
    assert len(secondary_structure_df) == 1
    secondary_structure_df = secondary_structure_df.iloc[0]
    self._validate_mutation(seasa_info['res_name'], mutation)
    secondary_structure = secondary_structure_df.ss_code

def get_secondary_structure(self):
    """Run `stride` to calculate protein secondary structure."""
    structure_file = self.get_structure_file(''.join(self.chain_ids))
    stride_results_file = op.join(
        op.dirname(structure_file),
        structure_tools.get_pdb_id(structure_file) + '_stride_results.txt'
    )
    system_command = 'stride {} -f{}'.format(structure_file, stride_results_file)
    logger.debug('stride system command: %s' % system_command)
    p = system_tools.run(system_command, cwd=self.working_dir)
    logger.debug('stride return code: %i' % p.returncode)
    logger.debug('stride result: %s' % p.stdout)
    logger.debug('stride error: %s' % p.stderr)
    # collect results
    with open(stride_results_file) as fh:
        file_data_df = pd.DataFrame(
            [[AAA_DICT[row.split()[1]], row.split()[2],
              row.split()[3], int(row.split()[4]), row.split()[5]]
             for row in fh.readlines() if row[:3] == 'ASG'],
            columns=['amino_acid', 'chain', 'resnum', 'idx', 'ss_code'])
    return file_data_df
