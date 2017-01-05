from elaspic.tools._abc import ToolError


class MSMSError(ToolError):
    pass


def __init__(self):
    # Solvent accessibility
    (seasa_by_chain_together, seasa_by_chain_separately,
     seasa_by_residue_together, seasa_by_residue_separately) = self.get_seasa()
    seasa_info = (
        seasa_by_residue_separately[
            (seasa_by_residue_separately['pdb_chain'] == chain_id) &
            (seasa_by_residue_separately['res_num'] == mutation[1:-1])
        ].iloc[0]
    )
    self._validate_mutation(seasa_info['res_name'], mutation)
    solvent_accessibility = seasa_info['rel_sasa']

def get_seasa(self):
    structure_file = self.get_structure_file(''.join(self.chain_ids))
    seasa_by_chain, seasa_by_residue = self._run_msms(structure_file)
    if len(self.chain_ids) > 1:
        seasa_by_chain_separately = []
        seasa_by_residue_separately = []
        for chain_id in self.chain_ids:
            structure_file = self.get_structure_file(chain_id)
            seasa_by_chain, seasa_by_residue = self._run_msms(structure_file)
            seasa_by_chain_separately.append(seasa_by_chain)
            seasa_by_residue_separately.append(seasa_by_residue)
        seasa_by_chain_separately = pd.concat(seasa_by_chain_separately, ignore_index=True)
        seasa_by_residue_separately = pd.concat(seasa_by_residue_separately, ignore_index=True)
        return [
            seasa_by_chain, seasa_by_chain_separately, seasa_by_residue,
            seasa_by_residue_separately
        ]
    else:
        return [None, seasa_by_chain, None, seasa_by_residue]

def _run_msms(self, filename):
    """.

    In the future, could add an option to measure residue depth
    using Bio.PDB.ResidueDepth().residue_depth()...
    """
    base_filename = op.splitext(filename)[0]

    # Convert pdb to xyz coordiates
    assert(os.path.isfile(op.join(self.working_dir, filename)))

    system_command = 'pdb_to_xyzrn {0}.pdb'.format(op.join(self.working_dir, base_filename))
    logger.debug('msms system command 1: %s' % system_command)
    p = system_tools.run(system_command, cwd=self.working_dir)
    if p.returncode != 0:
        logger.debug('msms 1 stdout:\n{}'.format(p.stdout))
        logger.debug('msms 1 stderr:\n{}'.format(p.stderr))
        logger.debug('msms 1 returncode:\n{}'.format(p.returncode))
        raise exc.MSMSError(p.stderr)
    else:
        tempfile_xyzrn = tempfile.NamedTemporaryFile('wt', delete=False)
        tempfile_xyzrn.write(p.stdout)
        tempfile_xyzrn.close()

    # Calculate SASA and SESA (excluded)
    probe_radius = 1.4
    system_command_template = """\
msms -probe_radius {probe_radius:.1f} -surface ases -if '{input_file}' -af '{area_file}' \
"""
    system_command = system_command_template.format(
        probe_radius=probe_radius,
        input_file=tempfile_xyzrn.name,
        area_file=op.join(self.working_dir, base_filename + '.area'))
    logger.debug('msms system command 2: %s' % system_command)
    p = system_tools.run(system_command, cwd=self.working_dir)
    number_of_tries = 0
    while p.returncode != 0 and number_of_tries < 5:
        logger.warning('MSMS exited with an error!')
        probe_radius -= 0.1
        logger.debug('Reducing probe radius to {}'.format(probe_radius))
        p = system_tools.run(system_command, cwd=self.working_dir)
        number_of_tries += 1
    if p.returncode != 0:
        logger.debug('msms stdout 2:\n{}'.format(p.stdout))
        logger.debug('msms stderr 2:\n{}'.format(p.stderr))
        logger.debug('msms returncode 2:\n{}'.format(p.returncode))
        raise exc.MSMSError(p.stderr)
    os.remove(tempfile_xyzrn.name)

    # Read and parse the output
    with open(op.join(self.working_dir, base_filename + '.area'), 'r') as fh:
        file_data = fh.readlines()
    file_data = [
        [l.strip() for l in line.split()] for line in file_data
    ]
    del file_data[0]

    msms_columns = [
        'atom_num', 'abs_sesa', 'abs_sasa', 'atom_id', 'res_name', 'res_num', 'pdb_chain'
    ]

    def msms_parse_row(row):
        parsed_row = [
            int(row[0]), float(row[1]), float(row[2]),
            row[3].split('_')[0].strip(),
            row[3].split('_')[1].strip(),
            row[3].split('_')[2],
            row[3].split('_')[3]
        ]
        return parsed_row

    file_data = [msms_parse_row(row) for row in file_data if row]
    seasa_df = pd.DataFrame(data=file_data, columns=msms_columns)
    seasa_df['atom_num'] = seasa_df['atom_num'].apply(lambda x: x + 1)
    seasa_df['rel_sasa'] = [
        x[0] / STANDARD_SASA.get(x[1], x[0]) * 100
        for x in zip(seasa_df['abs_sasa'], seasa_df['res_name'])
    ]

    seasa_gp_by_chain = seasa_df.groupby(['pdb_chain'])
    seasa_gp_by_residue = seasa_df.groupby(['pdb_chain', 'res_name', 'res_num'])
    seasa_by_chain = seasa_gp_by_chain.sum().reset_index()
    seasa_by_residue = seasa_gp_by_residue.sum().reset_index()

    return seasa_by_chain, seasa_by_residue
