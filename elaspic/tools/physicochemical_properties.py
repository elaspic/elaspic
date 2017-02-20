        # # PhysiChem
        # physchem, physchem_ownchain = self.get_physi_chem(chain_id, mutation)


def get_physi_chem(self, chain_id, mutation):
    """Return the atomic contact vector.

    Count how many interactions there are between charged, polar or "carbon" residues.
    The "carbon" interactions give you information about the Van der Waals packing
    of the residues. Comparing the wildtype vs. the mutant values is used in
    the machine learning algorithm.

    'mutation' is of the form: 'A16' where A is the chain identifier and 16
    the residue number (in pdb numbering) of the mutation
    chainIDs is a list of strings with the chain identifiers to be used
    if more than two chains are given, the chains not containing the mutation
    are considered as "opposing" chain
    """
    model = self.sp.structure[0]
    mutated_chain = model[chain_id]
    opposite_chains = [chain for chain in model.child_list if chain.id != chain_id]
    main_chain_atoms = ['CA', 'C', 'N', 'O']

    # Find the mutated residue (assuming mutation is in resnum)
    for residue in mutated_chain:
        if residue.resname in AMINO_ACIDS and not residue.id[0].strip():
            if str(residue.id[1]) == mutation[1:-1]:
                self._validate_mutation(residue.resname, mutation)
                mutated_residue = residue
                break
    mutated_atoms = [atom for atom in mutated_residue if atom.name not in main_chain_atoms]

    # Go through each atom in each residue in each partner chain...
    opposite_chain_contacts = {
        'equal_charge': [],
        'opposite_charge': [],
        'h_bond': [],
        'carbon_contact': []}
    for opposite_chain in opposite_chains:
        for partner_residue in opposite_chain:
            self._increment_vector(
                mutated_residue, mutated_atoms, partner_residue, opposite_chain_contacts)

    # Go through each atom in each residue in the mutated chain...
    same_chain_contacts = {
        'equal_charge': [],
        'opposite_charge': [],
        'h_bond': [],
        'carbon_contact': []}
    for partner_residue in mutated_chain:
        # Skipping the mutated residue...
        if partner_residue == mutated_residue:
            continue
        self._increment_vector(
            mutated_residue, mutated_atoms, partner_residue, same_chain_contacts)

    opposite_chain_contact_vector = [
        len(opposite_chain_contacts['equal_charge']),
        len(opposite_chain_contacts['opposite_charge']),
        len(opposite_chain_contacts['h_bond']),
        len(set(opposite_chain_contacts['carbon_contact']))]

    same_chain_contact_vector = [
        len(same_chain_contacts['equal_charge']),
        len(same_chain_contacts['opposite_charge']),
        len(same_chain_contacts['h_bond']),
        len(set(same_chain_contacts['carbon_contact']))]

    return opposite_chain_contact_vector, same_chain_contact_vector

def _increment_vector(
        self, mutated_residue, mutated_atoms, partner_residue, contact_features_dict):
    # For each residue each atom of the mutated residue has to be checked
    for mutated_atom in mutated_atoms:
        mutated_atom_type = self._get_atom_type(mutated_residue.resname, mutated_atom)
        # And each partner residue and partner atom
        for partner_atom in partner_residue:
            r = structure_tools.calculate_distance(
                mutated_atom, partner_atom, self.vdw_distance)
            if r is not None:
                partner_atom_type = self._get_atom_type(partner_residue.resname, partner_atom)
                if partner_atom_type == 'ignore':
                    continue
                if mutated_atom_type == 'carbon' and partner_atom_type == 'carbon':
                    # The Van der Waals packing should be determined
                    # nonredundant. Thus, the atomic coordinates are
                    # used to keep track of which interactions where
                    # already counted. Does not matter as much for others.
                    contact_features_dict['carbon_contact'].append(
                        tuple(partner_atom.coord))
                if r <= self.min_contact_distance:
                    if (mutated_atom_type == 'charged_plus' and
                            partner_atom_type == 'charged_plus'):
                        contact_features_dict['equal_charge'].append(
                            tuple(mutated_atom.coord))
                    if (mutated_atom_type == 'charged_minus' and
                            partner_atom_type == 'charged_plus'):
                        contact_features_dict['opposite_charge'].append(
                            tuple(mutated_atom.coord))
                    if (mutated_atom_type == 'charged_plus' and
                            partner_atom_type == 'charged_minus'):
                        contact_features_dict['opposite_charge'].append(
                            tuple(mutated_atom.coord))
                    if (mutated_atom_type == 'charged_minus' and
                            partner_atom_type == 'charged_minus'):
                        contact_features_dict['equal_charge'].append(
                            tuple(mutated_atom.coord))
                    if (mutated_atom_type == 'charged' and
                            partner_atom_type == 'polar'):
                        contact_features_dict['h_bond'].append(
                            tuple(mutated_atom.coord))
                    if (mutated_atom_type == 'polar' and
                            partner_atom_type == 'charged'):
                        contact_features_dict['h_bond'].append(
                            tuple(mutated_atom.coord))
                    if (mutated_atom_type == 'polar' and
                            partner_atom_type == 'polar'):
                        contact_features_dict['h_bond'].append(
                            tuple(mutated_atom.coord))

def _get_atom_type(self, residue, atom):
    """Get the type of atom we are dealing with (i.e. charged, polar, carbon).

    In order to see what "interaction" type two atoms are forming, check the
    individual label which every atom in every residue has (see pdb file
    convention for an explanation of the labels).
    With this label, one can determine which atom of the residue one is looking
    at, and hence, one can determine which "interaction" two atoms are forming.
    """
    # This is based on the naming convention for the atoms in crystalography
    charged_plus = ['NH1', 'NH2', 'NZ']
    # Label of negatively charged atoms
    charged_minus = ['OD1', 'OD2', 'OE1', 'OE2']
    # Label of polar atoms
    polar = [
        'OG', 'OG1', 'OD1', 'OD2', 'ND1', 'OE1', 'NE', 'NE1', 'NE2', 'ND1', 'ND2', 'SG',
        'OH', 'O', 'N'
    ]

    if residue.upper() in ['ARG', 'R', 'LYS', 'K']:
        if atom.name in charged_plus:
            return 'charged_plus'

    if residue.upper() in ['ASP', 'D', 'GLU', 'E']:
        if atom.name in charged_minus:
            return 'charged_minus'

    if atom.name in polar:
        return 'polar'

    if atom.name[0] == 'C' or atom.name == 'SD':
        return 'carbon'

    return 'ignore'
