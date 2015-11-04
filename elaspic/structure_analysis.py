# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import zip
from builtins import object

import os
import os.path as op
import time
import subprocess
import logging
from collections import OrderedDict

import six
import pandas as pd

from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser

from . import errors, helper, structure_tools

logger = logging.getLogger(__name__)


# %% CONSTANTS

# Standard accessibilities for a ALA-X-ALA tripeptide
# (obtained from NACCESS)
STANDARD_DATA = """
STANDARD ACCESSIBILITES FOR PROBE 1.40 AND RADII vdw.radii
ATOM S   2  ALA  107.95   0.0  69.41   0.0   0.00   0.0  69.41   0.0  38.54   0.0  71.38   0.0  36.58   0.0
ATOM S   2  CYS  134.28   0.0  96.75   0.0   0.00   0.0  96.75   0.0  37.53   0.0  97.93   0.0  36.35   0.0
ATOM S   2  ASP  140.39   0.0  48.00   0.0  54.69   0.0 102.69   0.0  37.70   0.0  49.24   0.0  91.15   0.0
ATOM S   2  GLU  172.25   0.0  59.10   0.0  75.64   0.0 134.74   0.0  37.51   0.0  60.29   0.0 111.96   0.0
ATOM S   2  PHE  199.48   0.0 164.11   0.0   0.00   0.0 164.11   0.0  35.37   0.0 165.25   0.0  34.23   0.0
ATOM S   2  GLY   80.10   0.0  32.33   0.0   0.00   0.0  32.33   0.0  47.77   0.0  37.55   0.0  42.55   0.0
ATOM S   2  HIS  182.88   0.0  96.01   0.0  51.07   0.0 147.08   0.0  35.80   0.0  97.15   0.0  85.73   0.0
ATOM S   2  ILE  175.12   0.0 137.96   0.0   0.00   0.0 137.96   0.0  37.16   0.0 139.14   0.0  35.98   0.0
ATOM S   2  LYS  200.81   0.0 115.38   0.0  47.92   0.0 163.30   0.0  37.51   0.0 116.57   0.0  84.24   0.0
ATOM S   2  LEU  178.63   0.0 141.12   0.0   0.00   0.0 141.12   0.0  37.51   0.0 142.31   0.0  36.32   0.0
ATOM S   2  MET  194.15   0.0 156.64   0.0   0.00   0.0 156.64   0.0  37.51   0.0 157.84   0.0  36.32   0.0
ATOM S   2  ASN  143.94   0.0  44.98   0.0  61.26   0.0 106.24   0.0  37.70   0.0  46.23   0.0  97.72   0.0
ATOM S   2  PRO  136.13   0.0 119.90   0.0   0.00   0.0 119.90   0.0  16.23   0.0 120.95   0.0  15.19   0.0
ATOM S   2  GLN  178.50   0.0  51.03   0.0  89.96   0.0 140.99   0.0  37.51   0.0  52.22   0.0 126.28   0.0
ATOM S   2  ARG  238.76   0.0  76.60   0.0 124.65   0.0 201.25   0.0  37.51   0.0  77.80   0.0 160.97   0.0
ATOM S   2  SER  116.50   0.0  46.89   0.0  31.22   0.0  78.11   0.0  38.40   0.0  48.55   0.0  67.95   0.0
ATOM S   2  THR  139.27   0.0  74.54   0.0  27.17   0.0 101.70   0.0  37.57   0.0  75.72   0.0  63.55   0.0
ATOM S   2  VAL  151.44   0.0 114.28   0.0   0.00   0.0 114.28   0.0  37.16   0.0 115.47   0.0  35.97   0.0
ATOM S   2  TRP  249.36   0.0 187.67   0.0  23.60   0.0 211.26   0.0  38.10   0.0 189.67   0.0  59.69   0.0
ATOM S   2  TYR  212.76   0.0 135.35   0.0  42.03   0.0 177.38   0.0  35.38   0.0 136.50   0.0  76.26   0.0
"""
STANDARD_SASA_ALL = [[l.strip() for l in line.split()] for line in STANDARD_DATA.strip().split('\n')[1:]]
STANDARD_SASA = {x[3]: float(x[4]) for x in STANDARD_SASA_ALL}


# %% Init
class AnalyzeStructure(object):
    """
    Runs the program pops to calculate the interface size of the complexes
    This is done by calculating the surface of the complex and the seperated parts.
    The interface is then given by the substracting.
    """

    def __init__(self, pdb_file, working_dir, vdw_distance=5.0, min_contact_distance=4.0):
        self.pdb_file = pdb_file
        #: Folder with all the binaries (i.e. ./analyze_structure)
        self.working_dir = working_dir
        self.vdw_distance = vdw_distance
        self.min_contact_distance = min_contact_distance

        self._prepare_temp_folder(self.working_dir)

        self.sp = structure_tools.StructureParser(pdb_file)
        self.sp.extract()
        self.sp.save_structure(output_dir=self.working_dir)

        self.chain_ids = self.sp.chain_ids


    def _prepare_temp_folder(self, temp_folder):
        # ./analyze_structure
        os.makedirs(temp_folder, exist_ok=True)



    #%%
    def __call__(self, chain_id, mutation, chain_id_other=None):
        """Calculate all properties.
        """
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

        # Contact distance
        contact_distance = None
        if chain_id_other is not None:
            try:
                contact_distance = self.get_interchain_distances(chain_id, mutation)
                contact_distance = contact_distance[chain_id][chain_id_other]
                logger.debug(
                    'The shortest interchain distance between chain {} and chain {} is {}'
                    .format(chain_id, chain_id_other, contact_distance))
                if not contact_distance:
                    raise ValueError()
            except (IndexError, KeyError, ValueError) as e:
                logger.error(
                    'Could not calculate the shortest contact distance between two chains!'
                )
                logger.error(e)
                logger.error(contact_distance)
                raise

        # PhysiChem
        physchem, physchem_ownchain = self.get_physi_chem(chain_id, mutation)

        # Compile results
        results = dict(
            solvent_accessibility = solvent_accessibility,
            secondary_structure = secondary_structure,
            contact_distance = contact_distance,
            physchem = physchem,
            physchem_ownchain = physchem_ownchain,
        )
        return results


    def get_structure_file(self, chains, ext='.pdb'):
        return op.join(self.working_dir, self.sp.pdb_id + chains + ext)



    #%%
    def get_physi_chem(self, chain_id, mutation):
        """
        Return the atomic contact vector, that is, counting how many interactions
        between charged, polar or "carbon" residues there are. The "carbon"
        interactions give you information about the Van der Waals packing of
        the residues. Comparing the wildtype vs. the mutant values is used in
        the machine learning algorithm.

        'mutation' is of the form: 'A16' where A is the chain identifier and 16
        the residue number (in pdb numbering) of the mutation
        chainIDs is a list of strings with the chain identifiers to be used
        if more than two chains are given, the chains not containing the mutation
        are considered as "opposing" chain
        """
        model = self.sp.structure[0]
        mutated_chain = model[chain_id]
        opposite_chains = [ chain for chain in model.child_list if chain.id != chain_id ]
        main_chain_atoms = ['CA', 'C', 'N', 'O']

        # Find the mutated residue (assuming mutation is in index coordinates)
#        mutation_pos = int(mutation[1:-1])
#        residue_counter = 0
#        for residue in mutated_chain:
#            if residue.resname in structure_tools.AMINO_ACIDS and not residue.id[0].strip():
#                residue_counter += 1
#                if residue_counter == mutation_pos:
#                    self._validate_mutation(residue.resname, mutation)
#                    mutated_residue = residue
#                    break
#        mutated_atoms = [atom for atom in mutated_residue if atom.name not in main_chain_atoms]

        # Find the mutated residue (assuming mutation is in resnum)
        for residue in mutated_chain:
            if residue.resname in structure_tools.AMINO_ACIDS and not residue.id[0].strip():
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


    def _validate_mutation(self, resname, mutation):
        valid_aa = [mutation[0].upper(), mutation[-1].upper()]
        if structure_tools.AAA_DICT.get(resname, resname) not in valid_aa:
            logger.error(resname)
            logger.error(mutation)
            logger.error(structure_tools.AAA_DICT[resname])
            raise errors.MutationMismatchError()


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
                        contact_features_dict['carbon_contact'].append(tuple(partner_atom.coord))
                    if r <= self.min_contact_distance:
                        if mutated_atom_type == 'charged_plus' and partner_atom_type == 'charged_plus':
                            contact_features_dict['equal_charge'].append(tuple(mutated_atom.coord))
                        if mutated_atom_type == 'charged_minus' and partner_atom_type == 'charged_plus':
                            contact_features_dict['opposite_charge'].append(tuple(mutated_atom.coord))
                        if mutated_atom_type == 'charged_plus' and partner_atom_type == 'charged_minus':
                            contact_features_dict['opposite_charge'].append(tuple(mutated_atom.coord))
                        if mutated_atom_type == 'charged_minus' and partner_atom_type == 'charged_minus':
                            contact_features_dict['equal_charge'].append(tuple(mutated_atom.coord))
                        if mutated_atom_type == 'charged' and partner_atom_type == 'polar':
                            contact_features_dict['h_bond'].append(tuple(mutated_atom.coord))
                        if mutated_atom_type == 'polar' and partner_atom_type == 'charged':
                            contact_features_dict['h_bond'].append(tuple(mutated_atom.coord))
                        if mutated_atom_type == 'polar' and partner_atom_type == 'polar':
                            contact_features_dict['h_bond'].append(tuple(mutated_atom.coord))


    def _get_atom_type(self, residue, atom):
        """
        Checks what type of atom it is (i.e. charged, polar, carbon)

        In order to see what "interaction" type two atoms are forming, check the
        individual label which every atom in every residue has (see pdb file
        convention for an explanation of the labels).
        With this label, one can determine which atom of the residue one is looking
        at, and hence, one can determine which "interaction" two atoms are forming.
        """
        # This is based on the naming convention for the atoms in crystalography
        charged_plus  = ['NH1', 'NH2', 'NZ']
        # Label of negatively charged atoms
        charged_minus = ['OD1', 'OD2', 'OE1', 'OE2']
        # Label of polar atoms
        polar = ['OG', 'OG1', 'OD1', 'OD2', 'ND1', 'OE1', 'NE', 'NE1', 'NE2', 'ND1', 'ND2', 'SG', 'OH', 'O', 'N']

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


    #%% SASA New
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
            return [seasa_by_chain, seasa_by_chain_separately, seasa_by_residue, seasa_by_residue_separately]
        else:
            return [None, seasa_by_chain, None, seasa_by_residue]


    def _run_msms(self, filename):
        """ In the future, could add an option to measure residue depth
        using Bio.PDB.ResidueDepth().residue_depth()...
        """
        base_filename = op.splitext(filename)[0]

        # Convert pdb to xyz coordiates
        assert(os.path.isfile(op.join(self.working_dir, filename)))

        system_command = 'pdb_to_xyzrn {0}.pdb'.format(op.join(self.working_dir, base_filename))
        logger.debug('msms system command 1: %s' % system_command)
        result, error_message, return_code = (
            helper.subprocess_check_output_locally(self.working_dir, system_command)
        )
        if return_code != 0:
            logger.debug('msms result 1:')
            logger.debug(result)
            logger.debug('msms error message 1:')
            logger.debug(error_message)
            logger.debug('naccess rc 1:')
            logger.debug(return_code)
            raise errors.MSMSError(error_message)
        else:
            with open(op.join(self.working_dir, base_filename + '.xyzrn'), 'w') as ofh:
                ofh.writelines(result)

        # Calculate SASA and SESA (excluded)
        probe_radius = 1.4
        system_command_string = (
            'msms '
            '-probe_radius {1:.1f} '
            '-surface ases '
            '-if {0}.xyzrn '
            '-af {0}.area')
        system_command = system_command_string.format(op.join(self.working_dir, base_filename), probe_radius)
        logger.debug('msms system command 2: %s' % system_command)
        result, error_message, return_code = (
            helper.subprocess_check_output_locally(self.working_dir, system_command)
        )
        number_of_tries = 0
        while return_code != 0 and number_of_tries < 5:
            logger.error('MSMS exited with an error!')
            probe_radius -= 0.1
            logger.debug('Reducing probe radius to {}'.format(probe_radius))
            system_command = system_command_string.format(op.join(self.working_dir, base_filename), probe_radius)
            result, error_message, return_code = (
                helper.subprocess_check_output_locally(self.working_dir, system_command)
            )
            number_of_tries += 1
        if return_code != 0:
            logger.debug('msms result 2:')
            logger.debug(result)
            logger.debug('msms error message 2:')
            logger.debug(error_message)
            logger.debug('naccess rc 2:')
            logger.debug(return_code)
            raise errors.MSMSError(error_message)

        # Read and parse the output
        with open(op.join(self.working_dir, base_filename + '.area'), 'r') as fh:
            file_data = fh.readlines()
        file_data = [ [l.strip() for l in line.split()] for line in file_data]
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


    #%% SASA Old
    def get_sasa(self, program_to_use='pops'):
        """Get Solvent Accessible Surface Area scores.

        .. note:: deprecated

            Use :fn:`get_seasa` instead.

        """
        if program_to_use == 'naccess':
            run_sasa_atom = self._run_naccess_atom
        elif program_to_use == 'pops':
            run_sasa_atom = self._run_pops_atom
        else:
            raise Exception('Unknown program specified!')
        sasa_score_splitchains = {}
        for chain_id in self.chain_ids:
            sasa_score_splitchains.update(run_sasa_atom(chain_id + '.pdb'))
        sasa_score_allchains = run_sasa_atom(''.join(self.chain_ids) + '.pdb')
        return [sasa_score_splitchains, sasa_score_allchains]


    def _run_naccess_atom(self, filename):
        """
        NACESS is not used anymore.
        It is not reliable enough on a PDB-wide scale.
        """
        # run naccess
        system_command = ('naccess ' + filename)
        logger.debug('naccess system command: %s' % system_command)
        assert(os.path.isfile(self.working_dir + filename))
        result, error_message, return_code = (
            helper.subprocess_check_output_locally(self.working_dir, system_command)
        )
        logger.debug('naccess result: {}'.format(result))
        logger.debug('naccess error: {}'.format(error_message))
        logger.debug('naccess rc: {}'.format(return_code))
        # Collect results
        sasa_scores = {}
        with open(self.working_dir + filename.split('.')[0] + '.rsa') as fh:
            for line in fh:
                row = line.split()
                if row[0] != 'RES':
                    continue
                try:
                    (line_id, res, chain, num, all_abs, all_rel,
                    sidechain_abs, sidechain_rel, mainchain_abs, mainchain_rel,
                    nonpolar_abs, nonpolar_rel, polar_abs, polar_rel) = row
                except ValueError as e:
                    print(e)
                    print(line)
                    print(row)
                sasa_scores.setdefault(chain, []).append(sidechain_rel) # percent sasa on sidechain
        return sasa_scores


    def _run_pops_atom(self, chain_id):
        # Use pops to calculate sasa score for the given chain
        termination, rc, e = self.__run_pops_atom(chain_id)
        if termination != 'Clean termination':
            logger.error('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
            logger.error(e)
            raise errors.PopsError(e, self.pdb_file, self.chain_ids)
        else:
            logger.warning('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
            logger.warning(e)

        # Read the sasa scores from a text file
        sasa_scores = self.__read_pops_atom(chain_id)
        return sasa_scores


    def __run_pops_atom(self, chain_id):
        system_command = ('pops --noHeaderOut --noTotalOut --atomOut --pdb {0}.pdb --popsOut {0}.out'.format(chain_id))
        result, error_message, return_code = (
            helper.subprocess_check_output_locally(self.working_dir, system_command)
        )
        # The returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        logger.debug('result: %s' % result)
        output = [ line for line in result.split('\n') if line != '' ]
        return output[-1], return_code, error_message


    def __read_pops_atom(self, chain_id):
        """
        Read pops sasa results atom by atom, ignoring all main chain atoms except for Ca
        """
        # The new way
        ignore = ['N', 'C', 'O']
        per_residue_sasa_scores = []
        current_residue_number = None
        total_sasa, total_sa = None, None # so pylint stops complaining
        with open(self.working_dir + chain_id + '.out', 'r') as fh:
            for line in fh:
                row = line.split()
                if len(row) != 11:
                    continue
                atom_number, atom_name, residue_name, chain, residue_number, sasa, __, __, __, __, sa = line.split()
                atom_number, residue_number, sasa, sa = int(atom_number), int(residue_number), float(sasa), float(sa)
                if atom_name in ignore:
                    continue
                if current_residue_number != residue_number:
                    if current_residue_number:
                        per_residue_sasa_scores.append(total_sasa / float(total_sa))
                    current_residue_number = residue_number
                    total_sasa = 0
                    total_sa = 0
                total_sasa += sasa
                total_sa += sa
            per_residue_sasa_scores.append(total_sasa / float(total_sa))
        return per_residue_sasa_scores


    #%% Secondary Structure
    def get_secondary_structure(self):
        return self.get_stride()


    def get_stride(self):
        structure_file = self.get_structure_file(''.join(self.chain_ids))
        stride_results_file = op.join(
            op.dirname(structure_file),
            structure_tools.get_pdb_id(structure_file) + '_stride_results.txt'
        )
        system_command = 'stride {} -f{}'.format(structure_file, stride_results_file)
        logger.debug('stride system command: %s' % system_command)
        result, error_message, return_code = (
            helper.subprocess_check_output_locally(self.working_dir, system_command)
        )
        logger.debug('stride return code: %i' % return_code)
        logger.debug('stride result: %s' % result)
        logger.debug('stride error: %s' % error_message)
        # collect results
        with open(stride_results_file) as fh:
            file_data_df = pd.DataFrame(
                [[structure_tools.AAA_DICT[row.split()[1]], row.split()[2],
                  row.split()[3], int(row.split()[4]), row.split()[5]]
                for row in fh.readlines() if row[:3] == 'ASG'],
                columns=['amino_acid', 'chain', 'resnum', 'idx', 'ss_code'])
        return file_data_df


    def get_dssp(self):
        """Not used because crashes on server.
        """
        n_tries = 0
        return_code = -1
        while return_code != 0 and n_tries < 5:
            if n_tries > 0:
                logger.debug('Waiting for 1 minute before trying again...')
                time.sleep(60)
            system_command = ('dssp -v -i ' + ''.join(self.chain_ids) + '.pdb' + ' -o ' + 'dssp_results.txt')
            logger.debug('dssp system command: %s' % system_command)
            result, error_message, return_code = (
                helper.subprocess_check_output_locally(self.working_dir, system_command)
            )
            logger.debug('dssp return code: %i' % return_code)
            logger.debug('dssp result: %s' % result)
            logger.debug('dssp error: %s' % error_message)
            n_tries += 1
        if return_code != 0:
            if 'boost::thread_resource_error' in error_message:
                system_command = "rsync {0}{1}.pdb /home/kimlab1/strokach/tmp/elaspic_bad_pdbs/"
                helper.subprocess_check_output_locally(self.working_dir, system_command)
                raise errors.ResourceError(error_message)
        # collect results
        dssp_ss = {}
        dssp_acc = {}
        start = False
        with open(op.join(self.working_dir, 'dssp_results.txt')) as fh:
            for l in fh:
                row = l.split()
                if not row or len(row) < 2:
                    continue
                if row[1] == "RESIDUE":
                    start = True # Start parsing from here
                    continue
                if not start:
                    continue
                if l[9] == ' ':
                    continue # Skip -- missing residue
                # resseq, icode, chainid, aa, ss = int(l[5:10]), l[10], l[11], l[13], l[16]
                chainid, ss = l[11], l[16]
                if ss == ' ':
                    ss = '-'
                try:
                    acc = int(l[34:38])
                    #phi = float(l[103:109])
                    #psi = float(l[109:115])
                except ValueError as e:
                    # DSSP output breaks its own format when there are >9999
                    # residues, since only 4 digits are allocated to the seq num
                    # field.  See 3kic chain T res 321, 1vsy chain T res 6077.
                    # Here, look for whitespace to figure out the number of extra
                    # digits, and shift parsing the rest of the line by that amount.
                    if l[34] != ' ':
                        shift = l[34:].find(' ')
                        acc = int((l[34+shift:38+shift]))
                        #phi = float(l[103+shift:109+shift])
                        #psi = float(l[109+shift:115+shift])
                    else:
                        raise e
                dssp_ss.setdefault(chainid, []).append(ss) # percent sasa on sidechain
                dssp_acc.setdefault(chainid, []).append(acc)
        for key in list(dssp_ss.keys()):
            dssp_ss[key] = ''.join(dssp_ss[key])
        return dssp_ss, dssp_acc


    #%%
    def get_interchain_distances(self, pdb_chain=None, pdb_mutation=None, cutoff=None):
        """
        """
        model = self.sp.structure[0]
        shortest_interchain_distances = {}
        # Chain 1
        for i, chain_1_id in enumerate(self.sp.chain_ids):
            shortest_interchain_distances[chain_1_id] = {}
            chain_1 = model[chain_1_id]
            if pdb_chain:
                if chain_1_id == pdb_chain:
                    continue
                chain_2_ids = [pdb_chain]
            else:
                chain_2_ids = self.sp.chain_ids[i+1:]
            # Chain 2
            for chain_2_id in chain_2_ids:
                chain_2 = model[chain_2_id]
                min_r = cutoff
                # Residue 1
                for residue_1 in chain_1:
                    if (residue_1.resname not in structure_tools.AMINO_ACIDS or
                        residue_1.id[0] != ' '):
                            continue
                    # Residue 2
                    for residue_2_idx, residue_2 in enumerate(chain_2):
                        if (residue_2.resname not in structure_tools.AMINO_ACIDS or
                            residue_2.id[0] != ' '):
                                continue
                        if pdb_mutation:
                            if str(residue_2.id[1]) != pdb_mutation[1:-1]:
                                continue
                            if (structure_tools.convert_aa(residue_2.resname) != pdb_mutation[0] and
                                structure_tools.convert_aa(residue_2.resname) != pdb_mutation[-1]):
                                    logger.debug(pdb_mutation)
                                    logger.debug(structure_tools.convert_aa(residue_2.resname))
                                    logger.debug(residue_2.id)
                                    raise errors.MutationMismatchError()
                        # Atom 1
                        for atom_1 in residue_1:
                            # Atom 2
                            for atom_2 in residue_2:
                                r = structure_tools.calculate_distance(atom_1, atom_2, min_r)
                                if min_r is None or (r is not None and r < min_r):
                                    min_r = r

                shortest_interchain_distances[chain_1_id][chain_2_id] = min_r


        if not shortest_interchain_distances:
            logger.error(
                'get_interchain_distances({pdb_chain}, {pdb_mutation}, {cutoff}) failed!'
                .format(pdb_chain=pdb_chain, pdb_mutation=pdb_mutation, cutoff=cutoff)
            )
            raise Exception()

        _shortest_interchain_distances_complement = {}
        for key in shortest_interchain_distances:
            for key_2, value in shortest_interchain_distances[key].items():
                _shortest_interchain_distances_complement.setdefault(key_2, dict())[key] = value
        shortest_interchain_distances.update(_shortest_interchain_distances_complement)

        all_chains = {key for key in shortest_interchain_distances}
        all_chains.update(
            {key_2 for key in shortest_interchain_distances
             for key_2 in shortest_interchain_distances[key]}
        )

        if set(all_chains) != set(self.sp.chain_ids):
            logger.error(
                'get_interchain_distances({pdb_chain}, {pdb_mutation}, {cutoff}) failed!'
                .format(pdb_chain=pdb_chain, pdb_mutation=pdb_mutation, cutoff=cutoff)
            )
            logger.error('Did not calculate chain distances for all chain pairs!')
            logger.error('all_chains: {}'.format(all_chains))
            logger.error('self.sp.chain_ids: {}'.format(self.sp.chain_ids))
            raise Exception()

        for key_1, value_1 in shortest_interchain_distances.items():
            logger.debug(
                'Calculated interchain distances between chain {} and chains {}.'
                .format(key_1, ', '.join(list(value_1.keys()))))

        return shortest_interchain_distances


    #%%
    def get_interface_area(self, chain_ids):
        """
        """
        assert len(chain_ids) == 2

        termination, rc, e = self.__run_pops_area(self.get_structure_file(''.join(chain_ids)))
        if rc != 0:
            if termination != 'Clean termination':
                logger.error('Pops error for pdb: %s:' % self.pdb_file)
                logger.error(e)
                return [None, None, None]
        result = self.__read_pops_area(self.get_structure_file(''.join(chain_ids)) + '.out')

        # Distinguish the surface area by hydrophobic, hydrophilic, and total
        for item in result:
            if item[0] == 'hydrophobic:':
                hydrophobic = float(item[1])
            elif item[0] == 'hydrophilic:':
                hydrophilic = float(item[1])
            elif item[0] == 'total:':
                total = float(item[1])
        sasa_complex = hydrophobic, hydrophilic, total

        # calculate SASA for chain, i.e. part one of the complex:
        termination, rc, e = self.__run_pops_area(self.get_structure_file(chain_ids[0]))
        if rc != 0:
            if termination != 'Clean termination':
                logger.error('Error in pops for pdb: %s:' % self.pdb_file)
                logger.error(e)
                return [None, None, None]
        result = self.__read_pops_area(self.get_structure_file(chain_ids[0]) + '.out')

        # Distinguish the surface area by hydrophobic, hydrophilic, and total
        for item in result:
            if item[0] == 'hydrophobic:':
                hydrophobic = float(item[1])
            elif item[0] == 'hydrophilic:':
                hydrophilic = float(item[1])
            elif item[0] == 'total:':
                total = float(item[1])
        sasa_chain = hydrophobic, hydrophilic, total

        # calculate SASA for oppositeChain, i.e. the second part of the complex:
        termination, rc, e = self.__run_pops_area(self.get_structure_file(chain_ids[1]))
        if rc != 0:
            if termination != 'Clean termination':
                logger.error('Error in pops for pdb: %s:' % self.pdb_file)
                logger.error(e)
                return [None, None, None]
        result = self.__read_pops_area(self.get_structure_file(chain_ids[1]) + '.out')

        for item in result:
            if item[0] == 'hydrophobic:':
                hydrophobic = float(item[1])
            elif item[0] == 'hydrophilic:':
                hydrophilic = float(item[1])
            elif item[0] == 'total:':
                total = float(item[1])
        sasa_oppositeChain = hydrophobic, hydrophilic, total

        sasa = [ 0, 0, 0 ]
        # hydrophobic
        sasa[0] = (sasa_chain[0] + sasa_oppositeChain[0] - sasa_complex[0]) / 2.0
        # hydrophilic
        sasa[1] = (sasa_chain[1] + sasa_oppositeChain[1] - sasa_complex[1]) / 2.0
        # total
        sasa[2] = (sasa_chain[2] + sasa_oppositeChain[2] - sasa_complex[2]) / 2.0

        return sasa


    def __run_pops_area(self, full_filename):
        system_command = ('pops --chainOut'
            ' --pdb ' + full_filename +
            ' --popsOut ' + full_filename + '.out')
        result, error_message, return_code = (
            helper.subprocess_check_output_locally(self.working_dir, system_command)
        )
        # The returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        output = [ line for line in result.split('\n') if line != '' ]
        logger.debug(system_command)
#        logger.debug('pops result: %s' % result) # Prints the entire POPs output
#        logger.debug('pops error: %s' % e)
        error_message_1 = 'Warning: Atom distance too short! Probably incorrect POPS results!'
        if error_message_1 in error_message:
            logger.error(error_message_1)
        logger.debug('pops rc: %s' % return_code)
        return output[-1], return_code, error_message


    def __read_pops_area(self, filename):
        """
        This function parses POPS output that looks like this::

            === MOLECULE SASAs ===

            hydrophobic:    5267.01
            hydrophilic:    4313.68
            total:          9580.69

        """
        keep = ['hydrophobic:', 'hydrophilic:', 'total:']
        with open(filename, 'r') as pops:
            result = [ x.split(' ') for x in pops.readlines() if x != '' and x.split(' ')[0] in keep ]
            result = [ [ x.strip() for x in item if x != '' ] for item in result ]
        if not result or len(result) != 3:
            result = self.__read_pops_area_new(filename)
        return result


    def __read_pops_area_new(self, filename):
        """
        This function parses POPS output that looks like this::

            === MOLECULE SASAs ===

            Phob/A^2		Phil/A^2		Total/A^2
               5267.01	   4313.68	   9580.69

        """
        use_next_line = False
        with open(filename, 'r') as ifh:
            for line in ifh:
                if line.strip() == 'Phob/A^2\t\tPhil/A^2\t\tTotal/A^2':
                    use_next_line = True
                    continue
                if use_next_line:
                    row = line.strip().split()
                    result = [['hydrophobic:', row[0]], ['hydrophilic:', row[1]], ['total:', row[2]]]
                    break
        return result

