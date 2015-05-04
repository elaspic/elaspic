# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import zip
from builtins import object

import os
import time
from collections import OrderedDict

import six
import pandas as pd

from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser

from . import errors
from . import pdb_template
from . import helper_functions as hf


#%% CONSTANTS

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
standard_sasa_all = [[l.strip() for l in line.split()] for line in STANDARD_DATA.strip().split('\n')[1:]]
STANDARD_SASA = {x[3]: float(x[4]) for x in standard_sasa_all}



#%% STANDALONE FUNCTIONS

def get_interactions_between_chains(model, pdb_chain_1, pdb_chain_2, r_cutoff=5):
    """
    Calculate interactions between residues in pdb_chain_1 and pdb_chain_2. An
    interaction is defines as a pair of residues where at least one pair of atom
    is closer than r_cutoff. The default value for r_cutoff is 5 Angstroms.
    """
    try:
        from Bio.PDB import NeighborSearch
    except ImportError as e:
        print('Importing Biopython NeighborSearch returned an error: {}'.format(e))
        print('Using the the slow version of the neighbour-finding algorithm...')
        return get_interactions_between_chains_slow(model, pdb_chain_1, pdb_chain_2, r_cutoff)
        
    # Extract the chains of interest from the model
    chain_1 = None
    chain_2 = None
    for child in model.get_list():
        if child.id == pdb_chain_1:
            chain_1 = child
        if child.id == pdb_chain_2:
            chain_2 = child
    if chain_1 is None or chain_2 is None:
        raise Exception('Chains %s and %s were not found in the model' % (pdb_chain_1, pdb_chain_2))

    ns = NeighborSearch(list(chain_2.get_atoms()))
    interactions_between_chains = OrderedDict()
    for idx, residue_1 in enumerate(chain_1):
        if residue_1.resname in pdb_template.amino_acids and residue_1.id[0] == ' ':
            resnum_1 = str(residue_1.id[1]) + residue_1.id[2].strip()
            resaa_1 = pdb_template.convert_aa(residue_1.get_resname())
            interacting_residues = set()
            for atom_1 in residue_1:
                interacting_residues.update(ns.search(atom_1.get_coord(), r_cutoff, 'R'))
            interacting_resids = []
            for residue_2 in interacting_residues:
                resnum_2 = str(residue_2.id[1]) + residue_2.id[2].strip()
                resaa_2 = pdb_template.convert_aa(residue_2.get_resname())
                if residue_2.resname in pdb_template.amino_acids and residue_2.id[0] == ' ':
                    interacting_resids.append((resnum_2, resaa_2,))
            if interacting_resids:
                interacting_resids.sort(key=lambda x: int(''.join([c for c in x[0] if c.isdigit()])))
                interactions_between_chains[(resnum_1, resaa_1)] = interacting_resids
    return interactions_between_chains


def get_interactions_between_chains_slow(model, pdb_chain_1, pdb_chain_2, r_cutoff=5):
    """
    Calculate interactions between residues in pdb_chain_1 and pdb_chain_2. An
    interaction is defines as a pair of residues where at least one pair of atom
    is closer than r_cutoff. The default value for r_cutoff is 5 Angstroms.
    """
    # Extract the chains of interest from the model
    chain_1 = None
    chain_2 = None
    for child in model.get_list():
        if child.id == pdb_chain_1:
            chain_1 = child
        if child.id == pdb_chain_2:
            chain_2 = child
    if chain_1 is None or chain_2 is None:
        raise Exception('Chains %s and %s were not found in the model' % (pdb_chain_1, pdb_chain_2))

    interactions_between_chains = OrderedDict()
    for idx, residue_1 in enumerate(chain_1):
        if residue_1.resname in pdb_template.amino_acids and residue_1.id[0] == ' ':
            resnum_1 = str(residue_1.id[1]) + residue_1.id[2].strip()
            resaa_1 = pdb_template.convert_aa(residue_1.get_resname())
            interacting_resids = []
            for residue_2 in chain_2:
                resnum_2 = str(residue_2.id[1]) + residue_2.id[2].strip()
                resaa_2 = pdb_template.convert_aa(residue_2.get_resname())
                r_min = None
                if residue_2.resname in pdb_template.amino_acids and residue_2.id[0] == ' ':
                    for atom_1 in residue_1:
                        for atom_2 in residue_2:
                            r = pdb_template.calculate_distance(atom_1, atom_2, r_cutoff)
                            if r is not None:
                                if r_min and r < r_min:
                                    r_min = r
                                elif not r_min:
                                    r_min = r
                if r_min:
                    interacting_resids.append((resnum_2, resaa_2, r_min,))
            if interacting_resids:
                interactions_between_chains[(resnum_1, resaa_1)] = interacting_resids
    return interactions_between_chains



#%%

class PhysiChem(object):

    def __init__(self, vdW, d, unique, logger):
        self.vdW_distance = float(vdW)
        self.contact_distance = float(d)
        self.logger = logger


    def __call__(self, pdb_filename, mutated_chain_id, mutation):
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
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        structure = parser.get_structure('ID', pdb_filename)
        model = structure[0]
        mutated_chain = model[mutated_chain_id]
        opposite_chains = [ chain for chain in model.child_list if chain.id != mutated_chain_id ]
        main_chain_atoms = ['CA', 'C', 'N', 'O']

        # Find the mutated residue
        mutation_position = int(mutation[1:-1])
        residue_counter = 0
        for residue in mutated_chain:
            if residue.resname in pdb_template.amino_acids and residue.id[0] == ' ':
                residue_counter += 1
                if residue_counter == mutation_position:
                    if (pdb_template.AAA_DICT[residue.resname] not in
                            [mutation[0].upper(), mutation[-1].upper()]):
                        self.logger.error(residue)
                        self.logger.error(pdb_template.AAA_DICT[residue.resname])
                        self.logger.error(mutation.upper())
                        self.logger.error(mutation_position)
                        raise Exception('PhysiChem detected mutated amino acid position mismatch!')
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
                r = pdb_template.calculate_distance(mutated_atom, partner_atom, self.vdW_distance)
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
                    if r <= self.contact_distance:
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



#%%

class AnalyzeStructure(object):
    """
    Runs the program pops to calculate the interface size of the complexes
    This is done by calculating the surface of the complex and the seperated parts.
    The interface is then given by the substracting
    """

    def __init__(self, data_path, working_path, pdb_file, chains, domain_defs, logger):
        self.data_path = data_path #: folder with the structures (modeller_path, foldx_path, etc,)
        self.working_path = working_path # analyze_structure path with all the binaries
        self.pdb_file = pdb_file
        self.chain_ids = chains
        self.domain_defs = domain_defs
        self.logger = logger
        self.structure = self.__split_pdb_into_chains()


    def __split_pdb_into_chains(self):
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        io = PDBIO()
        self.logger.debug('Saving parsed pdbs into the following working path:')
        self.logger.debug(self.data_path + self.pdb_file)

        # Save all chains together with correct chain letters
        self.logger.debug('Saving all chains together')
        structure = parser.get_structure('ID', self.data_path + self.pdb_file)
        if len(structure) > 1:
            # Delete all models except for the first (otherwise get errors with naccess)
            del structure[1:]
        model = structure[0]
        children = model.get_list()
#        for child in children:
#            self.logger.debug('child id before:' + child.id)
#            if child.id not in self.chain_ids:
#            self.logger.debug('child id after:' + child.id)
        if len(children) == 1 and (children[0].id == '' or children[0].id == ' ' or children[0].id == '0'):
            children[0].id = self.chain_ids[0]
        if len(children) == 2 and children[0].id == '0' and children[1].id == '0':
            children[0].id = 'A'
            children[1].id = 'B'
        io.set_structure(structure)
        outFile = self.working_path +  self.pdb_file
        io.save(outFile)

        if len(self.chain_ids) > 1:
            # Save a structure with only the chains of interest
            self.logger.debug('Saving only the chains of interest: %s' % ','.join(self.chain_ids))
            structure = parser.get_structure('ID', self.data_path + self.pdb_file)
            model = structure[0]
            for child in model.get_list():
                self.logger.debug('child id:' + child.id)
                if child.id not in self.chain_ids:
                    self.logger.debug('detaching chain')
                    model.detach_child(child.id)
            io.set_structure(structure)
            outFile = self.working_path + ''.join(self.chain_ids) + '.pdb'
            io.save(outFile)

        # Save a structure for each chain
        for chain_id in self.chain_ids:
            self.logger.debug('Saving chain %s separately' % chain_id)
            # save chain, i.e. part one of the complex:
            structure = parser.get_structure('ID', self.data_path + self.pdb_file)
            model = structure[0]
            for child in model.get_list():
                self.logger.debug('child id:' + child.id)
                if child.id != chain_id:
                    self.logger.debug('detaching chain %s' % child.id)
                    model.detach_child(child.id)
            io.set_structure(structure)
            outFile = self.working_path + chain_id + '.pdb'
            io.save(outFile)

        # The main class structure is the one that only has the chains of interest
        structure = parser.get_structure('ID', self.working_path + ''.join(self.chain_ids) + '.pdb')
        return structure


    #%%
    def get_seasa(self):
        seasa_by_chain_together, seasa_by_residue_together = self._run_msms(''.join(self.chain_ids) + '.pdb')
        if len(self.chain_ids) > 1:
            seasa_by_chain_separately = []
            seasa_by_residue_separately = []
            for chain_id in self.chain_ids:
                seasa_by_chain, seasa_by_residue = self._run_msms(chain_id + '.pdb')
                seasa_by_chain_separately.append(seasa_by_chain)
                seasa_by_residue_separately.append(seasa_by_residue)
            seasa_by_chain_separately = pd.concat(seasa_by_chain_separately, ignore_index=True)
            seasa_by_residue_separately = pd.concat(seasa_by_residue_separately, ignore_index=True)
            return [seasa_by_chain_together, seasa_by_chain_separately, seasa_by_residue_together, seasa_by_residue_separately]
        else:
            return [seasa_by_chain_together, seasa_by_chain_together, seasa_by_residue_together, seasa_by_residue_together]


    def _run_msms(self, filename):
        """ In the future, could add an option to measure residue depth
        using Bio.PDB.ResidueDepth().residue_depth()...
        """
        base_filename = filename[:filename.rfind('.')]

        # Convert pdb to xyz coordiates
        assert(os.path.isfile(self.working_path + filename))
        system_command = 'pdb_to_xyzrn {0}.pdb'.format(self.working_path + base_filename)
        self.logger.debug('msms system command 1: %s' % system_command)
        child_process = hf.run_subprocess_locally(self.working_path, system_command)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        return_code = child_process.returncode
        if return_code != 0:
            self.logger.debug('msms result 1:')
            self.logger.debug(result)
            self.logger.debug('msms error message 1:')
            self.logger.debug(error_message)
            self.logger.debug('naccess rc 1:')
            self.logger.debug(child_process.returncode)
            raise errors.MSMSError(error_message)
        else:
            with open(self.working_path + base_filename + '.xyzrn', 'w') as ofh:
                ofh.writelines(result)

        # Calculate SASA and SESA (excluded)
        probe_radius = 1.4
        system_command_string = (
            'msms '
            '-probe_radius {1:.1f} '
            '-surface ases '
            '-if {0}.xyzrn '
            '-af {0}.area')
        system_command = system_command_string.format(self.working_path + base_filename, probe_radius)
        self.logger.debug('msms system command 2: %s' % system_command)
        child_process = hf.run_subprocess_locally(self.working_path, system_command)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        return_code = child_process.returncode
        number_of_tries = 0
        while return_code != 0 and number_of_tries < 5:
            self.logger.error('MSMS exited with an error!')
            probe_radius -= 0.1
            self.logger.debug('Reducing probe radius to {}'.format(probe_radius))
            system_command = system_command_string.format(self.working_path + base_filename, probe_radius)
            child_process = hf.run_subprocess_locally(self.working_path, system_command)
            result, error_message = child_process.communicate()
            if six.PY3:
                result = str(result, encoding='utf-8')
                error_message = str(error_message, encoding='utf-8')
            return_code = child_process.returncode
            number_of_tries += 1
        if return_code != 0:
            self.logger.debug('msms result 2:')
            self.logger.debug(result)
            self.logger.debug('msms error message 2:')
            self.logger.debug(error_message)
            self.logger.debug('naccess rc 2:')
            self.logger.debug(child_process.returncode)
            raise errors.MSMSError(error_message)

        # Read and parse the output
        with open(self.working_path + base_filename + '.area', 'r') as fh:
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
            
#        seasa_df['chain'] = seasa_df['atom_num'].apply(lambda x: atom_to_chain.get(x, None))
#        seasa_df['res_name'] = seasa_df['atom_num'].apply(lambda x: atom_to_res_name.get(x, None))
#        seasa_df['res_num'] = seasa_df['atom_num'].apply(lambda x: atom_to_res_num.get(x, None))
#        seasa_df['atom_id'] = seasa_df['atom_num'].apply(lambda x: atom_to_atom_id.get(x, None))
#        seasa_df.dropna(inplace=True)

#        incorrect_assignments = \
#            [x for x in zip(seasa_df['res_name'], seasa_df['res_name_msms']) if x[0]!=x[1]] + \
#            [x for x in zip(seasa_df['atom_id'], seasa_df['atom_id_msms']) if x[0]!=x[1]]
#        if len(incorrect_assignments) > 5:
#            self.logger.error('Could not correctly assign msms output to chains!')
#            self.logger.error(incorrect_assignments)
#            raise errors.MSMSError('Could not correctly assign msms output to chains: ' + str(incorrect_assignments))
#        del seasa_df['atom_id_msms']
#        del seasa_df['res_name_msms']
#        del seasa_df['res_num_msms']

        seasa_gp_by_chain = seasa_df.groupby(['pdb_chain'])
        seasa_gp_by_residue = seasa_df.groupby(['pdb_chain', 'res_name', 'res_num'])
        seasa_by_chain = seasa_gp_by_chain.sum().reset_index()
        seasa_by_residue = seasa_gp_by_residue.sum().reset_index()

        return seasa_by_chain, seasa_by_residue


    #%%
    def get_sasa(self, program_to_use='naccess'):
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


    # NACESS is not used anymore
    # It is not reliable enough on a PDB-wide scale
#    def _run_naccess_atom(self, filename):
#        # run naccess
#        system_command = ('naccess ' + filename)
#        self.logger.debug('naccess system command: %s' % system_command)
#        assert(os.path.isfile(self.working_path + filename))
#        child_process = hf.run_subprocess_locally(self.working_path, system_command)
#        result, error_message = child_process.communicate()
#        if six.PY3:
#            result = str(result, encoding='utf-8')
#            error_message = str(error_message, encoding='utf-8')
#        return_code = child_process.returncode
#        self.logger.debug('naccess result: {}'.format(result))
#        self.logger.debug('naccess error: {}'.format(error_message))
#        self.logger.debug('naccess rc: {}'.format(return_code))
#        # Collect results
#        sasa_scores = {}
#        with open(self.working_path + filename.split('.')[0] + '.rsa') as fh:
#            for line in fh:
#                row = line.split()
#                if row[0] != 'RES':
#                    continue
#                try:
#                    (line_id, res, chain, num, all_abs, all_rel,
#                    sidechain_abs, sidechain_rel, mainchain_abs, mainchain_rel,
#                    nonpolar_abs, nonpolar_rel, polar_abs, polar_rel) = row
#                except ValueError as e:
#                    print(e)
#                    print(line)
#                    print(row)
#                sasa_scores.setdefault(chain, []).append(sidechain_rel) # percent sasa on sidechain
#        return sasa_scores


    def _run_pops_atom(self, chain_id):
        # Use pops to calculate sasa score for the given chain
        termination, rc, e = self.__run_pops_atom(chain_id)
        if termination != 'Clean termination':
            self.logger.error('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
            self.logger.error(e)
            raise errors.PopsError(e, self.data_path + self.pdb_file, self.chain_ids)
        else:
            self.logger.warning('Pops error for pdb: %s, chains: %s: ' % (self.pdb_file, ' '.join(self.chain_ids),) )
            self.logger.warning(e)

        # Read the sasa scores from a text file
        sasa_scores = self.__read_pops_atom(chain_id)
        return sasa_scores


    def __run_pops_atom(self, chain_id):
        system_command = ('pops --noHeaderOut --noTotalOut --atomOut --pdb {0}.pdb --popsOut {0}.out'.format(chain_id))
        child_process = hf.run_subprocess_locally(self.working_path, system_command)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        return_code = child_process.returncode
        # The returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        self.logger.debug('result: %s' % result)
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
        with open(self.working_path + chain_id + '.out', 'r') as fh:
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


    #%%
    def get_secondary_structure(self):
        return self.get_stride()


    def get_stride(self):
        system_command = 'stride ' + ''.join(self.chain_ids) + '.pdb ' + '-fstride_results.txt'
        self.logger.debug('stride system command: %s' % system_command)
        child_process = hf.run_subprocess_locally(self.working_path, system_command)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        return_code = child_process.returncode
        self.logger.debug('stride return code: %i' % return_code)
        self.logger.debug('stride result: %s' % result)
        self.logger.debug('stride error: %s' % error_message)
        # collect results
        with open(self.working_path + 'stride_results.txt') as fh:
            file_data_df = pd.DataFrame(
                [[pdb_template.AAA_DICT[row.split()[1]], row.split()[2], row.split()[3], int(row.split()[4]), row.split()[5]]
                for row in fh.readlines() if row[:3] == 'ASG'],
                columns=['amino_acid', 'chain', 'resnum', 'idx', 'ss_code'])
        return file_data_df


    def get_dssp(self):
        """ Not used because crashes on server
        """
        n_tries = 0
        return_code = -1
        while return_code != 0 and n_tries < 5:
            if n_tries > 0:
                self.logger.debug('Waiting for 1 minute before trying again...')
                time.sleep(60)
            system_command = ('dssp -v -i ' + ''.join(self.chain_ids) + '.pdb' + ' -o ' + 'dssp_results.txt')
            self.logger.debug('dssp system command: %s' % system_command)
            child_process = hf.run_subprocess_locally(self.working_path, system_command)
            result, error_message = child_process.communicate()
            if six.PY3:
                result = str(result, encoding='utf-8')
                error_message = str(error_message, encoding='utf-8')
            return_code = child_process.returncode
            self.logger.debug('dssp return code: %i' % return_code)
            self.logger.debug('dssp result: %s' % result)
            self.logger.debug('dssp error: %s' % error_message)
            n_tries += 1
        if return_code != 0:
            if 'boost::thread_resource_error' in error_message:
                system_command = "rsync {0}{1}.pdb /home/kimlab1/strokach/tmp/elaspic_bad_pdbs/"
                hf.run_subprocess_locally(self.working_path, system_command)
                raise errors.ResourceError(error_message)
        # collect results
        dssp_ss = {}
        dssp_acc = {}
        start = False
        with open(self.working_path + 'dssp_results.txt') as fh:
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
                resseq, icode, chainid, aa, ss = int(l[5:10]), l[10], l[11], l[13], l[16]
                if ss == ' ':
                    ss = '-'
                try:
                    acc = int(l[34:38])
                    phi = float(l[103:109])
                    psi = float(l[109:115])
                except ValueError as e:
                    # DSSP output breaks its own format when there are >9999
                    # residues, since only 4 digits are allocated to the seq num
                    # field.  See 3kic chain T res 321, 1vsy chain T res 6077.
                    # Here, look for whitespace to figure out the number of extra
                    # digits, and shift parsing the rest of the line by that amount.
                    if l[34] != ' ':
                        shift = l[34:].find(' ')
                        acc = int((l[34+shift:38+shift]))
                        phi = float(l[103+shift:109+shift])
                        psi = float(l[109+shift:115+shift])
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
        model = self.structure[0]
        shortest_interchain_distances = {}
        for chain_1 in [chain for chain in model]:
            if pdb_chain and chain_1.id != pdb_chain:
                continue # skip chains that we are not interested in
            shortest_interchain_distances[chain_1.id] = {}
            for idx, residue_1 in enumerate(chain_1):
                if residue_1.resname in pdb_template.amino_acids and residue_1.id[0] == ' ':
                    if pdb_mutation:
                        if str(residue_1.id[1]) != pdb_mutation[1:-1]:
                            continue # skip all residues that we are not interested in
                        if (pdb_template.convert_aa(residue_1.resname) != pdb_mutation[0] and
                            pdb_template.convert_aa(residue_1.resname) != pdb_mutation[-1]):
                                self.logger.debug(pdb_mutation)
                                self.logger.debug(pdb_template.convert_aa(residue_1.resname))
                                self.logger.debug(residue_1.id)
                                raise Exception
                    for chain_2 in [chain for chain in model if chain != chain_1]:
                        min_r = cutoff
                        for residue_2 in chain_2:
                            if (residue_1.resname not in pdb_template.amino_acids or
                                residue_2.resname not in pdb_template.amino_acids):
                                    continue
                            for atom_1 in residue_1:
                                for atom_2 in residue_2:
                                    r = pdb_template.calculate_distance(atom_1, atom_2, min_r)
                                    if r and (not min_r or min_r > r):
                                        min_r = r
                        shortest_interchain_distances[chain_1.id][chain_2.id] = min_r

        if (not shortest_interchain_distances or
            (pdb_chain and not shortest_interchain_distances[pdb_chain])):
                self.logger.error(
                    'Could not calculate the shortest interchain distance for chain {} residue {}'
                    .format(pdb_chain, pdb_mutation))
                raise Exception()

        for key_1, value_1 in shortest_interchain_distances.items():
            self.logger.debug(
                'Calculated interchain distances between chain {} and chains {}'
                .format(key_1, ', '.join(list(value_1.keys()))))

        return shortest_interchain_distances


    #%%
    def get_interface_area(self):

        termination, rc, e = self.__run_pops_area(self.working_path + ''.join(self.chain_ids) + '.pdb')
        if rc != 0:
            if termination != 'Clean termination':
                self.logger.error('Pops error for pdb: %s:' % self.pdb_file)
                self.logger.error(e)
                return '0', '0', '0'
        result = self.__read_pops_area(self.working_path + ''.join(self.chain_ids) + '.out')

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
        termination, rc, e = self.__run_pops_area(self.working_path + self.chain_ids[0] + '.pdb')
        if rc != 0:
            if termination != 'Clean termination':
                self.logger.error('Error in pops for pdb: %s:' % self.pdb_file)
                self.logger.error(e)
                return '0', '0', '0'
        result = self.__read_pops_area(self.working_path + self.chain_ids[0] + '.out')

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
        termination, rc, e = self.__run_pops_area(self.working_path + self.chain_ids[1] + '.pdb')
        if rc != 0:
            if termination != 'Clean termination':
                self.logger.error('Error in pops for pdb: %s:' % self.pdb_file)
                self.logger.error(e)
                return '0', '0', '0'
        result = self.__read_pops_area(self.working_path + self.chain_ids[1] + '.out')

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
            ' --popsOut ' + self.working_path + full_filename.split('/')[-1].replace('pdb', 'out'))
        child_process = hf.run_subprocess_locally(self.working_path, system_command)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        return_code = child_process.returncode
        # The returncode can be non zero even if pops calculated the surface
        # area. In that case it is indicated by "clean termination" written
        # to the output. Hence this check:
        # if output[-1] == 'Clean termination' the run should be OK
        output = [ line for line in result.split('\n') if line != '' ]
        self.logger.debug(system_command)
#        self.logger.debug('pops result: %s' % result) # Prints the entire POPs output
#        self.logger.debug('pops error: %s' % e)
        error_message_1 = 'Warning: Atom distance too short! Probably incorrect POPS results!'
        if error_message_1 in error_message:
            self.logger.error(error_message_1)
        self.logger.debug('pops rc: %s' % return_code)
        return output[-1], return_code, error_message


    def __read_pops_area(self, filename):
        # The old way
        keep = ['hydrophobic:', 'hydrophilic:', 'total:']
        with open(filename, 'r') as pops:
            result = [ x.split(' ') for x in pops.readlines() if x != '' and x.split(' ')[0] in keep ]
        return [ [ x.strip() for x in item if x != '' ] for item in result ]


#%%
if __name__ == '__main__':
    # Insert debug code here
    pass