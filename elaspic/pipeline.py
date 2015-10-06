# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import object

import os
import os.path as op
import json
import signal
import six
import logging 
from functools import wraps

from . import conf, errors, helper, structure_tools, structure_analysis

logger = logging.getLogger(__name__)
configs = conf.Configs()

sql_db = None
domain_alignment = None
domain_model = None
domain_mutation = None

MAX_DISTANCE_BETWEEN_INTERACTING_CHAINS = 6 # Angstrom
ELASPIC_LOGO = """

8888888888 888             d8888  .d8888b.  8888888b. 8888888 .d8888b.  
888        888            d88888 d88P  Y88b 888   Y88b  888  d88P  Y88b 
888        888           d88P888 Y88b.      888    888  888  888    888 
8888888    888          d88P 888  "Y888b.   888   d88P  888  888        
888        888         d88P  888     "Y88b. 8888888P"   888  888        
888        888        d88P   888       "888 888         888  888    888 
888        888       d8888888888 Y88b  d88P 888         888  Y88b  d88P 
8888888888 88888888 d88P     888  "Y8888P"  888       8888888 "Y8888P"  

"""

#%% 
class Pipeline(object):

    def __init__(self, configurations):
        """
        It should be possible to initialize one pipeline and call it in parallel using different
        mutations as input
        """       
        # Read the configuration file and set the variables
        if isinstance(configurations, six.string_types):
            conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            configs.update(**configurations)
        
        # TODO: remove so error message does not appear in a production release.  
        self._validate_temp_path(conf.configs)

        # Initialize a logger
        for line in ELASPIC_LOGO.split('\n'):
            logger.info(line)

        self.PWD = os.getcwd()
        

    def _validate_temp_path(self, configs):
        """
        Make sure that we are using a job specific temporary folder if we are on a cluster.
        """
        hostname = helper.get_hostname()
        no_job_specific_folder = configs['temp_dir'].startswith(configs['global_temp_dir'])
        on_node_with_manditory_job_specific_folder = (
            any([(x.lower() in hostname) for x in ['node', 'behemoth', 'grendel', 'beagle']])
        )
        if no_job_specific_folder and on_node_with_manditory_job_specific_folder:
            raise Exception('You should be using a temp folder that it specific to the particular job!')
    

    def run(self):
        raise NotImplementedError()
        
        
    def run_sequence(self, sequence):
        raise NotImplementedError()


    def _get_model(self):
        raise NotImplementedError()
    
    
    def run_model(self, model):
        raise NotImplementedError()


    def _get_mutation(self, model, mutation_data):
        raise NotImplementedError


    def run_mutation(self, mutation):
        raise NotImplementedError()




#%%
class StructurePipeline(Pipeline):

    def __init__(self, pdb_file, configurations=None):
        """
        TODO: Add an option to store provean results based on sequence hash.
        """
        super().__init__(configurations)
        
        # Input parameters
        self.pdb_id = op.basename(pdb_file).rstrip('.pdb')
        self.pdb_file = pdb_file
        
        logger.info('pdb_file: {}'.format(self.pdb_file))
        logger.info('pwd: {}'.format(self.PWD))
        
        ### Load PDB structure and extract required sequences and chains.
        #fix_pdb(self.pdb_file, self.pdb_file)
        self.sp = structure_tools.StructureParser(self.pdb_file)
        self.sp.extract()
        self.sp.save_structure()
        self.sp.save_sequences()
        
        self.pdb_chains = [chain.id for chain in self.sp.structure[0].child_list]
    

    def run(self):
        """Run all steps of the pipeline.
        
        Made this a separate method because it makes it easier to test
        :py:func:`compute_provean`, :py:func:`compute_model`, and :py:func:`compute_mutation`.
        """
        ...
        

    def calculate_provean(self, pdb_chain):
        """
        Raises
        -------
        errors.ProveanError
        errors.ProveanResourceError
        """
        protein_id, protein_name, provean_supset_file, chain_sequence = (
            self.get_provean_parameters(pdb_chain)
        )
        try:
            provean_supset_filename, provean_supset_length = (
                domain_alignment.build_provean_supporting_set(
                    protein_id, protein_name, chain_sequence, conf.configs, 
                    self.unique_temp_folder, self.provean_temp_path,
                    provean_supset_file)
            )
        except errors.ProveanError as e:
            logger.error(e)
            raise
        except errors.ProveanResourceError as e:
            logger.error(e)
            # Send the kill signal to the main process group, killing everything
            self._clear_provean_temp_files() # these won't get cleaned once the process dies
            logger.error('Killing group...')
            os.killpg(e.child_process_group_id, signal.SIGTERM)
            error_message = (
                'Provean ran out of resources. Everything has been killed. '
                'The user will probably not see this message...'
            )
            logger.critical(error_message)
            raise
        finally:
            self._clear_provean_temp_files()

        logger.info('Finished computing provean.')
        provean_dict = {
            'provean_supset_filename': provean_supset_filename,
            'provean_supset_length': provean_supset_length,
        }
        return provean_dict


    def get_provean_parameters(self, pdb_chain):
        """Information needed to run Provean.
        """
        protein_id = "{}{}".format(self.pdb_id, pdb_chain)
        protein_name = "{}.chain{}".format(op.basename(self.pdb_file), pdb_chain)
        provean_supset_file = op.join(
            self.pwd, 'sequence_conservation', "{}_provean_supset".format(protein_id))
        chain_sequence, chain_numbering_extended = self.sp.get_chain_sequence_and_numbering(pdb_chain)
        return protein_id, protein_name, provean_supset_file, chain_sequence
       
   
    def get_structure_sequence(self):
        """
        """
        structure_sequence = []
        # After running `sp.extract()`, all but the last chain are guaranteed to be without HATATMS
        for chain in self.sp.structure.child_list[0].child_list[:-1]:
            chain_sequence, chain_numbering_extended = self.sp.get_chain_sequence_and_numbering(chain.id)
            structure_sequence += [chain_sequence]
        
        # If the last chain starts with a HETATM, it should be entirely HETATMs.
        last_chain = self.sp.structure.child_list[0].child_list[-1]
        if last_chain.child_list[0].resname in structure_tools.AAA_DICT:
            chain_sequence, chain_numbering_extended = self.sp.get_chain_sequence_and_numbering(last_chain.id)
            structure_sequence += [chain_sequence]
        else:
            structure_sequence += ['.' * len(last_chain.child_list)]
            
        structure_sequence = '/'.join(structure_sequence)
        return structure_sequence
    
    
    def calculate_model(self):
        """
        Returns
        -------
        model_dict : dict
            Contains all important values calculated by modeller.
            
        Raises
        -------
        errors.ModellerError
        errors.PDBChainError
        errors.PDBEmptySequenceError
        errors.PDBNotFoundError
        """
        structure_sequence = self.get_structure_sequence()

        modeller_path = op.join(self.unique_temp_folder, 'modeller')
        os.makedirs(modeller_path, mode=0o777, exist_ok=True)

        modeller_target_id = self.pdb_id + '_modeller'
        modeller_template_id = self.pdb_id

        pir_alignment_filename = op.join(modeller_path, modeller_target_id + '.pir')
        with open(pir_alignment_filename, 'w') as ofh:
            domain_model.write_to_pir_alignment(ofh, 'sequence', modeller_target_id, structure_sequence)
            domain_model.write_to_pir_alignment(ofh, 'structure', modeller_template_id, structure_sequence)           
        
        new_chains = ''.join([chain.id for chain in self.sp.structure[0].child_list])
        model_norm_dope, model_file, model_knotted = domain_model.run_modeller(
            pir_alignment_filename, [modeller_target_id], [modeller_template_id], 
            self.unique_temp_folder, conf.configs['modeller_runs'],
            new_chains
        )
        
        model_dict = {
            'model_norm_dope': model_norm_dope,
            'model_file': model_file,
            'model_knotted': model_knotted,
        }
        return model_dict
        

    def calculate_mutation(self, pdb_chain, pdb_mutation, provean_results, modeller_results):
        """
        Raises
        ------
        errors.PDBError
        errors.FoldxError
        errors.ResourceError
        errors.FoldXAAMismatchError
        errors.MutationOutsideDomainError
        errors.MutationOutsideInterfaceError
        """
        mut_data = self.get_mutation_data(pdb_chain, pdb_mutation, provean_results, modeller_results)
        mutation.get_provean_score(mut_data)
        mutation_interactions = self.get_interactions(pdb_chain, pdb_mutation)
        
        return mutation_interactions

    
    def get_mutation_data(self, pdb_chain, pdb_mutation, provean_results, modeller_results):
        """
        Create a MutData class that holds all the information we need
        about a given mutation.
        """    
        protein_id = "{}{}".format(self.pdb_id, pdb_chain)
        logger.info(provean_results)
        provean_supset_filename = provean_results['provean_supset_filename']
        model_file = modeller_results['model_file']
        pdb_chains = [chain.id for chain in self.sp.structure[0].child_list]

        mutation_chain_pos = 1
        for aa, aa_pos in zip(*self.sp.get_chain_sequence_and_numbering(pdb_chain)):
            if aa_pos != pdb_mutation[1:-1]:
                mutation_chain_pos += 1
            else:
                break
        mutation_chain = pdb_mutation[0] + str(mutation_chain_pos) + pdb_mutation[-1]
        
        mutation_modeller = self.mutation_to_modeller(pdb_chain, pdb_mutation)
        mutation_modeller_pos = mutation_modeller[1:-1]

        ### Assign values to appropriate variables
        mut_data = domain_mutation.MutationData()
        
        # Mutation info
        mut_data.uniprot_id_1 = protein_id 
        mut_data.mutation = mutation_modeller
        mut_data.chains_modeller = pdb_chains
        mut_data.domain_sequences = [
            self.sp.get_chain_sequence_and_numbering(chain_id)[0]
            for chain_id in pdb_chains
        ]
        mut_data.path_to_provean_supset = op.join(
            self.unique_temp_folder, 'sequence_conservation', provean_supset_filename
        )    
        mut_data.save_path = self.unique_temp_folder
        mut_data.pdbFile_wt = op.join(self.unique_temp_folder, 'modeller', model_file)
        mut_data.position_domain = mutation_chain_pos
        mut_data.mutation_domain = mutation_chain
        mut_data.position_modeller = mutation_modeller_pos
        mut_data.mutation_modeller = mutation_modeller
        
        # Provean results
        mut_data.provean_mutation = None
        mut_data.provean_score = None

        # Validate results            
        mut_data.validate()
               
        # Print all attributes
#        for attr in dir(mut_data):
#            if not attr.startswith('_'):
#                logger.debug("{}: {}".format(attr, getattr(mut_data, attr)))

        return mut_data            
    
    
    def _mutation_to_modeller(self, pdb_chain, pdb_mutation):
        """
        Modeller numbers all residues from 1 to n_residues, for all chains. 
        Use this function to convert mutations from PDB coordinates to Modeller coordinates.
        """
        mutation_id = pdb_mutation[:-1]
        mutation_idx = 1
        
        for chain in self.sp.structure.child_list[0].child_list:
            if chain.id != pdb_chain:
                mutation_idx += len(chain.child_list)
            else:
                break
            
        for residue in chain.child_list:
            residue_id = (structure_tools.AAA_DICT[residue.resname] + str(residue.id[1]))
            if residue_id != mutation_id:
                mutation_idx += 1
            else:
                break
    
        if residue_id != mutation_id:
            return None
        else:
            mutation_modeller = "{}{}{}".format(pdb_mutation[0], mutation_idx, pdb_mutation[-1])
            return mutation_modeller
    
    
    
    def _get_interactions(self, pdb_chain, pdb_mutation):
        interactions = structure_analysis.get_interactions(
            self.sp.structure.child_list[0], pdb_chain, self.R_CUTOFF)

        mutation_contacts = []
        for chain in interactions:
            for aa_pos, aa in interactions[chain]:
                if aa + aa_pos == pdb_mutation[:-1]:
                    mutation_contacts.append(chain)

        return mutation_contacts



        
        
#%% PDB only
#class UserStructureMutation(Base):
#    __tablename__ = 'user_structure_mutation'
#    _indexes = [
#        ['pdb_id', 'pdb_chain', 'pdb_chain_2', 'pdb_mutation'],
#    ]
#    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
#    
#    id = sa.Column(sa.Integer, nullable=False, primary_key=True, autoincrement=True)
#    
#    ### Input
#    unique = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')), 
#        primary_key=True)
#    pdb_file = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    pdb_id = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    pdb_chains = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    
#    
#    ### Mutation    
#    # Original PDB coordinates
#    pdb_mutation = sa.Column(sa.String(SHORT), nullable=False)
#    # Where 1 is the start of the domain (i.e. chain)
#    domain_mutation = sa.Column(sa.String(SHORT), nullable=False)
#    # Where 1 is the start of the protein (i.e. chain A residue 1)
#    modeller_mutation = sa.Column(sa.String(SHORT), nullable=False)
#    
#    
#    ### Homology Model
#    # Core / Interface
#    model_filename = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    norm_dope = sa.Column(sa.Float)
#    sasa_score = sa.Column(sa.Text)
#        
#    # Interface
#    pdb_chain = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    pdb_chain_2 = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    interface_area_hydrophobic = sa.Column(sa.Float)
#    interface_area_hydrophilic = sa.Column(sa.Float)
#    interface_area_total = sa.Column(sa.Float)
#    interacting_aa_1 = sa.Column(sa.Text)
#    interacting_aa_2 = sa.Column(sa.Text)
#    
#    
#    
#    ### Mutation
#    mutation = sa.Column(sa.String(SHORT), nullable=False)
#    mutation_errors = sa.Column(sa.Text)
#    model_filename_wt = sa.Column(sa.String(MEDIUM))
#    model_filename_mut = sa.Column(sa.String(MEDIUM))
#    chain_modeller = sa.Column(sa.String(SHORT))
#    mutation_modeller = sa.Column(sa.String(SHORT))
#    analyse_complex_energy_wt = sa.Column(sa.Text)
#    stability_energy_wt = sa.Column(sa.Text)
#    analyse_complex_energy_mut = sa.Column(sa.Text)
#    stability_energy_mut = sa.Column(sa.Text)
#    physchem_wt = sa.Column(sa.Text)
#    physchem_wt_ownchain = sa.Column(sa.Text)
#    physchem_mut = sa.Column(sa.Text)
#    physchem_mut_ownchain = sa.Column(sa.Text)
#    matrix_score = sa.Column(sa.Float)
#    secondary_structure_wt = sa.Column(sa.Text)
#    solvent_accessibility_wt = sa.Column(sa.Float)
#    secondary_structure_mut = sa.Column(sa.Text)
#    solvent_accessibility_mut = sa.Column(sa.Float)
#    contact_distance_wt = sa.Column(sa.Float)
#    contact_distance_mut = sa.Column(sa.Float)
#    provean_score = sa.Column(sa.Float)
#    ddg = sa.Column(sa.Float, index=False)
#    mut_date_modified = sa.Column(
#        sa.DateTime, default=datetime.datetime.utcnow,
#        onupdate=datetime.datetime.utcnow, nullable=False)
#        
#    
#    ### Relationships
#    model = sa.orm.relationship(
#        UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
#        backref=sa.orm.backref('mutations', cascade='expunge')) # many to one
#
