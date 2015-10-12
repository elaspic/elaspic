# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import object

import os
import os.path as op
import re
import time
import subprocess
import tempfile
import atexit
import signal
import six
import logging 

from functools import wraps 

import pandas as pd

from Bio.PDB import PDBParser

from . import conf, errors, helper, structure, structure_analysis

logger = logging.getLogger(__name__)

sql_db = None
domain_alignment = None
domain_model = None
domain_mutation = None

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

from .database_tables import (
    Base, 
    Domain, DomainContact, UniprotSequence, Provean,
    UniprotDomain, UniprotDomainTemplate, UniprotDomainModel, UniprotDomainMutation,
    UniprotDomainPair, UniprotDomainPairTemplate, UniprotDomainPairModel, UniprotDomainPairMutation,
)


#%% 
class PipelineTemplate(object):

    def __init__(self, configurations):
        """
        It should be possible to initialize one pipeline and call it in parallel using different
        mutations as input
        """
        global database, sequence, model, mutation
        
        # Read the configuration file and set the variables
        if isinstance(configurations, six.string_types):
            conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            conf.configs = configurations.copy()
        
        # Can imort sql_db only after the configuration file has been read
        # Otherwise classes in `sql_db` won't be properly configured
        from . import database, sequence, model, mutation

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
        no_job_specific_folder = (
            configs['temp_path'].startswith(
                os.path.join(configs['global_temp_path'], configs['temp_path_suffix']))
        )
        on_node_with_manditory_job_specific_folder = (
            any([(x.lower() in hostname) for x in ['node', 'behemoth', 'grendel', 'beagle']])
        )
        if no_job_specific_folder and on_node_with_manditory_job_specific_folder:
            raise Exception('You should be using a temp folder that it specific to the particular job!')
    

    def run(self):
        raise NotImplementedError()
                
    def make_sequence(self):
        raise NotImplementedError()
        
    def mutate_sequence(self):
        raise NotImplementedError()
        
    def make_model(self):
        raise NotImplementedError()
        
    def mutate_model(self):
        raise NotImplementedError()
        
    def predict_ddg(self, protein_id, mutation, sequence, model):
        raise NotImplementedError()






#%%
class Pipeline(PipelineTemplate):

    def __init__(self, uniprot_id, mutations, configurations, run_type=1, number_of_tries=[], uniprot_domain_pair_ids=[]):
        """Run the main function of the program and parse errors.
        
        Parameters
        ----------
        uniprot_domain_pair_ids : list of integers
            List of uniprot_domain_pair_ids specifying which uniprot domain pairs to analyse.
        """
        super().__init__(configurations)
        
        self.uniprot_id = uniprot_id
        self.mutations = mutations.split(',')
        self.calculated_mutations = []
        self.run_type = run_type
        self.number_of_tries = number_of_tries
        self.uniprot_domain_pair_ids = uniprot_domain_pair_ids
        
        # TODO: change pipeline so that `self.unique_temp_folder` does not end in '/'            
        logger.info('uniprot_id: {}'.format(uniprot_id))
        logger.info('mutations: {}'.format(mutations))
        logger.info('run_type: {}'.format(run_type))
        logger.info('uniprot_domain_pair_ids: {}'.format(uniprot_domain_pair_ids))
        logger.info('unique_temp_folder: {}'.format(conf.configs['unique_temp_folder']))
        logger.info('db_schema: {}'.format(conf.configs.get('db_schema')))
        logger.info('global_temp_path: {global_temp_path}'.format(**conf.configs))
        logger.info('temp_path: {temp_path}'.format(**conf.configs))

        # Switch to the root of the unique tmp directory
        os.chdir(conf.configs['unique_temp_folder'])

        # Initialise the sql database for accessing all information
        self.db = database.MyDatabase()
        
        # Obtain all domains and interactions for a given uniprot
        if conf.configs['look_for_interactions'] == 2:
            logger.info('Skipping protein domain information...')
            self.uniprot_domains = []
        else:
            logger.info('Obtaining protein domain information...')
            self.uniprot_domains = self.db.get_uniprot_domain(self.uniprot_id, True)
            self._update_path_to_data(self.uniprot_domains)

        # Mutations results
        self.uniprot_mutations = []


    def _update_path_to_data(self, d_list):
        if not isinstance(d_list, list):
            d_list = [d_list]
        for d in d_list:
            if not d.path_to_data or any([len(x) > 255 for x in d.path_to_data.split('/')]):
                d.path_to_data = database.get_uniprot_base_path(d) + database.get_uniprot_domain_path(d)
                self.db.merge_row(d)
            os.makedirs(op.join(conf.configs['temp_archive_path'], d.path_to_data), exist_ok=True)


    def run(self):
        # Find provean
        if not self.uniprot_domains:
            logger.info('Warning! Uniprot {} has no pfam domains!'.format(self.uniprot_id))
        # You don't want to compute provean for a protein wihtout domains,
        # but you might want to evaluate mutations (if everything else is 
        # precalculated and you are running on a cluster).
        if self.run_type in [1, 5] and self.uniprot_domains:
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Computing provean...")
            provean = self.calculate_provean()
            if provean:
                if self.run_type == 1:
                    logger.info('Finished run_type {}'.format(self.run_type))
                    return
                # If provean was updated, we need to reload uniprot domain data
                # for all the other domains
                logger.info('\n\n\n')
                logger.info('Obtaining protein domain information...')
                self.uniprot_domains = self.db.get_uniprot_domain(self.uniprot_id, True)

        # Get interactions
        if conf.configs['look_for_interactions']:
            logger.info('Obtaining protein domain pair information...')
            self.uniprot_domain_pairs = self.db.get_uniprot_domain_pair(
                self.uniprot_id, True, self.uniprot_domain_pair_ids)
            self._update_path_to_data(self.uniprot_domain_pairs)

        # Make models
        if self.run_type in [2, 4, 5]:
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Building models...")
            for d in self.uniprot_domains + self.uniprot_domain_pairs:
                self.calculate_model()
            logger.info('Finished processing all models for {} {}'.format(self.uniprot_id, self.mutations))

        # Analyse mutations
        if self.run_type in [3, 4, 5] and self.mutations:
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Analyzing mutations...")
            for d in self.uniprot_domains + self.uniprot_domain_pairs:
                for mutation in self.mutations:
                    self.calculate_mutation(d, mutation)


    def calculate_provean(self):
        """
        """
        d = self.uniprot_domains[0] # All uniprot domains will have the same uniprot_sequence info
        self.__print_header(d)

        if (d.uniprot_sequence.provean and
            d.uniprot_sequence.provean.provean_supset_filename):
                path_to_provean_supset = (
                    conf.configs['temp_archive_path'] + database.get_uniprot_base_path(d) +
                    d.uniprot_sequence.provean.provean_supset_filename)
#                path_to_provean_supset = (
#                    conf.configs['path_to_archive'] + database.get_uniprot_base_path(d) +
#                    d.uniprot_sequence.provean.provean_supset_filename)
                if os.path.isfile(path_to_provean_supset):
                    if not conf.configs['remake_provean_supset']:
                        logger.debug('The provean supset has already been calculated. Done!\n')
                        return None
                    first_aa = d.uniprot_sequence.uniprot_sequence[0]
                    domain_mutation = '{0}1{0}'.format(first_aa)
                    result, error_message, return_code = \
                        domain_alignment.check_provean_supporting_set(
                            domain_mutation, d.uniprot_sequence.uniprot_sequence,
                            conf.configs, d.uniprot_id, path_to_provean_supset,
                            save_supporting_set=False, check_mem_usage=False)
                    if return_code == 0:
                        logger.debug('The provean supset has already been calculated. Done!\n')
                        return None
                    logger.debug('Provean supporting set caused an error:\n\n{}'.format(error_message))
                    logger.debug('Recompiling...')
                    provean = d.uniprot_sequence.provean
                else:
                    logger.error(
                        'Provean has been calculated but the file is missing from:\n{}\nRecompiling...'
                        .format(path_to_provean_supset))
                    provean = d.uniprot_sequence.provean
        elif d.uniprot_sequence.provean:
            logger.info('Provean supporting set has not been calculated previously. Computing...')
            provean = d.uniprot_sequence.provean
        else:
            provean = database.Provean()
            provean.uniprot_id = d.uniprot_id

        try:
            path_to_provean_supset = (
                domain_alignment.get_path_to_provean_supset(
                    d.uniprot_id, d.uniprot_sequence.uniprot_name, conf.configs)
            )
            provean.provean_supset_filename, provean.provean_supset_length = \
                domain_alignment.build_provean_supporting_set(
                    d.uniprot_id,
                    d.uniprot_sequence.uniprot_name,
                    d.uniprot_sequence.uniprot_sequence,
                    conf.configs, 
                    path_to_provean_supset)

        except (errors.ProveanError, errors.ProveanResourceError) as e:
            provean.provean_errors = self.__add_new_error(provean.provean_errors, e)
            logger.error(provean.provean_errors)
            logger.error(e.__str__())
            if isinstance(e, errors.ProveanResourceError):
                # Send the kill signal to the main process group, killing everything
                self._clear_provean_temp_files() # these won't get cleaned once the process dies
                logger.error('Killing group...')
                os.killpg(e.child_process_group_id, signal.SIGTERM)
                logger.critical(
                    'Provean ran out of resources. Everything has been killed. '
                    'The user will probably not see this message...')

        self._clear_provean_temp_files()
        self.db.merge_provean(provean, conf.configs['copy_data'] and database.get_uniprot_base_path(d))
        logger.info('Finished computing provean for {}\n'.format(d.uniprot_id))
        return provean


    def calculate_model(self, d):
        """
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        possible_model_errors = (
            errors.ModellerError,
            errors.PDBChainError,
            errors.PDBEmptySequenceError,
            errors.PDBNotFoundError,
            errors.ChainsNotInteractingError,
            errors.MutationOutsideDomainError,
            errors.MutationOutsideInterfaceError,
            errors.NoSequenceFound,
            errors.TcoffeeError,
        )

        # Go over all domains and domain pairs for a given protein
        self.__print_header(d)

        # Check if we should skip the model
        if (isinstance(d, database.UniprotDomain) and
                not (d.template and d.template.cath_id) ):
            logger.error('No structural template availible for this domain. Skipping...')
            return
        elif (isinstance(d, database.UniprotDomainPair) and
                not (d.template and d.template.cath_id_1 and d.template.cath_id_2) ):
            logger.error('No structural template availible for this domain pair. Skipping...')
            return
        elif d.template.model and d.template.model.model_filename:
            logger.info('Model already calculated. Skipping...')
            return
        elif (d.template.model and d.template.model.model_errors and
                (('Giving up' in d.template.model.model_errors) or
                (d.template.model.model_errors.count(';') > 10)) ):
            logger.info(
                'Previous model had unfixable errors: "{}". Skipping...'
                .format(d.template.model.model_errors))
            return

        # Make a homology model using the domain and template information
        try:
            self.model(d)

        except possible_model_errors as e:
            # Find domains that were used as a template and eventually led to
            # the error in the model, and add the error to their `domain_errors`
            # or `domain_contact_errors` fields.
            logger.error(str(e))
            if isinstance(d, database.UniprotDomain):
                if d.template.model == None:
                    d.template.model = sql_db.UniprotDomainModel()
                    d.template.model.uniprot_domain_id = d.uniprot_domain_id
                bad_domains = self.db.get_rows_by_ids(
                    sql_db.Domain, [sql_db.Domain.cath_id], [d.template.cath_id])
                assert len(bad_domains) == 1
                bad_domain = bad_domains[0]
                bad_domain.domain_errors = str(d.uniprot_domain_id) + ': ' + str(type(e))
                logger.error(
                    "Making a homology model failed!!!\n"
                    "Adding error '{0}' to the domain with cath_id {1}..."
                    .format(bad_domain.domain_errors, d.template.cath_id))
            elif isinstance(d, sql_db.UniprotDomainPair):
                if d.template.model == None:
                    d.template.model = sql_db.UniprotDomainPairModel()
                    d.template.model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
                bad_domains = self.db.get_rows_by_ids(
                    sql_db.DomainContact,
                    [sql_db.DomainContact.cath_id_1, sql_db.DomainContact.cath_id_2],
                    [d.template.cath_id_1, d.template.cath_id_2])
                if len(bad_domains) == 0:
                    bad_domains = self.db.get_rows_by_ids(
                        sql_db.DomainContact,
                        [sql_db.DomainContact.cath_id_1, sql_db.DomainContact.cath_id_2],
                        [d.template.cath_id_2, d.template.cath_id_1])
                assert len(bad_domains) == 1
                bad_domain = bad_domains[0]
                bad_domain.domain_contact_errors = str(d.uniprot_domain_pair_id) + ': ' + str(type(e))
                logger.error(
                    "Making a homology model failed!!!\n"
                    "Adding error '{0}' to the domain pair with cath_id_1 {1} "
                    "and cath_id_2 {2}..."
                    .format(bad_domain.domain_contact_errors, d.template.cath_id_1, d.template.cath_id_2))
            # Add the error type to the model_errors column
            d.template.model.model_errors = self.__add_new_error(d.template.model.model_errors, e)
            logger.error(d.template.model.model_errors)
            self.db.merge_row(bad_domain)
            # d.template.model = empty_model

        # Add either the empty model or the calculated model to the database
        logger.info('Adding model...')
        self.db.merge_model(d, conf.configs['copy_data'] and d.path_to_data)
        

    def get_mutation_data(self, d, uniprot_id_1, mutation):
        """
        """
        path_to_provean_supset = None
        #######################################################################
        # Core
        if isinstance(d, database.UniprotDomain):
            d_1 = d
            logger.debug("Analyzing core mutation for uniprot: %s" % uniprot_id_1)
            pdbfam_name = d.pdbfam_name
            if d.template.model.model_domain_def != None:
                logger.debug('Using model domain definitions')
                domain_start, domain_end = helper.decode_domain_def(d.template.model.model_domain_def)
            else:
                logger.debug('Using template domain definitions')
                domain_start, domain_end = helper.decode_domain_def(d.template.domain_def)
            alignment, __ = self.db.get_alignment(d.template.model, d.path_to_data)
            chains_pdb = [d.template.domain.pdb_chain, ]
            uniprot_sequences = [self.db.get_uniprot_sequence(uniprot_id_1), ]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],]
            chains_modeller = [d.template.model.chain, ]
            uniprot_domain_id = d.uniprot_domain_id
    
        #######################################################################
        # Interface
        elif isinstance(d, database.UniprotDomainPair):
            mutation_pos = int(mutation[1:-1])
            interacting_aa_1 = self._get_interacting_aa(d, 1)
            interacting_aa_2 = self._get_interacting_aa(d, 2)
            assert uniprot_id_1 in [d.uniprot_domain_1.uniprot_id, d.uniprot_domain_2.uniprot_id]
            
            if uniprot_id_1 == d.uniprot_domain_1.uniprot_id and mutation_pos in interacting_aa_1:
                # Mutation is inside the first domain
                uniprot_id_2 = d.uniprot_domain_2.uniprot_id
                d_1, d_2 = d.uniprot_domain_1, d.uniprot_domain_2
                #
                if (d.template.model.model_domain_def_1 != None and
                    d.template.model.model_domain_def_2 != None):
                        domain_start, domain_end = helper.decode_domain_def(d.template.model.model_domain_def_1)
                        domain_2_start, domain_2_end = helper.decode_domain_def(d.template.model.model_domain_def_2)
                else:
                        domain_start, domain_end = helper.decode_domain_def(d.uniprot_domain_1.template.domain_def)
                        domain_2_start, domain_2_end = helper.decode_domain_def(d.uniprot_domain_2.template.domain_def)
    
                alignment, __ = self.db.get_alignment(d.template.model, d.path_to_data)
                chains_pdb = [d.template.domain_1.pdb_chain, d.template.domain_2.pdb_chain]
                chains_modeller = [d.template.model.chain_1, d.template.model.chain_2]
    
            elif uniprot_id_1 == d.uniprot_domain_2.uniprot_id and mutation_pos in interacting_aa_2:
                # Mutation is inside the second domain
                logger.debug('Mutated uniprot is uniprot 2. Rearranging...')
                uniprot_id_2 = d.uniprot_domain_1.uniprot_id
                d_1, d_2 = d.uniprot_domain_2, d.uniprot_domain_1
    
                if (d.template.model.model_domain_def_1 != None and
                    d.template.model.model_domain_def_2 != None):
                        domain_start, domain_end = helper.decode_domain_def(d.template.model.model_domain_def_2)
                        domain_2_start, domain_2_end = helper.decode_domain_def(d.template.model.model_domain_def_1)
                else:
                        domain_start, domain_end = helper.decode_domain_def(d.uniprot_domain_2.template.domain_def)
                        domain_2_start, domain_2_end = helper.decode_domain_def(d.uniprot_domain_1.template.domain_def)
    
                __, alignment = self.db.get_alignment(d.template.model, d.path_to_data)
                chains_pdb = [d.template.domain_2.pdb_chain, d.template.domain_1.pdb_chain]
                chains_modeller = [d.template.model.chain_2, d.template.model.chain_1]
    
            else:
                # Mutation is outside the interface
                logger.error('Uniprot ID: {}\tMutation: {}'.format(
                    uniprot_id_1, mutation))
                logger.error('Uniprot ID 1: {}\tInteracting AA 1: {}'.format(
                    d.uniprot_domain_1.uniprot_id, interacting_aa_1))
                logger.error('Uniprot ID 2: {}\tInteracting AA 2: {}'.format(
                    d.uniprot_domain_2.uniprot_id, interacting_aa_2))
                raise errors.MutationOutsideInterfaceError('mutated residue not involved in the interaction')
    
    
            logger.debug('Analysing interface mutation between uniprots %s and %s' % (uniprot_id_1, uniprot_id_2,))
            uniprot_sequences = [self.db.get_uniprot_sequence(d_1.uniprot_id),
                                 self.db.get_uniprot_sequence(d_2.uniprot_id)]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],
                                uniprot_sequences[1][domain_2_start-1:domain_2_end]]
    
    
        #######################################################################
        # Common
        if int(mutation[1:-1]) < domain_start or int(mutation[1:-1]) > domain_end:
            raise errors.MutationOutsideDomainError('Mutation falls outside domain')
    
        save_path = conf.configs['temp_archive_path'] + d.path_to_data
        pdbFile_wt = d.template.model.model_filename
        pdbfam_name = d_1.pdbfam_name
        uniprot_domain_id = d_1.uniprot_domain_id
    
        # Provean
        if (d_1.uniprot_sequence.provean and
            d_1.uniprot_sequence.provean.provean_supset_filename):
                path_to_provean_supset = (
                    conf.configs['temp_archive_path'] + database.get_uniprot_base_path(d_1) +
                    d_1.uniprot_sequence.provean.provean_supset_filename )
                if not os.path.isfile(path_to_provean_supset):
                    error_message = (
                        'Provean supporting set sequence does not exist even though it should!\n{}'
                        .format(path_to_provean_supset)
                    )
                    logger.error(error_message)
                    logger.error('d_1: {}'.format(d_1))
                    logger.error('listdir: {}'.format(os.listdir(os.path.dirname(path_to_provean_supset))))
                    raise Exception(error_message)
    
    
        # Convert mutation from uniprot domain numbering to pdb domain numbering
        position_domain = int(mutation[1:-1]) - domain_start + 1 # +1 to have the numbering staring with 1
        mutation_domain = mutation[0] + str(position_domain) + mutation[-1]
    
        logger.debug('Modeller pdb file: {}'.format(save_path + pdbFile_wt))
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        structure = parser.get_structure('ID', save_path + pdbFile_wt)
        
        def convert_position_to_resid(chain, positions, domain_def_tuple=None):
            __, chain_numbering = structure_tools.get_chain_sequence_and_numbering(
                chain, domain_def_tuple)
            return [chain_numbering[p-1] for p in positions]
        
        position_modeller = convert_position_to_resid(
            structure[0][chains_modeller[0]], [position_domain])
        mutation_modeller = mutation[0] + position_modeller[0] + mutation[-1]
    
        # Save the results
        mut_data = mutation.MutationData()
    
        mut_data.uniprot_id_1 = uniprot_id_1
        mut_data.mutation = mutation
        mut_data.pdbfam_name = pdbfam_name
        mut_data.domain_start = domain_start
        mut_data.domain_end = domain_end
        mut_data.alignment = alignment
        mut_data.chains_pdb = chains_pdb
        mut_data.uniprot_sequences = uniprot_sequences
        mut_data.domain_sequences = domain_sequences
        mut_data.chains_modeller = chains_modeller
        mut_data.uniprot_domain_id = uniprot_domain_id
        mut_data.path_to_provean_supset = path_to_provean_supset
        mut_data.save_path = save_path
        mut_data.pdbFile_wt = pdbFile_wt
        mut_data.position_domain = position_domain
        mut_data.mutation_domain = mutation_domain
    #        mut_data.structure = structure
        mut_data.position_modeller = position_modeller
        mut_data.mutation_modeller = mutation_modeller
    
        for key, value in mut_data.__dict__.items():
            logger.debug(key + ':')
            logger.debug(value)
    
        # Run some sanity checks
        mutated_aa_uniprot = mut_data.uniprot_sequences[0][int(mutation[1:-1])-1]
        if mutated_aa_uniprot != mutation[0]:
            logger.error('Uniprot sequence: {}'.format(uniprot_sequences[0]))
            logger.error('Uniprot AA: {};\t Mutation AA: {}'.format(mutated_aa_uniprot, mutation[0]))
            raise Exception('Mutated amino acid was not found inside the specified uniprot!')
    
        mutated_aa_domain = str(mut_data.domain_sequences[0].seq)[int(mut_data.mutation_domain[1:-1])-1]
        if mutated_aa_domain != mut_data.mutation_domain[0]:
            logger.error('Domain sequence: {}'.format(str(mut_data.domain_sequences[0].seq)))
            logger.error('Domain AA: {};\t Mutation AA: {}'.format(mutated_aa_domain, mut_data.mutation_domain[0]))
            raise Exception('Mutated amino acid was not found inside the specified domain!')
    
        return mut_data            
     
     
    def _get_interacting_aa(self, d, domain_1or2=1):
        if domain_1or2 == 1:
            interacting_aa = d.template.model.interacting_aa_1
        elif domain_1or2 == 2:
            interacting_aa = d.template.model.interacting_aa_2
        else:
            raise Exception("`domain_1or2` should be either '1' or '2'!")
        return [int(uniprot_num) for uniprot_num in interacting_aa.split(',') if uniprot_num]
        

    def calculate_mutation(self, d, mutation):
        """
        """
        possible_mutation_errors = (
            errors.PDBError,
            errors.FoldxError,
            errors.ResourceError,
            errors.FoldXAAMismatchError,
        )
        self.get_mutation = domain_mutation.GetMutation(self.db)
        self.__print_header(d)

        if not self.mutations:
            logger.debug('Not evaluating mutations because no mutations specified...')
            return
        elif ( (isinstance(d, sql_db.UniprotDomain) and not d.template.domain) or
                (isinstance(d, sql_db.UniprotDomainPair) and not (d.template.domain_1 and d.template.domain_2)) ):
            logger.debug('Skipping because no structural template is availible...')
            return
        elif d.template.model == None or d.template.model.model_filename == None:
            logger.debug('d.template.model: {}'.format(d.template.model))
            logger.debug('d.template.model.model_filename: {}'.format(
                getattr(d.template.model, 'model_filename', 'does not exist!')))
            logger.debug('Skipping because no model...')
            return
        elif d.template.model.model_errors != None:
            logger.debug('Skipping because the model has errors: {}!'.format(d.template.model.model_errors))
            return

        logger.debug('Going over all mutations: {}'.format(self.mutations))

        logger.debug('-' * 80)
        logger.debug('Mutation: {}'.format(mutation))
        mutation_prototype = re.compile("^[A-z][0-9]+[A-z]$")
        if not mutation_prototype.match(mutation) or int(mutation[1:-1]) == 0:
            logger.error('The mutation {} is not a supported format! Skiping!'.format(mutation))
            return
        # Check to see if we have a precalculated mutation. Skip if all
        # parameters have been calculated; otherwise analyse the remaining
        # parameters. Create an empty mutation if the mutation has not
        # been precalculated.
        precalculated_mutation = self.db.get_uniprot_mutation(d, mutation, self.uniprot_id, True)
        logger.info('Have the following precalculated mutation: {}'.format(precalculated_mutation))
        if (precalculated_mutation and
                (precalculated_mutation.provean_score and
                precalculated_mutation.stability_energy_wt and
                precalculated_mutation.ddg != None)): 
            self.calculated_mutations.append(precalculated_mutation)
            logger.info('Mutation has already been completely evaluated. Skipping...')
            return
        elif precalculated_mutation:
            uniprot_mutation = precalculated_mutation
        # Construct empty models that will be used if we have errors
        elif isinstance(d, sql_db.UniprotDomain):
            uniprot_mutation = sql_db.UniprotDomainMutation()
            uniprot_mutation.uniprot_domain_id = d.uniprot_domain_id
            uniprot_mutation.uniprot_id = self.uniprot_id
            uniprot_mutation.mutation = mutation
        elif isinstance(d, sql_db.UniprotDomainPair):
            uniprot_mutation = sql_db.UniprotDomainPairMutation()
            uniprot_mutation.uniprot_domain_pair_id = d.uniprot_domain_pair_id
            uniprot_mutation.uniprot_id = self.uniprot_id
            uniprot_mutation.mutation = mutation

        try:
            mut_data = self.get_mutation.get_mutation_data(d, self.uniprot_id, mutation)
            uniprot_mutation = self.get_mutation.evaluate_mutation(d, mut_data, uniprot_mutation)
        except (errors.MutationOutsideDomainError,
                errors.MutationOutsideInterfaceError) as e:
            logger.debug('{}: {}; OK'.format(type(e), e))
            return
        except possible_mutation_errors as e:
            uniprot_mutation.mutation_errors = '{}: {}'.format(type(e), e)
            logger.debug(uniprot_mutation.mutation_errors)

        logger.info('Adding mutation {}'.format(mutation))
        logger.debug("Mutation attributes:")
        for attr in dir(uniprot_mutation):
            if attr.startswith('_'):
                return
            attr_field = getattr(uniprot_mutation, attr)
            attr_type = type(attr_field)
            logger.debug(attr)
            logger.debug(attr_field)
            logger.debug(attr_type)
            if (six.PY2 and
                (isinstance(attr_field, six.binary_type) or
                 isinstance(attr_field, six.text_type))):
                    logger.debug(
                         'Changing attribute {} from {} to {}...'
                         .format(attr, attr_type, str(attr_field)))
                    setattr(uniprot_mutation, attr, str(attr_field))
            if six.PY3 and isinstance(attr_field, six.binary_type):
                logger.debug(
                    'Changing attribute {} from {} to {}...'
                    .format(attr, attr_type, type(attr_field.decode())))
                setattr(uniprot_mutation, attr, attr_field.decode())
        self.calculated_mutations.append(uniprot_mutation)
        self.db.merge_mutation(uniprot_mutation, conf.configs['copy_data'] and d.path_to_data)
        self.uniprot_mutations.append(uniprot_mutation)
        logger.info('Finished processing all mutations for {} {}'.format(self.uniprot_id, self.mutations))


    def get_mutation_features(self, d, mut, row_idx=0):
        """
        Returns a dataframe that contains all features for the given mutation that are relevant for
        machine learning.
    
        Parameters
        ----------
        d : sqlalchemy orm object
           Contains domain information for the given mutation
        mut : sqlalchemy orm object
           Contains information about the mutation in question
        """
        feature_dict = {key: value for (key, value) in mut.__dict__.items() if not key.startswith('_')}
    
        feature_dict.update({
            # Header columns
            # 'uniprot_id': mut.uniprot_id,
            # 'mutation': mut.mutation,
            't_date_modified': d.template.t_date_modified,
            'm_date_modified': d.template.model.m_date_modified,
            # 'mut_date_modified': mut.mut_date_modified,
            #
            'norm_dope': d.template.model.norm_dope,
        })
    
        if hasattr(d, 'uniprot_domain_id'):
            feature_dict.update({
                # Header columns
                'uniprot_domain_id': d.uniprot_domain_id,
                'cath_id': d.template.cath_id,
                'pfam_name': d.pdbfam_name,
                'clan_name': d.pfam_clan,
                #
                'alignment_identity': d.template.alignment_identity,
                'alignment_coverage': d.template.alignment_coverage,
                'alignment_score': d.template.alignment_score,
            })
    
        elif hasattr(d, 'uniprot_domain_pair_id'):
            feature_dict.update({
                # Header columns
                'uniprot_domain_pair_id': d.uniprot_domain_pair_id,
                'cath_id_1': d.template.cath_id_1,
                'cath_id_2': d.template.cath_id_2,
                'pfam_name_1': d.uniprot_domain_1.pdbfam_name,
                'pfam_name_2': d.uniprot_domain_2.pdbfam_name,
                'clan_name_1': d.uniprot_domain_1.pfam_clan,
                'clan_name_2': d.uniprot_domain_2.pfam_clan,
                # Feature columns
                'alignment_identity': d.template.identical_1 + d.template.identical_2,
                'identical_1': d.template.identical_1,
                'identical_2': d.template.identical_2,
                'interface_area_hydrophobic': d.template.model.interface_area_hydrophobic,
                'interface_area_hydrophilic': d.template.model.interface_area_hydrophilic,
                'interface_area_total': d.template.model.interface_area_total,
                'score_total': d.template.score_total,
            })
    
        feature_df = pd.DataFrame(feature_dict, index=[row_idx])
    
        return feature_df


    def __get_unique_id(self, d):
        if isinstance(d, sql_db.UniprotDomain):
            return ('uniprot_domain_id', d.uniprot_domain_id)
        else:
            return('uniprot_domain_pair_id', d.uniprot_domain_pair_id)


    def __print_header(self, d):
        # logger.info('Domain or domain pair number: {}'.format(d_idx))
        logger.info('=' * 77)
        unique_id = self.__get_unique_id(d)
        logger.info('{}: {}'.format(*unique_id))


    def __add_new_error(self, d_error_log, e):
        if d_error_log is None:
            return str(type(e))
        else:
            return '{}; {}: {}'.format(d_error_log, type(e), str(e))



#%%
def lock(fn):
    """
    Allow only a single instance of function `fn`, 
    and save results to a lock file.
    """
    @wraps(fn)
    def locked_fn(self, *args, **kwargs):
        """
        
        Returns
        -------
        lock_filename : str
            Lock file that contains function output in json format.
        
        """
        # Get the lock filename
        if fn.__name__ == 'calculate_provean':
            lock_filename = '{}{}_provean.json'.format(self.pdb_id, args[0])
        elif fn.__name__ == 'calculate_model':
            lock_filename = '{}_modeller.json'.format(self.pdb_id)
        elif fn.__name__ == 'calculate_mutation':
            lock_filename = '{}{}_mutation_{}.json'.format(self.pdb_id, args[0], args[1])
        else:
            raise Exception("Function {} is not supported!".format(fn))
        
        # Make sure that we can get exclusive rights on the lock
        try:
            lock = open(lock_filename, 'x')
        except FileExistsError:
            try:
                results = json.load(open(lock_filename, 'r'))
                info_message = (
                    "Results have already been calculated and are in file: '{}'.\n"
                    .format(lock_filename, results)
                )
                logger.info(info_message)
                return lock_filename, results
            except ValueError:
                info_message = (
                    "Another process is currently running this function.\n"
                    "If you believe this is an error, delete lock file '{}' and try again."
                    .format(lock_filename)
                )
                logger.info(info_message)
                return lock_filename, None
        
        # Run the function and write results
        try:
            results = fn(self, *args, **kwargs)
            json.dump(results, lock)
            lock.close()
            return lock_filename, results
        except:
            lock.close()
            os.remove(lock.name)
            raise

    return locked_fn


class PipelineStructure(PipelineTemplate):

    def __init__(self, pdb_file, pdb_mutations, configurations):
        """
        TODO: Add an option to store provean results based on sequence hash.
        """
        super().__init__(configurations)
        
        # Input parameters
        self.R_CUTOFF = 6

        self.pdb_id = op.basename(pdb_file).rstrip('.pdb')
        self.pdb_file = pdb_file
        self.mutations = pdb_mutations.split(',')
        
        logger.info('pdb_file: {}'.format(self.pdb_file))
        logger.info('pwd: {}'.format(self.PWD))
        
        ### Load PDB structure and extract required sequences and chains.
        #fix_pdb(self.pdb_file, self.pdb_file)
        self.sp = structure.StructureParser(self.pdb_file)
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
        

    @lock
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
        if last_chain.child_list[0].resname in structure.AAA_DICT:
            chain_sequence, chain_numbering_extended = self.sp.get_chain_sequence_and_numbering(last_chain.id)
            structure_sequence += [chain_sequence]
        else:
            structure_sequence += ['.' * len(last_chain.child_list)]
            
        structure_sequence = '/'.join(structure_sequence)
        return structure_sequence
    
    
    @lock
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
        

    #@lock
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
            residue_id = (structure.AAA_DICT[residue.resname] + str(residue.id[1]))
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
class UserStructureMutation(Base):
    __tablename__ = 'user_structure_mutation'
    _indexes = [
        ['pdb_id', 'pdb_chain', 'pdb_chain_2', 'pdb_mutation'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
    
    id = sa.Column(sa.Integer, nullable=False, primary_key=True, autoincrement=True)
    
    ### Input
    unique = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')), 
        primary_key=True)
    pdb_file = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
    pdb_id = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
    pdb_chains = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
    
    
    ### Mutation    
    # Original PDB coordinates
    pdb_mutation = sa.Column(sa.String(SHORT), nullable=False)
    # Where 1 is the start of the domain (i.e. chain)
    domain_mutation = sa.Column(sa.String(SHORT), nullable=False)
    # Where 1 is the start of the protein (i.e. chain A residue 1)
    modeller_mutation = sa.Column(sa.String(SHORT), nullable=False)
    
    
    ### Homology Model
    # Core / Interface
    model_filename = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
    norm_dope = sa.Column(sa.Float)
    sasa_score = sa.Column(sa.Text)
        
    # Interface
    pdb_chain = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
    pdb_chain_2 = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
    interface_area_hydrophobic = sa.Column(sa.Float)
    interface_area_hydrophilic = sa.Column(sa.Float)
    interface_area_total = sa.Column(sa.Float)
    interacting_aa_1 = sa.Column(sa.Text)
    interacting_aa_2 = sa.Column(sa.Text)
    
    
    
    ### Mutation
    mutation = sa.Column(sa.String(SHORT), nullable=False)
    mutation_errors = sa.Column(sa.Text)
    model_filename_wt = sa.Column(sa.String(MEDIUM))
    model_filename_mut = sa.Column(sa.String(MEDIUM))
    chain_modeller = sa.Column(sa.String(SHORT))
    mutation_modeller = sa.Column(sa.String(SHORT))
    analyse_complex_energy_wt = sa.Column(sa.Text)
    stability_energy_wt = sa.Column(sa.Text)
    analyse_complex_energy_mut = sa.Column(sa.Text)
    stability_energy_mut = sa.Column(sa.Text)
    physchem_wt = sa.Column(sa.Text)
    physchem_wt_ownchain = sa.Column(sa.Text)
    physchem_mut = sa.Column(sa.Text)
    physchem_mut_ownchain = sa.Column(sa.Text)
    matrix_score = sa.Column(sa.Float)
    secondary_structure_wt = sa.Column(sa.Text)
    solvent_accessibility_wt = sa.Column(sa.Float)
    secondary_structure_mut = sa.Column(sa.Text)
    solvent_accessibility_mut = sa.Column(sa.Float)
    contact_distance_wt = sa.Column(sa.Float)
    contact_distance_mut = sa.Column(sa.Float)
    provean_score = sa.Column(sa.Float)
    ddg = sa.Column(sa.Float, index=False)
    mut_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow,
        onupdate=datetime.datetime.utcnow, nullable=False)
        
    
    ### Relationships
    model = sa.orm.relationship(
        UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('mutations', cascade='expunge')) # many to one

