# -*- coding: utf-8 -*-
import os
import os.path as op
import re
import json
import shutil
import logging

import six
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser

from . import (
    conf, errors, helper, structure_tools, sequence, model, predictor, 
    database, database_tables
)
from .pipeline import Pipeline, execute_and_remember

logger = logging.getLogger(__name__)
configs = conf.Configs()



#%%
class DatabasePipeline(Pipeline):

    def __init__(
            self, uniprot_id, mutations, configurations=None, 
            run_type=5, number_of_tries=[], uniprot_domain_pair_ids=[]):
        """Run the main function of the program and parse errors.
        
        Parameters
        ----------
        uniprot_domain_pair_ids : list of integers
            List of uniprot_domain_pair_ids specifying which uniprot domain pairs to analyse.
        """
        super().__init__(configurations)
        
        self.uniprot_id = uniprot_id
        self.calculated_mutations = []
        self.run_type = run_type
        self.number_of_tries = number_of_tries
        self.uniprot_domain_pair_ids = uniprot_domain_pair_ids
        
        if mutations is None:
            self.mutations = []
        elif isinstance(mutations, str):
            self.mutations = mutations.split(',')
        else:
            self.mutations = mutations
        
        logger.info(helper.color.UNDERLINE + 'Input parameters:' + helper.color.END)
        logger.info('uniprot_id: {}'.format(uniprot_id))
        logger.info('mutations: {}'.format(mutations))
        logger.info('run_type: {}'.format(run_type))
        logger.info('uniprot_domain_pair_ids: {}'.format(uniprot_domain_pair_ids))
        logger.info('unique_temp_dir: {}'.format(configs['unique_temp_dir']))
        logger.info('db_schema: {}'.format(conf.configs.get('db_schema')))
        logger.info('global_temp_dir: {global_temp_dir}'.format(**conf.configs))
        logger.info('temp_dir: {temp_dir}'.format(**conf.configs))
        
        # Switch to the root of the unique tmp directory
        os.chdir(configs['unique_temp_dir'])

        # Initialise the sql database for accessing all information
        self.db = database.MyDatabase()
        
        # Obtain all domains and interactions for a given uniprot
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
                d.path_to_data = (
                    database.get_uniprot_base_path(d) + database.get_uniprot_domain_path(d)
                )
                self.db.merge_row(d)
            os.makedirs(op.join(configs['archive_temp_dir'], d.path_to_data), exist_ok=True)


    def run(self):
        if not self.uniprot_domains:
            logger.info('Warning! Uniprot {} has no pfam domains!'.format(self.uniprot_id))
        
        # Find provean
        if self.run_type in [1, 5] and self.uniprot_domains and not self.sequences:
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Computing provean...")
            self.get_sequence(self.uniprot_domains[0])

        # Get interactions
        if configs['look_for_interactions']:
            logger.info('Obtaining protein domain pair information...')
            self.uniprot_domain_pairs = self.db.get_uniprot_domain_pair(
                self.uniprot_id, True, self.uniprot_domain_pair_ids)
            self._update_path_to_data(self.uniprot_domain_pairs)

        # Make models
        if self.run_type in [2, 4, 5]:
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Building models...")
            for d in self.uniprot_domains + self.uniprot_domain_pairs:
                self.get_model(d)
            logger.info(
                'Finished processing all models for {} {}'
                .format(self.uniprot_id, self.mutations)
            )

        # Analyse mutations
        if self.run_type in [3, 4, 5] and self.mutations:
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Analyzing mutations...")
            for d in self.uniprot_domains + self.uniprot_domain_pairs:
                for mutation in self.mutations:
                    self.get_mutation_score(d, mutation)


    def get_sequence(self, d):
        """
        """
        logger.debug('-' * 80)
        logger.debug('get_sequence({})'.format(d))
        return PrepareSequence(d, self.db)


    def get_model(self, d):
        """
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        logger.debug('-' * 80)
        logger.debug('get_model({})'.format(d))
        return PrepareModel(d, self.db)


    def get_mutation_score(self, d, mutation):
        logger.debug('-' * 80)
        logger.debug('get_mutation_score({}, {})'.format(d, mutation))
        if isinstance(d, database_tables.UniprotDomain):
            sequence = self.get_sequence(d)
        elif isinstance(d, database_tables.UniprotDomainPair):
            if self.uniprot_id == d.uniprot_domain_1.uniprot_id:
                sequence = self.get_sequence(d.uniprot_domain_1)
            elif self.uniprot_id == d.uniprot_domain_2.uniprot_id:
                sequence = self.get_sequence(d.uniprot_domain_2)
            else:
                raise Exception()
        else:
            raise Exception()
        model = self.get_model(d)
        return PrepareMutation(d, mutation, self.uniprot_id, sequence, model, self.db)



#%%
import random

@execute_and_remember
class Foo:
    def __init__(self):
        self.info = ['init {}'.format(random.randint(0,1000))]
        print(self.info)
        
    def __bool__(self):
        self.info.append('__bool__')
        print(self.info)
        return True

    def __enter__(self):
        self.info.append('__bool__')
        print(self.info)
        
    def run(self):
        self.info.append('__bool__')
        print(self.info)
        
    def __exit__(self, *exc):
        self.info.append('__bool__')
        print(self.info)
        print(exc)
        return False
            

#%%
@execute_and_remember
class PrepareSequence:
    
    def __init__(self, d, db):
        self.d = d
        self.db = db
        
        self.sequence = None
        self.skip = False

        if (d.uniprot_sequence and 
            d.uniprot_sequence.provean and
            d.uniprot_sequence.provean.provean_supset_filename):
            # Provean supset has already been calculated
            self.provean_supset_file = op.join(
                configs['archive_dir'],
                database.get_uniprot_base_path(d),
                d.uniprot_sequence.provean.provean_supset_filename
            )
        else:
            self.provean_supset_file = None


    def __bool__(self):
        return not self.skip
        
        
    def __enter__(self):
        d = self.d
        self.sequence_file = op.join(configs['sequence_dir'], d.uniprot_id + '.fasta')
        self.seqrecord = SeqRecord(
            id = d.uniprot_id,
            seq = Seq(d.uniprot_sequence.uniprot_sequence)
        )
        with open(self.sequence_file, 'w') as ofh:
            SeqIO.write(self.seqrecord, ofh, 'fasta')
        assert op.isfile(self.sequence_file)
    
    
    def run(self):
        self.sequence = sequence.Sequence(self.sequence_file, self.provean_supset_file)
#            if provean:
#                if self.run_type == 1:
#                    logger.info('Finished run_type {}'.format(self.run_type))
#                    return
#                # If provean was updated, we need to reload uniprot domain data
#                # for all the other domains
#                logger.info('\n\n\n')
#                logger.info('Obtaining protein domain information...')
#                self.uniprot_domains = self.db.get_uniprot_domain(self.uniprot_id, True)

        
    def __exit__(self, exc_type, exc_value, traceback):
        d = self.d

        if exc_type is not None:
            return False
        
        provean = database.Provean()
        provean.uniprot_id = d.uniprot_id
        provean.provean_supset_filename = op.basename(self.sequence.provean_supset_file)
        provean.provean_supset_length = self.sequence.provean_supset_length
        self.db.merge_provean(
            provean, 
            self.sequence.provean_supset_file,
            configs['copy_data'] and database.get_uniprot_base_path(d)
        )

    
    @property
    def result(self):
        return self.sequence



#%%
@execute_and_remember
class PrepareModel:
    
    handled_errors = (
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

    
    def __init__(self, d, db):
        """
        """
        print_header(d)
        self.d = d
        self.db = db
        self.skip = False
        self.model = None
        self.modeller_results_file = None
        
        # Check if we should skip the model
        if (isinstance(d, database.UniprotDomain) and not (d.template and d.template.cath_id)):
            logger.error('No structural template availible for this domain. Skipping...')
            self.skip = True
        elif (isinstance(d, database.UniprotDomainPair) and not 
                (d.template and d.template.cath_id_1 and d.template.cath_id_2)):
            logger.error('No structural template availible for this domain pair. Skipping...')
            self.skip = True
        elif d.template.model and d.template.model.model_filename:
            logger.info('Model already calculated.')
            self.modeller_results_file = self._get_modeller_results_file(d)
        elif (d.template.model and d.template.model.model_errors and
                (('Giving up' in d.template.model.model_errors) or
                (d.template.model.model_errors.count(';') > 10)) ):
            logger.info(
                'Previous model had unfixable errors: "{}". Skipping...'
                .format(d.template.model.model_errors))
            self.skip = True
            

    def _get_modeller_results_file(self, d):
        """Save relevant precalculated data to a json file.
        """
        modeller_results = dict()
        if isinstance(d, database.UniprotDomain):
            unique_id = d.uniprot_domain_id
            modeller_results['model_file'] = op.join(
                configs['archive_dir'],
                d.path_to_data,
                d.template.model.model_filename
            )
            modeller_results['alignment_files'] = [
                op.join(
                    configs['archive_dir'],
                    d.path_to_data,
                    d.template.model.alignment_filename,
                )
            ]
            modeller_results['norm_dope'] = d.template.model.norm_dope
        elif isinstance(d, database.UniprotDomainPair):
            unique_id = d.uniprot_domain_pair_id
            modeller_results['model_file'] = op.join(
                configs['archive_dir'],
                d.path_to_data,
                d.template.model.model_filename
            )
            modeller_results['alignment_files'] = [
                op.join(
                    configs['archive_dir'],
                    d.path_to_data,
                    d.template.model.alignment_filename_1,
                ),
                op.join(
                    configs['archive_dir'],
                    d.path_to_data,
                    d.template.model.alignment_filename_2,
                ),
            ]
            modeller_results['norm_dope'] = d.template.model.norm_dope

        modeller_results_file = op.join(
            configs['model_dir'],
            '{}_modeller_results.json'.format(unique_id)
        )
        with open(modeller_results_file, 'w') as ofh:
            json.dump(modeller_results, ofh)
        return modeller_results_file
        
        
    def __bool__(self):
        return not self.skip
    
    
    def __enter__(self):
        d = self.d
        ### Load data from database ORM to model input
        if isinstance(d, database.UniprotDomain):
            protein_ids = [d.uniprot_id]
            if d.template.model and d.template.model.model_domain_def:
                protein_domain_defs = [d.template.model.model_domain_def]
            else:
                protein_domain_defs = [d.template.domain_def]
            protein_sequences = [d.uniprot_sequence.uniprot_sequence]
            pdb_id = d.template.domain.pdb_id
            pdb_chains = [d.template.domain.pdb_chain]
            pdb_domain_defs = [d.template.domain.pdb_domain_def]

        elif isinstance(d, database.UniprotDomainPair):
            protein_ids = [
                d.uniprot_domain_1.uniprot_id,
                d.uniprot_domain_2.uniprot_id
            ]
            if (d.template.model and 
                d.template.model.model_domain_def_1 and 
                d.template.model.model_domain_def_2):
                protein_domain_defs = [
                    d.template.model.model_domain_def_1,
                    d.template.model.model_domain_def_2
                ]
            else:
                protein_domain_defs = [
                    d.uniprot_domain_1.template.domain_def,
                    d.uniprot_domain_2.template.domain_def
                ]
            protein_sequences = [
                d.uniprot_domain_1.uniprot_sequence.uniprot_sequence,
                d.uniprot_domain_2.uniprot_sequence.uniprot_sequence
            ]
            pdb_id = d.template.domain_1.pdb_id
            pdb_chains = [
                d.template.domain_1.pdb_chain, 
                d.template.domain_2.pdb_chain
            ]
            pdb_domain_defs = [
                d.template.domain_1.pdb_domain_def,
                d.template.domain_2.pdb_domain_def
            ]

        self.sequence_file, self.sequence_seqrecords = (
            self._write_domain_sequence_file(protein_ids, protein_domain_defs, protein_sequences)
        )
        
        self.structure_file, self.structure_seqrecords = (
            self._write_domain_structure_file(pdb_id, pdb_chains, pdb_domain_defs)
        )
        
    
    def run(self):
        self.model = model.Model(
            self.sequence_file, self.structure_file, self.modeller_results_file
        )
        
        
    def __exit__(self, exc_type, exc_value, traceback):
        d = self.d
        if exc_type is not None and exc_type in self.handled_errors:
            # Find domains that were used as a template and eventually led to
            # the error in the model, and add the error to their `domain_errors`
            # or `domain_contact_errors` fields.
            logger.error(exc_value)
            if isinstance(d, database.UniprotDomain):
                if d.template.model == None:
                    d.template.model = database_tables.UniprotDomainModel()
                    d.template.model.uniprot_domain_id = d.uniprot_domain_id
                bad_domains = self.db.get_rows_by_ids(
                    database_tables.Domain, [database_tables.Domain.cath_id], [d.template.cath_id])
                assert len(bad_domains) == 1
                bad_domain = bad_domains[0]
                bad_domain.domain_errors = str(d.uniprot_domain_id) + ': ' + str(exc_type)
                logger.error(
                    "Making a homology model failed!!!\n"
                    "Adding error '{0}' to the domain with cath_id {1}..."
                    .format(bad_domain.domain_errors, d.template.cath_id))
            elif isinstance(d, database_tables.UniprotDomainPair):
                if d.template.model == None:
                    d.template.model = database_tables.UniprotDomainPairModel()
                    d.template.model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
                bad_domains = self.db.get_rows_by_ids(
                    database_tables.DomainContact,
                    [database_tables.DomainContact.cath_id_1, 
                     database_tables.DomainContact.cath_id_2],
                    [d.template.cath_id_1, d.template.cath_id_2])
                if len(bad_domains) == 0:
                    bad_domains = self.db.get_rows_by_ids(
                        database_tables.DomainContact,
                        [database_tables.DomainContact.cath_id_1, 
                         database_tables.DomainContact.cath_id_2],
                        [d.template.cath_id_2, d.template.cath_id_1])
                assert len(bad_domains) == 1
                bad_domain = bad_domains[0]
                bad_domain.domain_contact_errors = (
                    str(d.uniprot_domain_pair_id) + ': ' + str(exc_type)
                )
                logger.error(
                    "Making a homology model failed!!!\n"
                    "Adding error '{0}' to the domain pair with cath_id_1 {1} "
                    "and cath_id_2 {2}..."
                    .format(bad_domain.domain_contact_errors, d.template.cath_id_1, 
                            d.template.cath_id_2))
                    
            # Add the error type to the model_errors column
            d.template.model.model_errors = (
                self.__add_new_error(d.template.model.model_errors, exc_value)
            )
            logger.error(d.template.model.model_errors)
            self.db.merge_row(bad_domain)
            # d.template.model = empty_model
            return True
        elif exc_type is not None:
            # Raise any exceptions that were not handled above...
            return False


        ### Upload model data to the database
        if self.modeller_results_file:
            # Model has already been calculated
            # TODO: there must be a more elegant solution...
            return
            
        d = self.d
        model = self.model
        
        ### Domains
        if isinstance(d, database.UniprotDomain):
            if d.template.model == None:
                d.template.model = database_tables.UniprotDomainModel()
                d.template.model.uniprot_domain_id = d.uniprot_domain_id
            
            d.template.model.model_filename = op.basename(model.modeller_results['model_file'])
            d.template.model.norm_dope = model.modeller_results['norm_dope']
            d.template.model.chain = model.chain_ids[0]
            d.template.model.alignment_filename = (
                op.basename(self.model.modeller_results['alignment_files'][0])
            )
            
            sasa_score = model.relative_sasa_scores[model.chain_ids[0]]
            d.template.model.sasa_score = ','.join('{:.2f}'.format(x) for x in sasa_score)
        
            d.template.model.model_domain_def = (
                self._truncate_domain_defs(
                    d.template.domain_def, 
                    model.modeller_results['domain_def_offsets'][0])
            )
        
        ### Domain pairs
        elif isinstance(d, database.UniprotDomainPair):
            if d.template.model == None:
                d.template.model = database.UniprotDomainPairModel()
                d.template.model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
            d.template.model.model_filename = op.basename(model.modeller_results['model_file'])
            d.template.model.norm_dope = model.modeller_results['norm_dope']
            d.template.model.chain_1 = model.chain_ids[0]
            d.template.model.chain_2 = model.chain_ids[1]
            d.template.model.alignment_filename_1 = (
                op.basename(self.model.modeller_results['alignment_files'][0])
            )
            d.template.model.alignment_filename_2 = (
                op.basename(self.model.modeller_results['alignment_files'][1])
            )
    
            model_domain_def_1 = self._truncate_domain_defs(
                d.uniprot_domain_1.template.domain_def,
                model.modeller_results['domain_def_offsets'][0])
    
            model_domain_def_2 = self._truncate_domain_defs(
                d.uniprot_domain_2.template.domain_def,
                model.modeller_results['domain_def_offsets'][1])
    
            # Convert interacting AA from indexes to uniprot numbers
            domain_start_1 = int(model_domain_def_1.split(':')[0])
            domain_start_2 = int(model_domain_def_2.split(':')[0])
            chain_1_interacting_uninum = [i + domain_start_1 for i in model.interacting_residues_1]
            chain_2_interacting_uninum = [i + domain_start_2 for i in model.interacting_residues_2]
    
            d.template.model.interacting_aa_1 = (
                ','.join([str(uniprot_num) for uniprot_num in chain_1_interacting_uninum])
            )
            d.template.model.interacting_aa_2 = (
                ','.join([str(uniprot_num) for uniprot_num in chain_2_interacting_uninum])
            )
            
            # Get interacting amino acids and interface area
            d.template.model.interface_area_hydrophobic = model.interface_area_hydrophobic
            d.template.model.interface_area_hydrophilic = model.interface_area_hydrophilic
            d.template.model.interface_area_total = model.interface_area_total
            
            # Save model_domain_defs, which might be truncated compared to uniprot_domain_template 
            # domain defs
            d.template.model.model_domain_def_1 = model_domain_def_1
            d.template.model.model_domain_def_2 = model_domain_def_2
    
    
        ### Values common for single domains and interactions
        model_errors =', '.join('{}'.format(e) for e in model.errors)
        if model_errors != '':
            d.template.model.model_errors = model_errors


        # TODO: Copy alignment and model file to the archive folder

        # Add either the empty model or the calculated model to the database
        logger.info('Adding model...')
        
        self.db.merge_model(d, model.modeller_results, configs['copy_data'] and d.path_to_data)
        
        
    def _truncate_domain_defs(self, domain_def, domain_def_offset):
        if (domain_def_offset is None or 
            domain_def_offset[0] is None or 
            domain_def_offset[1] is None):
            return domain_def
        n_gaps_start, n_gaps_end = domain_def_offset
        domain_def_new = (
            [str(int(domain_def.split(':')[0]) + n_gaps_start)] +
            domain_def.split(':')[1:-1] +
            [str(int(domain_def.split(':')[-1]) - n_gaps_end)])
        domain_def_new = ':'.join(domain_def_new)
        return domain_def_new


    def _write_domain_sequence_file(self, protein_ids, protein_domain_defs, protein_sequences):
        """Write a fasta file containing target domain sequences.
        """
        sequence_seqrecords = []
        for protein_id, protein_domain_def, protein_sequence in zip(
                protein_ids, protein_domain_defs, protein_sequences):
            #
            sequence_id = '{}.{}'.format(protein_id, protein_domain_def)
            domain_def = helper.decode_domain_def(
                protein_domain_def, merge=True, return_string=False)
            seqrecord = SeqRecord(
                id=protein_id, 
                name=sequence_id,
                seq=Seq(protein_sequence[domain_def[0]-1:domain_def[1]])
            )
            sequence_seqrecords.append(seqrecord)
        
        sequence_filename = '_'.join(seqrec.id for seqrec in sequence_seqrecords)
        sequence_file = op.join(configs['unique_temp_dir'], sequence_filename + '.fasta')
        with open(sequence_file, 'w') as ofh:
            SeqIO.write(sequence_seqrecords, ofh, 'fasta')
        
        return sequence_file, sequence_seqrecords
    
    
    def _write_domain_structure_file(self, pdb_id, pdb_chains, pdb_domain_defs):
        """Write a pdb file containing template domain chains (cut to domain bounaries).
        """
        pdb_file = structure_tools.get_pdb_file(pdb_id, configs['pdb_dir'], 'ent')
        sp = structure_tools.StructureParser(pdb_file, pdb_chains, pdb_domain_defs)
        sp.extract()
        sp.save_structure(configs['unique_temp_dir'])
        sp.save_sequences(configs['unique_temp_dir'])
    
        structure_file = op.join(
            configs['unique_temp_dir'], sp.pdb_id + ''.join(sp.chain_ids) + '.pdb'
        )
        
        structure_seqrecords = []
        for pdb_chain in pdb_chains:
            chain_sequence = sp.chain_sequence_dict[pdb_chain]
            seqrecord = SeqRecord(id=pdb_id+pdb_chain, seq=Seq(chain_sequence))
            structure_seqrecords.append(seqrecord)
            
        return structure_file, structure_seqrecords


    def __add_new_error(self, d_error_log, e):
        if d_error_log is None:
            return str(type(e))
        else:
            return '{}; {}: {}'.format(d_error_log, type(e), str(e))

    @property
    def result(self):
        return self.model
        


#%%
@execute_and_remember
class PrepareMutation:
    
    handled_exceptions = (
        errors.PDBError,
        errors.FoldxError,
        errors.ResourceError,
        errors.FoldXAAMismatchError,
    )
    
    
    def __init__(self, d, mutation, uniprot_id, sequence, model, db):
        print_header(d)
        
        self.sequence = sequence
        self.model = model
        self.db = db

        self.d = d
        self.uniprot_id = uniprot_id
        self.mutation = mutation
        
        self.skip = False
        
        if not mutation:
            logger.debug('Not evaluating mutations because no mutations specified...')
            self.skip = True
        elif ((isinstance(d, database_tables.UniprotDomain)
                and not d.template.domain) or
              (isinstance(d, database_tables.UniprotDomainPair)
                and not (d.template.domain_1 and d.template.domain_2))):
            logger.debug('Skipping because no structural template is availible...')
            self.skip = True
        elif d.template.model == None or d.template.model.model_filename == None:
            logger.debug('d.template.model: {}'.format(d.template.model))
            logger.debug('d.template.model.model_filename: {}'.format(
                getattr(d.template.model, 'model_filename', 'does not exist!')))
            logger.debug('Skipping because no model...')
            self.skip = True
        elif d.template.model.model_errors != None:
            logger.debug('Skipping because the model has errors: {}!'
                .format(d.template.model.model_errors))
            self.skip = True

        mutation_prototype = re.compile("^[A-z][0-9]+[A-z]$")
        if not mutation_prototype.match(mutation) or int(mutation[1:-1]) == 0:
            logger.error('The mutation {} is not a supported format! Skiping!'.format(mutation))
            self.skip = True

        # Check to see if we have a precalculated mutation. Skip if all
        # parameters have been calculated; otherwise analyse the remaining
        # parameters. Create an empty mutation if the mutation has not
        # been precalculated.
        precalculated_mutation = self.db.get_uniprot_mutation(d, mutation, self.uniprot_id, True)
        logger.info('Have the following precalculated mutation: {}'.format(precalculated_mutation))
        if (precalculated_mutation and
            precalculated_mutation.provean_score and
            precalculated_mutation.stability_energy_wt and
            precalculated_mutation.ddg != None):
            #
            logger.info('Mutation has already been completely evaluated. Skipping...')
            self.mut = precalculated_mutation
            self.skip = True
        elif precalculated_mutation:
            self.mut = precalculated_mutation
        elif isinstance(d, database_tables.UniprotDomain):
            # Construct empty models that will be used if we have errors
            self.mut = database_tables.UniprotDomainMutation()
            self.mut.uniprot_domain_id = d.uniprot_domain_id
            self.mut.uniprot_id = self.uniprot_id
            self.mut.mutation = self.mutation
        elif isinstance(d, database_tables.UniprotDomainPair):
            self.mut = database_tables.UniprotDomainPairMutation()
            self.mut.uniprot_domain_pair_id = d.uniprot_domain_pair_id
            self.mut.uniprot_id = self.uniprot_id
            self.mut.mutation = self.mutation


    def __bool__(self):
        return not self.skip


    def __enter__(self):
        pass


    def run(self):
        d = self.d
        uniprot_mutation = self.mut
        mut_data = self.get_mutation_data()

        ### Sequence features
        results = self.sequence.mutate(self.mutation)
        uniprot_mutation.provean_score = results['provean_score']
        uniprot_mutation.matrix_score = results['matrix_score']        
        
        ### Structure features
        results = self.model.mutate(mut_data.mutation_idx, mut_data.mutation_domain)
        uniprot_mutation.uniprot_id = self.uniprot_id #mut_data.uniprot_id_1
        uniprot_mutation.mutation = self.mutation #mut_data.mutation
        uniprot_mutation.chain_modeller = results['chain_modeller']
        uniprot_mutation.mutation_modeller = results['mutation_modeller']
        uniprot_mutation.model_filename_wt = (
            '{}_{}/{}'.format(
                self.uniprot_id, self.mutation, 
                op.basename(results['model_filename_wt']))
        )
        uniprot_mutation.model_filename_mut = (
            '{}_{}/{}'.format(
                self.uniprot_id, self.mutation, 
                op.basename(results['model_filename_mut']))
        )
        # Copy models to a 'archive_temp_dir`, so that `database.merge_mutation` knows where to 
        # find them.
        archive_model_file_wt = (
            op.join(
                configs['archive_temp_dir'], 
                d.path_to_data, 
                uniprot_mutation.model_filename_wt)
        )
        archive_model_file_mut = (
            op.join(
                configs['archive_temp_dir'], 
                d.path_to_data, 
                uniprot_mutation.model_filename_mut)
        )
        os.makedirs(op.dirname(archive_model_file_wt), exist_ok=True)
        os.makedirs(op.dirname(archive_model_file_mut), exist_ok=True)
        shutil.copy(results['model_filename_wt'], archive_model_file_wt) 
        shutil.copy(results['model_filename_mut'], archive_model_file_mut)
        #
        uniprot_mutation.stability_energy_wt = results['stability_energy_wt']
        uniprot_mutation.stability_energy_mut = results['stability_energy_mut']
        #
        uniprot_mutation.physchem_wt = (
            '{},{},{},{}'.format(*results['physchem_wt'])
        )
        uniprot_mutation.physchem_wt_ownchain = (
            '{},{},{},{}'.format(*results['physchem_ownchain_wt'])
        )
        uniprot_mutation.physchem_mut = (
            '{},{},{},{}'.format(*results['physchem_mut'])
        )
        uniprot_mutation.physchem_mut_ownchain = (
            '{},{},{},{}'.format(*results['physchem_ownchain_mut'])
        )
        #
        uniprot_mutation.secondary_structure_wt = results['secondary_structure_wt']
        uniprot_mutation.solvent_accessibility_wt = results['solvent_accessibility_wt']
        uniprot_mutation.secondary_structure_mut = results['secondary_structure_mut']
        uniprot_mutation.solvent_accessibility_mut = results['solvent_accessibility_mut']
        #
        if isinstance(d, database_tables.UniprotDomainPair):
            uniprot_mutation.analyse_complex_energy_wt = results['analyse_complex_energy_wt']
            uniprot_mutation.analyse_complex_energy_mut = results['analyse_complex_energy_mut']
            uniprot_mutation.contact_distance_wt = results['contact_distance_wt']
            uniprot_mutation.contact_distance_mut = results['contact_distance_mut']

        ### Machine learning
        pred = predictor.Predictor()
        row_idx = 0
        df = self.get_mutation_features(d, uniprot_mutation, row_idx=row_idx)
        uniprot_mutation.ddg = pred.score(df, mut_data.core_or_interface)
        logger.debug('Predicted ddg: {}'.format(uniprot_mutation.ddg))
        df.loc[row_idx, 'ddg'] = uniprot_mutation.ddg
        df.to_json()
      

    def get_mutation_data(self):
        """
        """
        d = self.d
        uniprot_id_1 = self.uniprot_id
        mutation = self.mutation
        
        if isinstance(d, database.UniprotDomain):
            logger.debug("Analyzing core mutation for uniprot: %s" % uniprot_id_1)
            logger.debug('model_domain_def: {}'.format(d.template.model.model_domain_def))
            domain_start, domain_end = (
                helper.decode_domain_def(d.template.model.model_domain_def)
            )
            mutation_idx = 0
            core_or_interface = 'core'
            
            if int(mutation[1:-1]) < domain_start or int(mutation[1:-1]) > domain_end:
                raise errors.MutationOutsideDomainError(
                    'Mutation {} falls outside the domain with definitions {}.'
                    .format(mutation, d.template.model.model_domain_def)
                )
            
        elif isinstance(d, database.UniprotDomainPair):
            mutation_pos = int(mutation[1:-1])
            interacting_aa_1 = self._get_interacting_aa(d, 1)
            interacting_aa_2 = self._get_interacting_aa(d, 2)
            assert uniprot_id_1 in [d.uniprot_domain_1.uniprot_id, d.uniprot_domain_2.uniprot_id]
            
            if (uniprot_id_1 == d.uniprot_domain_1.uniprot_id and 
                mutation_pos in interacting_aa_1):
                logger.debug('Mutation is inside the first domain.')
                logger.debug('model_domain_def: {}'.format(d.template.model.model_domain_def_1))
                domain_start, domain_end = helper.decode_domain_def(
                    d.template.model.model_domain_def_1)
                mutation_idx = 0
                core_or_interface = 'interface'
                
                
            elif (uniprot_id_1 == d.uniprot_domain_2.uniprot_id and
                  mutation_pos in interacting_aa_2):
                logger.debug('Mutation is inside the second domain.')
                logger.debug('model_domain_def: {}'.format(d.template.model.model_domain_def_2))
                domain_start, domain_end = helper.decode_domain_def(
                    d.template.model.model_domain_def_2)
                mutation_idx = 1
                core_or_interface = 'interface'
                
            else:
                # Mutation is outside the interface
                logger.error('Uniprot ID: {}\tMutation: {}'.format(
                    uniprot_id_1, mutation))
                logger.error('Uniprot ID 1: {}\tInteracting AA 1: {}'.format(
                    d.uniprot_domain_1.uniprot_id, interacting_aa_1))
                logger.error('Uniprot ID 2: {}\tInteracting AA 2: {}'.format(
                    d.uniprot_domain_2.uniprot_id, interacting_aa_2))
                raise errors.MutationOutsideInterfaceError(
                    'Mutated residue {} not involved in the interaction!'
                    .format(mutation[:-1])
                )    
    
        position_domain = int(mutation[1:-1]) - domain_start + 1 
        mutation_domain = mutation[0] + str(position_domain) + mutation[-1]

        class MutationData:
            def __init__(self, **kwargs):
                self.__dict__ = kwargs
        
        mut_data = MutationData(
            mutation_idx = mutation_idx,
            mutation_domain = mutation_domain,
            core_or_interface = core_or_interface,
        )
        return mut_data
     

    def _get_interacting_aa(self, d, domain_1or2=1):
        if domain_1or2 == 1:
            interacting_aa = d.template.model.interacting_aa_1
        elif domain_1or2 == 2:
            interacting_aa = d.template.model.interacting_aa_2
        else:
            raise Exception("`domain_1or2` should be either '1' or '2'!")
        return [int(uniprot_num) for uniprot_num in interacting_aa.split(',') if uniprot_num]

        
    def get_mutation_features(self, d, mut, row_idx=0):
        """
        Returns a dataframe that contains all features for the given mutation that are relevant 
        for machine learning.
    
        Parameters
        ----------
        d : sqlalchemy orm object
           Contains domain information for the given mutation
        mut : sqlalchemy orm object
           Contains information about the mutation in question
        """
        feature_dict = {
            key: value for (key, value) in mut.__dict__.items() if not key.startswith('_')
        }
    
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
    
        logger.debug('feature_dict: {}'.format(feature_dict))
        logger.debug('index: {}'.format(row_idx))
        feature_df = pd.DataFrame(feature_dict, index=[row_idx])
    
        return feature_df
    
    
    def __exit__(self, exc_type, exc_value, traceback):
        d = self.d
        mutation = self.mutation
        uniprot_mutation = self.mut
        
        if exc_type is not None and (exc_type in 
              (errors.MutationOutsideDomainError, 
               errors.MutationOutsideInterfaceError)):
            logger.debug('{}: {}; OK'.format(exc_type, exc_value))
        elif exc_type is not None and exc_type in self.handled_exceptions:
            uniprot_mutation.mutation_errors = '{}: {}'.format(exc_type, exc_value)
            logger.debug(uniprot_mutation.mutation_errors)
        elif exc_type is not None:
            return False
        
        logger.info('Adding mutation {}'.format(mutation))
        #logger.debug("Mutation attributes: {}".format(dir(uniprot_mutation)))
        for attr in dir(uniprot_mutation):
            if attr.startswith('_'):
                continue
            attr_field = getattr(uniprot_mutation, attr)
            attr_type = type(attr_field)
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
        logger.debug('Mergin mutation data with the database...')
        self.db.merge_mutation(uniprot_mutation, configs['copy_data'] and d.path_to_data)
        logger.debug('Done merging mutation data!')
        return True
        
        
    @property
    def result(self):
        return self


#%%
def get_unique_id(d):
    if isinstance(d, database_tables.UniprotDomain):
        return ('uniprot_domain_id', d.uniprot_domain_id)
    else:
        return('uniprot_domain_pair_id', d.uniprot_domain_pair_id)


def print_header(d):
    # logger.info('Domain or domain pair number: {}'.format(d_idx))
    logger.info('=' * 80)
    unique_id = get_unique_id(d)
    logger.info('{}: {}'.format(*unique_id))


def add_new_error(d_error_log, e):
    if d_error_log is None:
        return str(type(e))
    else:
        return '{}; {}: {}'.format(d_error_log, type(e), str(e))


