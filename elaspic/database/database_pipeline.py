import json
import logging
import os
import os.path as op
import re
import shutil

import numpy as np
import pandas as pd
import six
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from kmtools import structure_tools

import elaspic
import elaspic.database

logger = logging.getLogger(__name__)


class DatabasePipeline(elaspic.pipeline._Pipeline):

    def __init__(
            self, uniprot_id, mutations, configurations=None,
            run_type='5', number_of_tries=[], uniprot_domain_pair_ids=[]):
        """Run the main function of the program and parse elaspic.exc.

        Parameters
        ----------
        uniprot_domain_pair_ids : list of integers
            List of uniprot_domain_pair_ids specifying which uniprot domain pairs to analyse.
        run_type : str
            1 / sequence: Calculate provean supporting set.
            2 / model: Calculate all homology models.
            3 / mutations: Calculate mutations using precalculated provean and models
                (not working)
            4 : Calculate mutations using precalculated provean (calculating mutations on the fly)
                (not working)
            5 : Calculate mutations, calculating provean and homology models as required.
        """
        super().__init__(configurations)

        if isinstance(run_type, (int, float)):
            run_type = str(int(run_type))
        self.uniprot_id = uniprot_id
        self.mutations = self._split_mutations(mutations)
        self.calculated_mutations = []
        self.run_type = self._validate_run_type(run_type)
        assert not set(self.run_type.split('.')) - {'sequence', 'model', 'mutation'}

        self.number_of_tries = number_of_tries
        self.uniprot_domain_pair_ids = uniprot_domain_pair_ids

        logger.info('=' * 80)
        logger.info('## Input parameters')
        logger.info('uniprot_id: {}'.format(uniprot_id))
        logger.info('mutations: {}'.format(mutations))
        logger.info('run_type: {}'.format(run_type))
        logger.info('uniprot_domain_pair_ids: {}'.format(uniprot_domain_pair_ids))
        logger.info('unique_temp_dir: {}'.format(elaspic.CONFIGS['unique_temp_dir']))
        logger.info('db_schema: {}'.format(elaspic.CONFIGS.get('db_schema')))
        logger.info('temp_dir: {temp_dir}'.format(**elaspic.CONFIGS))
        logger.info('=' * 80)
        logger.info('')

        # Switch to the root of the unique tmp directory
        os.chdir(elaspic.CONFIGS['unique_temp_dir'])

        # Initialise the sql database for accessing all information
        self.db = elaspic.database.Database()

        # Obtain all domains and interactions for a given uniprot
        logger.info('Obtaining protein domain information...')
        self.uniprot_domains = self.db.get_uniprot_domain(self.uniprot_id, True)
        self._update_path_to_data(self.uniprot_domains)
        logger.info("Found {} domains.".format(len(self.uniprot_domains)))

        # Mutations results
        self.uniprot_mutations = []

    def _update_path_to_data(self, d_list):
        if not isinstance(d_list, list):
            d_list = [d_list]
        for d in d_list:
            if not d.path_to_data or any([len(x) > 255 for x in d.path_to_data.split('/')]):
                logger.debug("Updating 'path_to_data' from '{}'...".format(d.path_to_data))
                d.path_to_data = (
                    self.db.get_uniprot_base_path(d) +
                    self.db.get_uniprot_domain_path(d)
                )
                logger.debug("to '{}'...".format(d.path_to_data))
                self.db.merge_row(d)
            os.makedirs(
                op.join(elaspic.CONFIGS['archive_temp_dir'], d.path_to_data), exist_ok=True)

    def run(self):
        if not self.uniprot_domains:
            logger.info('Warning! Uniprot {} has no pfam domains!'.format(self.uniprot_id))
            return

        # Find provean
        if ((self.run_type in ['1', '5', '6'] or 'sequence' in self.run_type) and
                self.uniprot_domains):
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Computing provean...")
            self.get_sequence(self.uniprot_domains[0])

        # Get interactions
        if elaspic.CONFIGS['look_for_interactions']:
            logger.info('\n\n\n' + '*' * 110)
            logger.info('Obtaining protein domain pair information...')
            self.uniprot_domain_pairs = self.db.get_uniprot_domain_pair(
                self.uniprot_id, True, self.uniprot_domain_pair_ids)
            self._update_path_to_data(self.uniprot_domain_pairs)
            logger.info("Found {} interfaces.".format(len(self.uniprot_domain_pairs)))

        # Make models
        if self.run_type in ['2', '4', '5', '6'] or 'model' in self.run_type:
            logger.info('\n\n\n' + '*' * 110)
            logger.info("Building models...")
            for d in self.uniprot_domains + self.uniprot_domain_pairs:
                self.get_model(d)
            logger.info(
                'Finished processing all models for {} {}'
                .format(self.uniprot_id, self.mutations)
            )

        # Analyse mutations
        if self.run_type in ['3', '4', '5'] or 'mutation' in self.run_type:
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
        """Make a homology model for each uniprot domain that has a template in pdbfam."""
        logger.debug('-' * 80)
        logger.debug('get_model({})'.format(d))
        return PrepareModel(d, self.db)

    def get_mutation_score(self, d, mutation):
        logger.debug('-' * 80)
        logger.debug('get_mutation_score({}, {})'.format(d, mutation))
        if isinstance(d, elaspic.database.UniprotDomain):
            sequence = self.get_sequence(d)
        elif isinstance(d, elaspic.database.UniprotDomainPair):
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


class _PrepareSequence:

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
                elaspic.CONFIGS['archive_temp_dir'],
                self.db.get_uniprot_base_path(d),
                d.uniprot_sequence.provean.provean_supset_filename
            )
        else:
            self.provean_supset_file = None

    def __bool__(self):
        return not self.skip

    def __enter__(self):
        d = self.d
        self.sequence_file = op.join(elaspic.CONFIGS['sequence_dir'], d.uniprot_id + '.fasta')
        self.seqrecord = SeqRecord(
            id=d.uniprot_id,
            seq=Seq(d.uniprot_sequence.uniprot_sequence)
        )
        with open(self.sequence_file, 'w') as ofh:
            SeqIO.write(self.seqrecord, ofh, 'fasta')
        assert op.isfile(self.sequence_file)

    def run(self):
        self.sequence = elaspic.Sequence(self.sequence_file, self.provean_supset_file)
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

        if self.provean_supset_file is None:
            provean = elaspic.database.Provean()
            provean.uniprot_id = d.uniprot_id
            provean.provean_supset_filename = op.basename(self.sequence.provean_supset_file)
            provean.provean_supset_length = self.sequence.provean_supset_length
            self.db.merge_provean(
                provean,
                self.sequence.provean_supset_file,
                elaspic.CONFIGS['copy_data'] and self.db.get_uniprot_base_path(d)
            )

    @property
    def result(self):
        return self.sequence


PrepareSequence = DatabasePipeline.execute_and_remember(_PrepareSequence)


class _PrepareModel:

    handled_errors = (
        elaspic.exc.ModellerError,
        elaspic.exc.ChainsNotInteractingError,
        elaspic.exc.MutationOutsideDomainError,
        elaspic.exc.MutationOutsideInterfaceError,
    )
    bad_errors = (
        structure_tools.exc.PDBChainError,
        structure_tools.exc.PDBEmptySequenceError,
        structure_tools.exc.PDBNotFoundError,
        elaspic.exc.NoSequenceFound,
        elaspic.exc.TcoffeeError,
        elaspic.exc.InterfaceMismatchError,
    )

    def __init__(self, d, db, new_model=False):
        """
        """
        print_header(d)
        self.d = d
        self.db = db
        self.skip = False
        self.model = None
        self.modeller_results_file = None

        # Check if we should skip the model
        if new_model:
            pass
        elif (isinstance(d, elaspic.database.UniprotDomain) and
                not (d.template and d.template.cath_id)):
            logger.error('No structural template availible for this domain. Skipping...')
            self.skip = True
        elif (isinstance(d, elaspic.database.UniprotDomainPair) and not
                (d.template and d.template.cath_id_1 and d.template.cath_id_2)):
            logger.error('No structural template availible for this domain pair. Skipping...')
            self.skip = True
        elif d.template.model and d.template.model.model_filename:
            logger.info('Model already calculated.')
            self.modeller_results_file = self._get_modeller_results_file(d)
        elif (d.template.model and d.template.model.model_errors and
                (('Giving up' in d.template.model.model_errors) or
                 (d.template.model.model_elaspic.exc.count(';') >= 2))):
            logger.info(
                'Previous model had unfixable errors: "{}". Skipping...'
                .format(d.template.model.model_errors))
            self.skip = True

    # def domain_info(self):
    #     # Load data from database ORM to model input
    #     if isinstance(self.d, elaspic.database.UniprotDomain):
    #         protein_ids = [self.d.uniprot_id]
    #         if self.d.template.model and self.d.template.model.model_domain_def:
    #             protein_domain_defs = [self.d.template.model.model_domain_def]
    #         else:
    #             protein_domain_defs = [self.d.template.domain_def]
    #         protein_sequences = [self.d.uniprot_sequence.uniprot_sequence]
    #         pdb_id = self.d.template.domain.pdb_id
    #         pdb_chains = [self.d.template.domain.pdb_chain]
    #         pdb_domain_defs = [self.d.template.domain.pdb_domain_def]
    #
    #     elif isinstance(self.d, elaspic.database.UniprotDomainPair):
    #         protein_ids = [
    #             self.d.uniprot_domain_1.uniprot_id,
    #             self.d.uniprot_domain_2.uniprot_id
    #         ]
    #         if (self.d.template.model and
    #                 self.d.template.model.model_domain_def_1 and
    #                 self.d.template.model.model_domain_def_2):
    #             protein_domain_defs = [
    #                 self.d.template.model.model_domain_def_1,
    #                 self.d.template.model.model_domain_def_2
    #             ]
    #         else:
    #             protein_domain_defs = [
    #                 self.d.uniprot_domain_1.template.domain_def,
    #                 self.d.uniprot_domain_2.template.domain_def
    #             ]
    #         protein_sequences = [
    #             self.d.uniprot_domain_1.uniprot_sequence.uniprot_sequence,
    #             self.d.uniprot_domain_2.uniprot_sequence.uniprot_sequence
    #         ]
    #         pdb_id = self.d.template.domain_1.pdb_id
    #         pdb_chains = [
    #             self.d.template.domain_1.pdb_chain,
    #             self.d.template.domain_2.pdb_chain
    #         ]
    #         pdb_domain_defs = [
    #             self.d.template.domain_1.pdb_domain_def,
    #             self.d.template.domain_2.pdb_domain_def
    #         ]
    #
    #     return protein_domain_defs, protein_sequences, pdb_id, pdb_chains, pdb_domain_defs

    def _get_modeller_results_file(self, d):
        """Save relevant precalculated data to a json file."""
        modeller_results = dict()
        if isinstance(d, elaspic.database.UniprotDomain):
            unique_id = d.uniprot_domain_id
            modeller_results['model_file'] = op.join(
                elaspic.CONFIGS['archive_temp_dir'],
                d.path_to_data,
                d.template.model.model_filename
            )
            modeller_results['alignment_files'] = [
                op.join(
                    elaspic.CONFIGS['archive_temp_dir'],
                    d.path_to_data,
                    d.template.model.alignment_filename,
                )
            ]
            modeller_results['norm_dope'] = d.template.model.norm_dope
            modeller_results['domain_def_offsets'] = [(0, 0)]
            # modeller_results['model_domain_defs'] = [
            #     [int(i) for i in d.template.model.model_domain_def.split(':')],
            # ]
        elif isinstance(d, elaspic.database.UniprotDomainPair):
            unique_id = d.uniprot_domain_pair_id
            modeller_results['model_file'] = op.join(
                elaspic.CONFIGS['archive_temp_dir'],
                d.path_to_data,
                d.template.model.model_filename
            )
            modeller_results['alignment_files'] = [
                op.join(
                    elaspic.CONFIGS['archive_temp_dir'],
                    d.path_to_data,
                    d.template.model.alignment_filename_1,
                ),
                op.join(
                    elaspic.CONFIGS['archive_temp_dir'],
                    d.path_to_data,
                    d.template.model.alignment_filename_2,
                ),
            ]
            modeller_results['norm_dope'] = d.template.model.norm_dope
            modeller_results['domain_def_offsets'] = [(0, 0), (0, 0)]

            # modeller_results['model_domain_defs'] = [
            #     [int(i) for i in d.template.model.model_domain_def_1.split(':')],
            #     [int(i) for i in d.template.model.model_domain_def_2.split(':')],
            # ]

        modeller_results['model_domain_defs'] = [
            (0, int(x[1]) - int(x[0]) + 1) for x in self.d.protein_domain_defs
        ]
        modeller_results_file = op.join(
            elaspic.CONFIGS['model_dir'],
            '{}_modeller_results.json'.format(unique_id)
        )
        with open(modeller_results_file, 'w') as ofh:
            json.dump(modeller_results, ofh)
        return modeller_results_file

    def __bool__(self):
        return not self.skip

    def __enter__(self):
        self.sequence_file, self.sequence_seqrecords = self._write_domain_sequence_file()
        self.structure_file, self.structure_seqrecords = self._write_domain_structure_file()

    def run(self):
        self.model = elaspic.Model(
            self.sequence_file, self.structure_file, self.modeller_results_file
        )

    def __exit__(self, exc_type, exc_value, traceback):
        d = self.d
        if (exc_type is not None and
                exc_type in self.bad_errors and
                self.modeller_results_file is not None):
            logger.error('{}: {}'.format(exc_type, exc_value))
            logger.error("Removing model due to an error and trying again...")
            self.db.remove_model(self.d)
            self.d.template.model.model_filename = None
            self.model = PrepareModel(self.d, self.db, new_model=True)
            return True
        elif (exc_type is not None and exc_type in self.handled_errors):
            # Find domains that were used as a template and eventually led to
            # the error in the model, and add the error to their `domain_errors`
            # or `domain_contact_errors` fields.
            logger.error(exc_value)
            if isinstance(d, elaspic.database.UniprotDomain):
                if d.template.model == None:  # noqa
                    d.template.model = elaspic.database.UniprotDomainModel()
                    d.template.model.uniprot_domain_id = d.uniprot_domain_id
                bad_domains = self.db.get_rows_by_ids(
                    elaspic.database.Domain,
                    [elaspic.database.Domain.cath_id],
                    [d.template.cath_id])
                assert len(bad_domains) == 1
                bad_domain = bad_domains[0]
                bad_domain.domain_errors = str(d.uniprot_domain_id) + ': ' + str(exc_type)
                logger.error(
                    "Making a homology model failed!!!\n"
                    "Adding error '{0}' to the domain with cath_id {1}..."
                    .format(bad_domain.domain_errors, d.template.cath_id))
            else:  # UniprotDomainPair
                if d.template.model == None:  # noqa
                    d.template.model = elaspic.database.UniprotDomainPairModel()
                    d.template.model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
                bad_domains = self.db.get_rows_by_ids(
                    elaspic.database.DomainContact,
                    [elaspic.database.DomainContact.cath_id_1,
                     elaspic.database.DomainContact.cath_id_2],
                    [d.template.cath_id_1, d.template.cath_id_2])
                if len(bad_domains) == 0:
                    bad_domains = self.db.get_rows_by_ids(
                        elaspic.database.DomainContact,
                        [elaspic.database.DomainContact.cath_id_1,
                         elaspic.database.DomainContact.cath_id_2],
                        [d.template.cath_id_2, d.template.cath_id_1])
                if len(bad_domains) == 1:
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
                else:
                    logger.error("Found 0 or too many templates in the `domain_contact` table!")
                    logger.error("bad_domains: {}".format(bad_domains))
                    return True

            # Add the error type to the model_errors column
            d.template.model.model_errors = (
                self.__add_new_error(d.template.model.model_errors, exc_value)
            )
            logger.error(d.template.model.model_errors)
            logger.info('Adding bad domain and model...')
            self.db.merge_row(bad_domain)
            self.db.merge_model(d)
            # d.template.model = empty_model
            return True
        elif exc_type is not None:
            # Raise any exceptions that were not handled above...
            return False

        # Upload model data to the database
        if self.modeller_results_file:
            # Model has already been calculated
            # TODO: there must be a more elegant solution...
            return

        # Domains
        if isinstance(d, elaspic.database.UniprotDomain):
            if d.template.model == None:  # noqa
                d.template.model = elaspic.database.UniprotDomainModel()
                d.template.model.uniprot_domain_id = d.uniprot_domain_id

            d.template.model.model_filename = op.basename(
                self.model.modeller_results['model_file'])
            d.template.model.norm_dope = self.model.modeller_results['norm_dope']
            d.template.model.chain = self.model.chain_ids[0]
            d.template.model.alignment_filename = (
                op.basename(self.model.modeller_results['alignment_files'][0])
            )

            sasa_score = self.model.relative_sasa_scores[self.model.chain_ids[0]]
            d.template.model.sasa_score = ','.join('{:.2f}'.format(x) for x in sasa_score)

            d.template.model.model_domain_def = (
                self._truncate_domain_defs(
                    d.template.domain_def,
                    self.model.modeller_results['domain_def_offsets'][0])
            )

        # Domain pairs
        elif isinstance(d, elaspic.database.UniprotDomainPair):
            if d.template.model == None:  # noqa
                d.template.model = elaspic.database.UniprotDomainPairModel()
                d.template.model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
            d.template.model.model_filename = op.basename(
                self.model.modeller_results['model_file'])
            d.template.model.norm_dope = self.model.modeller_results['norm_dope']
            d.template.model.chain_1 = self.model.chain_ids[0]
            d.template.model.chain_2 = self.model.chain_ids[1]
            d.template.model.alignment_filename_1 = (
                op.basename(self.model.modeller_results['alignment_files'][0])
            )
            d.template.model.alignment_filename_2 = (
                op.basename(self.model.modeller_results['alignment_files'][1])
            )

            model_domain_def_1 = self._truncate_domain_defs(
                d.uniprot_domain_1.template.domain_def,
                self.model.modeller_results['domain_def_offsets'][0])

            model_domain_def_2 = self._truncate_domain_defs(
                d.uniprot_domain_2.template.domain_def,
                self.model.modeller_results['domain_def_offsets'][1])

            # Convert interacting AA from indexes to uniprot numbers
            domain_start_1 = int(model_domain_def_1.split(':')[0])
            domain_start_2 = int(model_domain_def_2.split(':')[0])
            chain_1_interacting_uninum = [
                i + domain_start_1 - 1 for i in self.model.interacting_aa_1]
            chain_2_interacting_uninum = [
                i + domain_start_2 - 1 for i in self.model.interacting_aa_2]

            d.template.model.interacting_aa_1 = (
                ','.join([str(uniprot_num) for uniprot_num in chain_1_interacting_uninum])
            )
            d.template.model.interacting_aa_2 = (
                ','.join([str(uniprot_num) for uniprot_num in chain_2_interacting_uninum])
            )

            # Get interacting amino acids and interface area
            d.template.model.interface_area_hydrophobic = self.model.interface_area_hydrophobic
            d.template.model.interface_area_hydrophilic = self.model.interface_area_hydrophilic
            d.template.model.interface_area_total = self.model.interface_area_total

            # Save model_domain_defs, which might be truncated compared to uniprot_domain_template
            # domain defs
            d.template.model.model_domain_def_1 = model_domain_def_1
            d.template.model.model_domain_def_2 = model_domain_def_2

        # Values common for single domains and interactions
        model_errors = ', '.join('{}'.format(e) for e in self.model.errors)
        if model_errors != '':
            d.template.model.model_errors = model_errors

        logger.info('Adding model...')
        self.db.merge_model(d, self.model.modeller_results)

    def _truncate_domain_defs(self, domain_def, domain_def_offset):
        n_gaps_start, n_gaps_end = domain_def_offset
        if n_gaps_start is None:
            n_gaps_start = 0
        if n_gaps_end is None:
            n_gaps_end = 0
        domain_def_new = (
            str(int(domain_def.split(':')[0]) + n_gaps_start) +
            ':' +
            str(int(domain_def.split(':')[-1]) - n_gaps_end)
        )
        return domain_def_new

    def _write_domain_sequence_file(self):
        """Write a fasta file containing target domain sequences."""
        sequence_seqrecords = []
        for protein_id, protein_domain_def, protein_sequence in zip(
                self.d.protein_ids, self.d.protein_domain_defs, self.d.protein_sequences):
            #
            sequence_id = '{}.{}'.format(protein_id, protein_domain_def)
            seqrecord = SeqRecord(
                id=sequence_id,
                name=protein_id,
                seq=Seq(protein_sequence[protein_domain_def])
            )
            sequence_seqrecords.append(seqrecord)

        sequence_filename = '_'.join(seqrec.id for seqrec in sequence_seqrecords)
        sequence_file = op.join(elaspic.CONFIGS['unique_temp_dir'], sequence_filename + '.fasta')
        with open(sequence_file, 'w') as ofh:
            SeqIO.write(sequence_seqrecords, ofh, 'fasta')

        return sequence_file, sequence_seqrecords

    def _write_domain_structure_file(self):
        """Write a pdb file containing template domain chains (cut to domain bounaries).

        TODO: Switch to using ``biounit=True``.
        """
        pdb_id = self.d.pdb_id
        pdb_chains = self.d.pdb_chains
        pdb_domain_defs = self.d.pdb_domain_defs
        model_id = 0

        s = structure_tools.fetch_structure(
            pdb_id, pdb_type='cif', biounit=False, pdb_mirror=elaspic.CONFIGS['pdb_dir'])
        ps = structure_tools.process_structure(s)
        pdb_domain_defs = _format_domain_defs_legacy(ps, pdb_chains, pdb_domain_defs)
        es = structure_tools.extract(ps, chain_ids=pdb_chains, domain_defs=pdb_domain_defs)

        structure_filename = op.join(
            elaspic.CONFIGS['unique_temp_dir'], pdb_id + ''.join(pdb_chains) + '.pdb')
        structure_tools.save_structure(es, structure_filename)

        structure_seqrecords = []
        for chain_id in pdb_chains:
            chain_sequence = structure_tools.get_chain_sequence(es[model_id][chain_id])
            seqrecord = SeqRecord(id=self.d.pdb_id + self.d.pdb_chains, seq=Seq(chain_sequence))
            structure_seqrecords.append(seqrecord)

        return structure_filename, structure_seqrecords

    def __add_new_error(self, d_error_log, e):
        if d_error_log is None:
            return str(type(e))
        else:
            return '{}; {}: {}'.format(d_error_log, type(e), str(e))

    @property
    def result(self):
        return self.model


PrepareModel = DatabasePipeline.execute_and_remember(_PrepareModel)


def _format_domain_defs_legacy(structure, pdb_chains, pdb_domain_defs):
    new_domain_defs = []
    logger.debug(sorted(structure.xtra['mapping'].keys()))
    logger.debug(sorted(structure.xtra['mapping'].values()))
    for chain_id, domain_def in zip(pdb_chains, pdb_domain_defs):
        _, _, new_domain_start = structure.xtra['mapping'][(0, chain_id, str(domain_def.start))]
        _, _, new_domain_stop = structure.xtra['mapping'][(0, chain_id, str(domain_def.stop))]
        new_domain_def = slice(new_domain_start, new_domain_stop)
        new_domain_defs.append(new_domain_def)
    return new_domain_defs


class _PrepareMutation:

    handled_exceptions = (
        structure_tools.exc.PDBError,
        elaspic.exc.FoldxError,
        elaspic.exc.ResourceError,
        elaspic.exc.FoldXAAMismatchError,
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
        elif ((isinstance(d, elaspic.database.UniprotDomain) and not
               d.template.domain) or
              (isinstance(d, elaspic.database.UniprotDomainPair) and not
               (d.template.domain_1 and d.template.domain_2))):
            logger.debug('Skipping because no structural template is availible...')
            self.skip = True
        elif (d.template.model == None or  # noqa
                d.template.model.model_filename == None):
            logger.debug('d.template.model: {}'.format(d.template.model))
            logger.debug('d.template.model.model_filename: {}'.format(
                getattr(d.template.model, 'model_filename', 'does not exist!')))
            logger.debug('Skipping because no model...')
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
                precalculated_mutation.ddg != None):  # noqa
            logger.info('Mutation has already been completely evaluated. Skipping...')
            self.mut = precalculated_mutation
            self.skip = True
        elif precalculated_mutation:
            self.mut = precalculated_mutation
        elif isinstance(d, elaspic.database.UniprotDomain):
            # Construct empty models that will be used if we have errors
            self.mut = elaspic.database.UniprotDomainMutation()
            self.mut.uniprot_domain_id = d.uniprot_domain_id
            self.mut.uniprot_id = self.uniprot_id
            self.mut.mutation = self.mutation
        elif isinstance(d, elaspic.database.UniprotDomainPair):
            self.mut = elaspic.database.UniprotDomainPairMutation()
            self.mut.uniprot_domain_pair_id = d.uniprot_domain_pair_id
            self.mut.uniprot_id = self.uniprot_id
            self.mut.mutation = self.mutation

    def __bool__(self):
        return not self.skip

    def __enter__(self):
        pass

    def run(self):
        d = self.d
        _domain_tables = (
            elaspic.database.UniprotDomain,
            elaspic.database.UniprotDomainPair
        )
        assert isinstance(d, _domain_tables)

        mut_data = self.get_mutation_data()

        # Sequence features
        results = self.sequence.mutate(self.mutation)
        self.mut.provean_score = results['provean_score']
        self.mut.matrix_score = results['matrix_score']

        # Structure features
        results = self.model.mutate(mut_data.mutation_idx, mut_data.mutation_domain)
        self.mut.uniprot_id = self.uniprot_id  # mut_data.uniprot_id_1
        self.mut.mutation = self.mutation  # mut_data.mutation
        self.mut.chain_modeller = results['chain_modeller']
        self.mut.mutation_modeller = results['mutation_modeller']
        self.mut.model_filename_wt = (
            '{}_{}/{}'.format(
                self.uniprot_id, self.mutation,
                op.basename(results['model_file_wt']))
        )
        self.mut.model_filename_mut = (
            '{}_{}/{}'.format(
                self.uniprot_id, self.mutation,
                op.basename(results['model_file_mut']))
        )
        # Copy models to a 'archive_temp_dir`, so that `database.merge_mutation` knows where to
        # find them.
        archive_model_file_wt = (
            op.join(
                elaspic.CONFIGS['archive_temp_dir'],
                d.path_to_data,
                self.mut.model_filename_wt)
        )
        archive_model_file_mut = (
            op.join(
                elaspic.CONFIGS['archive_temp_dir'],
                d.path_to_data,
                self.mut.model_filename_mut)
        )
        os.makedirs(op.dirname(archive_model_file_wt), exist_ok=True)
        os.makedirs(op.dirname(archive_model_file_mut), exist_ok=True)
        shutil.copy(results['model_file_wt'], archive_model_file_wt)
        shutil.copy(results['model_file_mut'], archive_model_file_mut)
        #
        self.mut.stability_energy_wt = results['stability_energy_wt']
        self.mut.stability_energy_mut = results['stability_energy_mut']
        #
        self.mut.physchem_wt = (
            '{},{},{},{}'.format(*results['physchem_wt'])
        )
        self.mut.physchem_wt_ownchain = (
            '{},{},{},{}'.format(*results['physchem_ownchain_wt'])
        )
        self.mut.physchem_mut = (
            '{},{},{},{}'.format(*results['physchem_mut'])
        )
        self.mut.physchem_mut_ownchain = (
            '{},{},{},{}'.format(*results['physchem_ownchain_mut'])
        )
        #
        self.mut.secondary_structure_wt = results['secondary_structure_wt']
        self.mut.solvent_accessibility_wt = results['solvent_accessibility_wt']
        self.mut.secondary_structure_mut = results['secondary_structure_mut']
        self.mut.solvent_accessibility_mut = results['solvent_accessibility_mut']
        #
        if isinstance(d, elaspic.database.UniprotDomainPair):
            self.mut.analyse_complex_energy_wt = results['analyse_complex_energy_wt']
            self.mut.analyse_complex_energy_mut = results['analyse_complex_energy_mut']
            self.mut.contact_distance_wt = results['contact_distance_wt']
            self.mut.contact_distance_mut = results['contact_distance_mut']

        # Machine learning
        if isinstance(d, elaspic.database.UniprotDomain):
            pred = elaspic.Predictor()
        else:
            pred = elaspic.Predictor()
        pred.load(elaspic.CACHE_DIR)
        row_idx = 0
        df = self.get_mutation_features(d, self.mut, row_idx=row_idx)
        self.mut.ddg = pred.score(df)[0]
        logger.debug('Predicted ddg: {}'.format(self.mut.ddg))
        df.loc[row_idx, 'ddg'] = self.mut.ddg
        # df.to_json()  # this seems to cause segfaults!?

    def get_mutation_data(self):
        """
        """
        d = self.d
        uniprot_id_1 = self.uniprot_id
        mutation = self.mutation

        if isinstance(d, elaspic.database.UniprotDomain):
            logger.debug("Analyzing core mutation for uniprot: %s" % uniprot_id_1)
            logger.debug('model_domain_def: {}'.format(d.template.model.model_domain_def))
            domain_start, domain_end = d.protein_domain_defs[0]
            mutation_idx = 0
            core_or_interface = 'core'
            domain_start = domain_start or 1
            domain_end = domain_end or len(self.sequence.sequence)

            if int(mutation[1:-1]) < domain_start or int(mutation[1:-1]) > domain_end:
                raise elaspic.exc.MutationOutsideDomainError(
                    'Mutation {} falls outside the domain with definitions {}.'
                    .format(mutation, d.template.model.model_domain_def)
                )

        elif isinstance(d, elaspic.database.UniprotDomainPair):
            mutation_pos = int(mutation[1:-1])
            interacting_aa_1 = self._get_interacting_aa(d, 1)
            interacting_aa_2 = self._get_interacting_aa(d, 2)
            assert uniprot_id_1 in [d.uniprot_domain_1.uniprot_id, d.uniprot_domain_2.uniprot_id]
            logger.debug('sequence: {}'.format(self.sequence.sequence))

            if (uniprot_id_1 == d.uniprot_domain_1.uniprot_id and
                    mutation_pos in interacting_aa_1):
                logger.debug('Mutation is inside the first domain.')
                logger.debug('model_domain_def: {}'.format(d.template.model.model_domain_def_1))
                domain_start, domain_end = d.protein_domain_defs[0]
                if domain_start is None:
                    domain_start = 1
                mutation_idx = 0
                core_or_interface = 'interface'
                domain_start = domain_start or 1  # if None
                domain_end = domain_end or len(self.sequence.sequence)

            elif (uniprot_id_1 == d.uniprot_domain_2.uniprot_id and
                    mutation_pos in interacting_aa_2):
                logger.debug('Mutation is inside the second domain.')
                logger.debug('model_domain_def: {}'.format(d.template.model.model_domain_def_2))
                domain_start, domain_end = d.protein_domain_defs[1]
                mutation_idx = 1
                core_or_interface = 'interface'
                domain_start = domain_start or 1
                domain_end = domain_end or len(self.sequence.sequence)

            else:
                # Mutation is outside the interface
                error_message = (
                    'Mutated residue {} not involved in the interaction!'
                    .format(mutation[:-1])
                )
                logger.error(error_message)
                logger.error('Uniprot ID: {}\tMutation: {}'.format(
                    uniprot_id_1, mutation))
                logger.error('Uniprot ID 1: {}\tInteracting AA 1: {}'.format(
                    d.uniprot_domain_1.uniprot_id, interacting_aa_1))
                logger.error('Uniprot ID 2: {}\tInteracting AA 2: {}'.format(
                    d.uniprot_domain_2.uniprot_id, interacting_aa_2))
                raise elaspic.exc.MutationOutsideInterfaceError(error_message)

        position_domain = int(mutation[1:-1]) - domain_start + 1
        mutation_domain = mutation[0] + str(position_domain) + mutation[-1]

        class MutationData:
            def __init__(self, **kwargs):
                self.__dict__ = kwargs

        mut_data = MutationData(
            mutation_idx=mutation_idx,
            mutation_domain=mutation_domain,
            core_or_interface=core_or_interface,
        )
        return mut_data

    def _get_interacting_aa(self, d, domain_1or2=1):
        if domain_1or2 == 1:
            interacting_aa = d.template.model.interacting_aa_1
        elif domain_1or2 == 2:
            interacting_aa = d.template.model.interacting_aa_2
        else:
            raise Exception("`domain_1or2` should be either '1' or '2'!")
        return (
            [int(uniprot_num) for uniprot_num in interacting_aa.split(',') if uniprot_num]
            if interacting_aa else []
        )

    def get_mutation_features(self, d, mut, row_idx=0):
        """Return a DataFrame containing features relevant for machine learning.

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
                'alignment_identity': np.sqrt(d.template.identical_1 * d.template.identical_2),
                'alignment_coverage': np.sqrt(d.template.coverage_1 * d.template.coverage_2),
                'alignment_score': np.sqrt(d.template.score_1 * d.template.score_2),
                #
                'interface_area_hydrophobic': d.template.model.interface_area_hydrophobic,
                'interface_area_hydrophilic': d.template.model.interface_area_hydrophilic,
                'interface_area_total': d.template.model.interface_area_total,
            })

        logger.debug('feature_dict: {}'.format(feature_dict))
        logger.debug('index: {}'.format(row_idx))
        feature_df = pd.DataFrame(feature_dict, index=[row_idx])

        return feature_df

    def __exit__(self, exc_type, exc_value, traceback):
        d = self.d
        mutation = self.mutation
        if (exc_type is not None and
                exc_type in (elaspic.exc.MutationOutsideDomainError,
                             elaspic.exc.MutationOutsideInterfaceError)):
            logger.debug('{}: {}; OK'.format(exc_type, exc_value))
            return True
        elif exc_type is not None and exc_type in self.handled_exceptions:
            self.mut.mutation_errors = '{}: {}'.format(exc_type, exc_value)
            logger.debug(self.mut.mutation_errors)
        elif exc_type is not None:
            logger.error('Unhandled error occured: {}'.format(exc_type))
            return False

        logger.info('Adding mutation {}'.format(mutation))
        # logger.debug("Mutation attributes: {}".format(dir(self.mut)))
        for attr in dir(self.mut):
            if attr.startswith('_'):
                continue
            attr_field = getattr(self.mut, attr)
            attr_type = type(attr_field)
            if isinstance(attr_field, six.binary_type):
                logger.debug(
                    'Changing attribute {} from {} to {}...'
                    .format(attr, attr_type, type(attr_field.decode())))
                setattr(self.mut, attr, attr_field.decode())
        logger.debug('Mergin mutation data with the database...')
        self.db.merge_mutation(self.mut, elaspic.CONFIGS['copy_data'] and d.path_to_data)
        logger.debug('Done merging mutation data!')

    @property
    def result(self):
        return self


PrepareMutation = DatabasePipeline.execute_and_remember(_PrepareMutation)


def get_unique_id(d):
    if isinstance(d, elaspic.database.UniprotDomain):
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
