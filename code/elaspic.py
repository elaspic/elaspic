# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 18:37:36 2012

@author: niklas
"""
import os
import logging
import subprocess
import tempfile
import argparse
import atexit
import signal

from modeller import ModellerError

import errors
import sql_db

import domain_template
import domain_model
import domain_mutation

from Bio.SubsMat import MatrixInfo
from ConfigParser import SafeConfigParser
import helper_functions as hf

blacklisted_uniprots = [
    'Q8WZ42', # Titin has a large number of domains that are joined in messed up ways
]


class ProteinDomainData(object):

    def __init__(self, domain, template, model, mutations=[]):
        self.d = domain
        self.t = template
        self.m = model
        self.mut = mutations


class Pipeline(object):

    def __init__(self, configFile):
        #######################################################################
        # Read the configuration file and set the variables
        configParser = SafeConfigParser(
            defaults={
                'global_temp_path': '/tmp/',
                'temp_path': '/tmp/elaspic/',
                'DEBUG': True,
                'path_to_archive': '/home/kimlab1/database_data/elaspic/',
                'look_for_interactions': 'True',
                'n_cores': 1,
                'db_type': sql_db.sql_flavor,
                'db_path': hf.path_to_pipeline_code() + '/../db/pipeline.db',
                'web_server': False,
            })
        configParser.read(configFile)

        # From [DEFAULT]
        self.global_temp_path = configParser.get('DEFAULT', 'global_temp_path')
        temp_path = configParser.get('DEFAULT', 'temp_path').strip('/') + '/' # So it never starts with '/';
        self.debug = configParser.getboolean('DEFAULT', 'DEBUG')
        self.path_to_archive = configParser.get('DEFAULT', 'path_to_archive')
        self.look_for_interactions = configParser.getboolean('DEFAULT', 'look_for_interactions')
        self.n_cores = configParser.get('DEFAULT', 'n_cores')
        self.db_type = configParser.get('DEFAULT', 'db_type')
        self.db_path = configParser.get('DEFAULT', 'db_path')
        self.web_server = configParser.get('DEFAULT', 'web_server')

        # From [SETTINGS]
        self.blast_db_path = configParser.get('SETTINGS', 'blast_db_path')
        self.pdb_path = configParser.get('SETTINGS', 'pdb_path')
        self.output_path = configParser.get('SETTINGS', 'output_path')
        self.bin_path = configParser.get('SETTINGS', 'bin_path')
        self.max_runtime = configParser.get('SETTINGS', 'max_runtime')

        # From [GET_TEMPLATE]

        # From [GET_MODEL]
        self.modeller_runs = configParser.getint('GET_MODEL', 'modeller_runs')

        # From [GET_MUTATION]
        self.foldx_water = configParser.get('GET_MUTATION', 'foldx_water')
        self.foldx_num_of_runs = configParser.get('GET_MUTATION', 'foldx_num_of_runs')
        self.matrix_type = configParser.get('GET_MUTATION', 'matrix_type')
        self.gap_start = configParser.getint('GET_MUTATION', 'gap_start')
        self.gap_extend = configParser.getint('GET_MUTATION', 'gap_extend')
        self.matrix = getattr(MatrixInfo, self.matrix_type)

        #######################################################################
        # If a TMPDIR is given as environment variable the tmp directory
        # is created relative to that. This is useful when running on banting
        # (the cluster in the ccbr) and also on Scinet (I might have set the
        # environment variable on Scinet myself though..). Make sure that it
        # points to '/dev/shm/' on Scinet
        tmpdir_cluster = hf.get_echo('$TMPDIR')
        if tmpdir_cluster:
            self.temp_path = os.path.join(tmpdir_cluster, temp_path)
        else:
            self.temp_path = os.path.join(self.global_temp_path, temp_path)
        subprocess.check_call('mkdir -p ' + self.temp_path, shell=True)

        #######################################################################
        # Copy the PDBFam database file
        if self.db_type == 'sqlite_file' and not os.path.isfile(self.temp_path + '/' + self.db_path.strip().split('/')[-1]):
            print self.temp_path
            print self.db_path
            system_command = 'cp -u ' + self.db_path + ' ' + self.temp_path + '/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, error_message = childProcess.communicate()
            if childProcess.returncode != 0:
                print result
                print error_message
                raise Exception('Could not copy the sqlite3 database!')

        #######################################################################
        # If running on the cluster copy the database to the tmp DIR for
        # speedup and to avoid killing the network. BLAST is very IO intensive
        # and you don't want that to be run over the network!
        # I can distinguish where the pipeline is running by checking the username
        # you will have to adjust that!
        # my usernames are:
        # local: niklas
        # banting: nberliner
        # Scinet: joan
        username  = hf.get_username()
        hostname = hf.get_hostname()
        if username == 'strokach' and (('beagle' in hostname) or ('grendel' in hostname)):
            # Copy the blast database if running on beagle or Grendel
            system_command = (
                'mkdir -p ' + self.global_temp_path + 'blast/pdbaa_db/ && ' +
                'mkdir -p ' + self.global_temp_path + 'blast/db/ && ' +
                'rsync -rzu ' + self.blast_db_path + 'pdbaa_db/ ' + self.global_temp_path + 'blast/pdbaa_db/ && ' +
                'rsync -rzu ' + self.blast_db_path + 'db/ ' + self.global_temp_path + 'blast/db/')
        elif username == 'strokach' or username == 'alexey':
            # Use a symbolic link to the blast database and use the home folder for temporary storage
            system_command = (
                'rm -rf ' + self.global_temp_path + 'blast/ && ' +
                'mkdir -p ' + self.global_temp_path + 'blast/ && ' +
                'ln -sf ' + self.blast_db_path + 'pdbaa_db ' + self.global_temp_path + 'blast/ && ' +
                'ln -sf ' + self.blast_db_path + 'db ' +  self.global_temp_path + 'blast/')
        elif username == 'witvliet':
            system_command = (
                'mkdir -p ' + self.temp_path + 'blast && ' +
                'cd ' + self.temp_path + 'blast && ' +
                'ln -sf /home/witvliet/working/bin/ncbi-blast-2.2.28+/pdbaa_db')
        elif username == 'joan':
            # For scinet, blast is already installed, but the database needs to be copied
            system_command = (
                'mkdir -p ' + self.global_temp_path + 'blast && ' +
                'cp -ru $HOME/niklas-pipeline/blastdb/pdbaa_db ' +
                self.global_temp_path + 'blast/')
        print system_command
        rc = 1
        n_tries = 0
        while rc != 0 and n_tries < 10:
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, error_message = childProcess.communicate()
            rc = childProcess.returncode
        if rc != 0:
            print result
            print error_message
            raise Exception(
                'Couldn\'t copy the blast database on {hostname} as {username}!'
                .format(hostname=hostname, username=username))



    def __call__(self, uniprot_id, mutations, run_type=1, n_cores=None, number_of_tries=[]):
        """ Run the main function of the program and parse errors
        """
        # Initialize the logger
        logger = logging.getLogger(__name__)
        if self.debug:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.handlers = []
        logger.addHandler(handler)
        self.log = logger

        #
        self.uniprot_id = uniprot_id
        self.mutations = mutations
        self.run_type = run_type
        if n_cores is not None:
            self.n_cores = n_cores
        self.number_of_tries = number_of_tries
        self.unique = tempfile.mkdtemp(prefix='', dir=self.temp_path).split('/')[-1]
        self.log.info(self.uniprot_id)
        self.log.info(self.unique)
        self.unique_temp_folder = self.temp_path + self.unique + '/'
        # Have to make a temporary folder for provean in my home directory because
        # there usually is not enough space on these nodes in the cluster
        username  = hf.get_username()
        hostname = hf.get_hostname()
        if username == 'strokach' and not (('beagle' in hostname) or ('grendel' in hostname)):
            self.log.debug('Using a temp folder on kimstg for provean temp files. This may lead to poor performace...')
            self.provean_temp_path = '/home/kimlab1/strokach/tmp/elaspic/' + self.unique
            subprocess.check_call('mkdir -p ' + self.provean_temp_path, shell=True)
        else:
            self.provean_temp_path = self.unique_temp_folder
        atexit.register(self.__clear_provean_temp_files)

        # Create temporary folders
        self.__prepare_temp_folder()

        # Initialise the sql database for accessing all information
        self.log.debug('Initializing the database...')
        self.db = sql_db.MyDatabase(
            path_to_sqlite_db=self.temp_path+'pipeline.db', sql_flavor=self.db_type,
            path_to_temp=self.temp_path, path_to_archive=self.path_to_archive,
            is_immutable=False, log=self.log)

        # Switch to the root of the unique tmp directory
        os.chdir(self.unique_temp_folder)

        # Initialize classes used for calculating templates, models and mutations
        self.get_template = domain_template.GetTemplate(
            self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
            self.db, self.log, self.n_cores,
            self.provean_temp_path)

        self.get_model = domain_model.GetModel(
            self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
            self.db, self.log, self.n_cores, self.modeller_runs)

        self.get_mutation = domain_mutation.GetMutation(
            self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
            self.db, self.log, self.n_cores, self.foldx_water,
            self.foldx_num_of_runs, self.matrix, self.gap_start, self.gap_extend)

        #######################################################################
        # Find all uniprot_domain and uniprot_domain_pairs for a given unirpot,
        # and add path_to_data variable to the domain if it does not exist
        self.log.debug('Obtaining protein domain information...')
        protein_domains = self.db.get_uniprot_domain(self.uniprot_id)
        for d, t, m in protein_domains:
            if not d.path_to_data:
                d.path_to_data = (
                    '{organism_name}/{uniprot_id_part1}/{uniprot_id_part2}/{uniprot_id_part3}/{pfam_name}*{envelope_def}/'
                    .format(
                        organism_name = d.uniprot_name.split('_')[-1].lower(),
                        uniprot_id_part1 = d.uniprot_id[:3],
                        uniprot_id_part2 = d.uniprot_id[3:5],
                        uniprot_id_part3 = d.uniprot_id,
                        pfam_name = d.pfam_name,
                        envelope_def = d.envelope_def.replace(':','-'),))
                self.db.add_domain(d)
            subprocess.check_call('mkdir -p ' + self.temp_path + d.path_to_data, shell=True)
        protein_domain_template_model = protein_domains
        self.log.debug('Obtaining protein domain pair information...')
        if self.look_for_interactions:
            protein_domain_pairs = self.db.get_uniprot_domain_pair(self.uniprot_id)
            for d, t, m in protein_domain_pairs:
                if not d.path_to_data:
                    d.path_to_data = (
                    '{organism_name}/{uniprot_id_1_part1}/{uniprot_id_1_part2}/{uniprot_id_1_part3}/'
                    '{pfam_name_1}*{envelope_def_1}/{pfam_name_2}*{envelope_def_2}/{uniprot_id_2}/'
                    .format(
                        organism_name = d.uniprot_domain_1.uniprot_name.split('_')[-1].lower(),
                        uniprot_id_1_part1 = d.uniprot_domain_1.uniprot_id[:3],
                        uniprot_id_1_part2 = d.uniprot_domain_1.uniprot_id[3:5],
                        uniprot_id_1_part3 = d.uniprot_domain_1.uniprot_id,
                        pfam_name_1 = d.uniprot_domain_1.pfam_name,
                        envelope_def_1 = d.uniprot_domain_1.envelope_def.replace(':','-'),
                        pfam_name_2 = d.uniprot_domain_2.pfam_name,
                        envelope_def_2 = d.uniprot_domain_2.envelope_def.replace(':','-'),
                        uniprot_id_2 = d.uniprot_domain_2.uniprot_id,))
                self.db.add_domain(d)
                subprocess.check_call('mkdir -p ' + self.temp_path + d.path_to_data, shell=True)
            # Protein domains have to come before protein domain pairs in order for Provean to work correctly
            protein_domain_template_model = protein_domains + protein_domain_pairs
        # Number of tries that we attempted to make a homology model
        # for each protein domain or domain contact
        if not self.number_of_tries:
            self.number_of_tries = [0] * len(protein_domain_template_model)
        if not protein_domain_template_model:
            self.log.error('Uniprot {} has no pfam domains'.format(self.uniprot_id))

        # Convert to a single object which holds domain, template, model, and mutations
        # (so that I can modify each of those parameters "in place" in a for loop)
        self.protein_definitions = []
        for d, t, m in protein_domain_template_model:
            self.protein_definitions.append(ProteinDomainData(d, t, m))

        # Do the calculations
        self.log.info('Choosing which type of analysis to perform...')
        if run_type in [1, 4, 5]:
            self.log.info("Finding templates...")
            self._compute_templates()
            self.log.info("Building models...")
            self._compute_models()
        if run_type in [3, 4, 5]:
            self.log.info("Analyzing mutations...")
            self._compute_mutations()
        if run_type in [2, 5]:
            self.log.info("Computing provean...")
            self._compute_provean()

        return self.protein_definitions


    def _compute_templates(self):
        """
        Align uniprot domains to structural domains in pdbfam, find the best
        template, and expand domain boundaries
        """
        possible_template_errors = (
            errors.NoStructuralTemplates,
            errors.NoSequenceFound,
            errors.TcoffeeError,
            errors.ProveanError,
            errors.NoTemplatesFound,)

        model_errors_to_redo = [
            "<class '_modeller.ModellerError'>",
            "<class '_modeller.SequenceMismatchError'>",]

        # Go over all domains and domain pairs for a given protein
        for p in self.protein_definitions:
            self.__print_header(p)

            if not p.t:
                # Construct empty templates that will be used if we have errors,
                # but only do this if we don't have a template already.
                # Otherwide we could replace existing templates with None if existing
                # templates have been blacklisted.
                if isinstance(p.d, sql_db.UniprotDomain):
                    p.t = sql_db.UniprotDomainTemplate()
                    p.t.uniprot_domain_id = p.d.uniprot_domain_id
                elif isinstance(p.d, sql_db.UniprotDomainPair):
                    p.t = sql_db.UniprotDomainPairTemplate()
                    p.t.uniprot_domain_pair_id = p.d.uniprot_domain_pair_id
            elif (p.m and p.m.model_errors and ('Giving up' not in p.m.model_errors) and
                    any([(error_string in p.m.model_errors) for error_string in model_errors_to_redo])):
                self.log.debug('Recalculating templates because the model had errors...')
            elif (isinstance(p.t, sql_db.UniprotDomainPairTemplate) and
                    p.t.domain_1 and
                    p.t.domain_2 and
                    (p.t.domain_1 == p.t.domain_2) and
                    (p.t.domain_1.pdb_chain == p.t.domain_2.pdb_chain)):
                self.log.debug('Recalculating templates because both interacting partners mapped to the same pdb chain...')
            else:
                self.log.info('Skipping template...')
                continue

            try: # Find templates and catch errors
                p.t = self.get_template(p.d)
            except errors.NoTemplatesFound as e:
                # This part is added so that you never end up without a template
                # if you had a template before. May not be necessary
                self.log.error('{}; Skipping because we previously had a template...'.format(type(e)))
                continue
            except possible_template_errors as e:
                p.t.template_errors = (
                    '{}'.format(type(e))
                    if (p.t.template_errors == None)
                    else '{}; {}'.format(p.t.template_errors, type(e)))
                self.log.error(p.t.template_errors)
            self.log.info('adding template\n\n\n')
            self.db.add_uniprot_template(p.t, p.d.path_to_data)
        self.log.info('Finished processing all templates for {} {}\n'.format(self.uniprot_id, self.mutations))


    def _compute_models(self):
        """
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        possible_model_errors = (
            ModellerError,
            errors.PDBChainError,
            errors.ModellerFailure,
            errors.ChainsNotInteracting,)

        # Go over all domains and domain pairs for a given protein
        for p_idx, p in enumerate(self.protein_definitions):
            self.__print_header(p)

            # Check if we should skip the model
            if p.t == None:
                self.log.error('No template information. Skipping...')
                continue
            elif (isinstance(p.d, sql_db.UniprotDomain) and not p.t.domain_def):
                self.log.error('Domain definitions missing. Skipping...')
                continue
            elif (isinstance(p.d, sql_db.UniprotDomainPair) and not
                    (p.t.domain_def_1 and p.t.domain_def_2)):
                self.log.error('Domain definitions missing. Skipping...')
                continue
            elif p.m and p.m.model_filename:
                self.log.info('Model already calculated. Skipping...')
                continue
            elif p.m and ('Giving up' in p.m.model_errors):
                self.log.info(
                    'Previous model had unfixable errors {}. Skipping...'
                    .format(p.m.model_errors))
                continue

            # Make a homology model using the domain and template information
            try:
                p.m = self.get_model(p.d, p.t)
            except possible_model_errors as e:
                # Find domains that were used as a template and eventually led to
                # the error in the model
                if isinstance(p.d, sql_db.UniprotDomain):
                    empty_model = sql_db.UniprotDomainModel()
                    empty_model.uniprot_domain_id = p.t.uniprot_domain_id
                    bad_domains = self.db._get_rows_by_ids(
                        sql_db.Domain, [sql_db.Domain.cath_id], [p.t.cath_id])
                    assert len(bad_domains) == 1
                    bad_domains[0].domain_errors = str(p.d.uniprot_domain_id) + ': ' + str(type(e))
                    error_handling_message = (
                        'Adding error to the domain with cath_id {cath_id} and running the pipeline again...'
                        .format(cath_id=p.t.cath_id))
                elif isinstance(p.d, sql_db.UniprotDomainPair):
                    empty_model = sql_db.UniprotDomainPairModel()
                    empty_model.uniprot_domain_pair_id = p.t.uniprot_domain_pair_id
                    bad_domains = self.db._get_rows_by_ids(
                        sql_db.DomainContact,
                        [sql_db.DomainContact.cath_id_1, sql_db.DomainContact.cath_id_2],
                        [p.t.cath_id_1, p.t.cath_id_2])
                    if len(bad_domains) == 0:
                        bad_domains = self.db._get_rows_by_ids(
                            sql_db.DomainContact,
                            [sql_db.DomainContact.cath_id_1, sql_db.DomainContact.cath_id_2],
                            [p.t.cath_id_2, p.t.cath_id_1])
                    assert len(bad_domains) == 1
                    bad_domains[0].domain_contact_errors = str(p.d.uniprot_domain_pair_id) + ': ' + str(type(e))
                    error_handling_message = (
                        'Adding error to the domain pair with cath_id_1 {cath_id_1} '
                        'and cath_id_2 {cath_id_2} and running the pipeline again...'
                        .format(cath_id_1=p.t.cath_id_1, cath_id_2=p.t.cath_id_2))
                # Add an error message to the model, while keeping previous error messages
                empty_model.model_errors = (
                    '{}'.format(type(e))
                    if (not p.m or not p.m.model_errors)
                    else '{}; {}'.format(p.m.model_errors, type(e)))
                # Run the pipeline again if the error could be fixed with a better template,
                # or otherwise move on to the next model while saving the error message
                if (self.number_of_tries[p_idx] < 3 and
                        ((isinstance(e, ModellerError) and 'exceeded max_molpdf' in e.message) or
                        (isinstance(e, errors.ChainsNotInteracting)))):
                    self.db._update_rows(bad_domains)
                    self.db.add_uniprot_model(empty_model, p.d.path_to_data)
                    self.log.error(error_handling_message)
                    self.log.error(empty_model.model_errors)
                    self.log.debug(
                        'Obtained a bad model on try no {}. Rerunning the pipeline to find better templates.\n\n\n'
                        .format(self.number_of_tries[p_idx]))
                    self.log.debug('*' * 80)
                    self.number_of_tries[p_idx] += 1
                    return self.__call__(self.uniprot_id, self.mutations, self.run_type, self.n_cores, self.number_of_tries)
                elif self.number_of_tries[p_idx] >= 3:
                    empty_model.model_errors = empty_model.model_errors + '; Could not resolve error after 3 tries. Giving up'
                self.log.error(empty_model.model_errors)
                p.m = empty_model
            # Add either the empty model or the calculated model to the database
            self.log.info('Adding model.\n\n\n')
            self.db.add_uniprot_model(p.m, p.d.path_to_data)
        self.log.info('Finished processing all models for {} {}\n'.format(self.uniprot_id, self.mutations))


    def _compute_provean(self):
        """
        Align uniprot domains to structural domains in pdbfam, find the best
        template, and expand domain boundaries
        """

        # Go over all domains and domain pairs for a given protein
        for p in self.protein_definitions:
            self.__print_header(p)

            # Make a provean supset for the template
            if ( p.t == None or
                    (p.t.alignment_filename == None
                    if isinstance(p.t, sql_db.UniprotDomainTemplate)
                    else (p.t.alignment_filename_1 == None or p.t.alignment_filename_2 == None)) ):
                self.log.info('No template availible. Skipping...')
                continue
            if isinstance(p.t, sql_db.UniprotDomainTemplate) and not p.t.provean_supset_filename:
                try:
                    p.t.provean_supset_filename, p.t.provean_supset_length \
                        = self.get_template.build_provean_supporting_set(p.d, p.t)
                except errors.ProveanResourceError as e:
                    self.log.critical(e.message)
                    self.__clear_provean_temp_files()
                    os.killpg(e.child_process_group_id, signal.SIGTERM) # Send the signal to all the process groupss
                    self.log.critical('You should be dead now...')
                except errors.ProveanError as e:
                    p.t.template_errors = (
                        '{}'.format(type(e))
                        if (p.t.template_errors == None)
                        else '{}; {}'.format(p.t.template_errors, type(e)))
                    self.log.error(p.t.template_errors)
                    self.log.error('Skipping...')
                    continue
                finally:
                    self.__clear_provean_temp_files()
                self.log.info('Adding template with provean info...')
                self.db.add_uniprot_template(p.t, p.d.path_to_data)

            # Get a provean score for each mutation
            if not self.mutations:
                self.log.info('No mutations specified.')
                continue
            if not p.m:
                self.log.error('No model availible. Skipping...')
                continue
            for mutation in self.mutations.split(','):
                precalculated_mutations = self.db.get_uniprot_mutation(p.m, self.uniprot_id, mutation, p.d.path_to_data)
                self.log.info(precalculated_mutations)
                for precalculated_mutation in precalculated_mutations:
                    if (precalculated_mutation and
                            precalculated_mutation.provean_score != None):
                        self.log.info('Mutation already has a provean score. Skipping...')
                        continue
                    try:
                        mut_data = self.get_mutation.get_mutation_data(p.d, p.t, p.m, self.uniprot_id, mutation)
                    except (errors.MutationOutsideDomain,
                            errors.MutationOutsideInterface) as e:
                        self.log.debug('{}: {}; OK'.format(type(e), e.message))
                        continue
                    if mut_data.path_to_provean_supset == None:
                        self.log.error(
                            'No provean supset pre-calculated! Skipping... '
                            'But normally this should not happen! Check the log for provean errors.')
                        continue
                    provean_mutation, provean_score = \
                        self.get_mutation.get_provean_score(
                            mut_data.uniprot_domain_id,
                            mut_data.mutation_domain,
                            mut_data.domain_sequences[0],
                            mut_data.path_to_provean_supset)
                    precalculated_mutation.provean_score = provean_score
                    self.log.debug('Provean mutation: {}'.format(provean_mutation))
                    self.log.debug('Provean score: {}'.format(provean_score))
                    self.log.info('Adding provean to mutation {}...\n\n\n'.format(mutation))
                    self.db.add_uniprot_mutation(precalculated_mutation, p.d.path_to_data)
        self.log.info(
            'Finished adding provean info to all templates and mutations in: -u {} -m {}'
            .format(self.uniprot_id, self.mutations))


    def _compute_mutations(self):
        """
        """
        possible_mutation_errors = (
            errors.PDBError,
            errors.FoldXError,
            errors.ResourceError,)

        for p in self.protein_definitions:
            self.__print_header(p)
            if not self.mutations:
                self.log.debug('Not evaluating mutations because no mutations specified...')
                continue
            if p.m == None or p.m.model_filename == None:
                self.log.debug('Skipping because no model...')
                continue
            self.log.debug(self.mutations)
            for mutation in self.mutations.split(','):
                self.log.debug(mutation)
                self.log.debug('-' * 80)
                precalculated_mutation = self.db.get_uniprot_mutation(p.m, self.uniprot_id, mutation, p.d.path_to_data)
                self.log.info(precalculated_mutation)
                if precalculated_mutation:
                    p.mut.append(precalculated_mutation)
                    self.log.info('Mutation already calculated. Skipping...')
                    continue

                # Construct empty models that will be used if we have errors
                if (isinstance(p.d, sql_db.UniprotDomain) and
                        p.t.domain_def == None):
                    uniprot_mutation = sql_db.UniprotDomainMutation()
                    uniprot_mutation.uniprot_domain_id = p.d.uniprot_domain_id
                    uniprot_mutation.uniprot_id = self.uniprot_id
                    uniprot_mutation.mutation = mutation
                elif (isinstance(p.d, sql_db.UniprotDomainPair) and
                        (p.t.domain_def_1 == None or p.t.domain_def_2 == None)):
                    uniprot_mutation = sql_db.UniprotDomainPairMutation()
                    uniprot_mutation.uniprot_domain_pair_id = p.d.uniprot_domain_pair_id
                    uniprot_mutation.uniprot_id = self.uniprot_id
                    uniprot_mutation.mutation = mutation

                try:
                    mut_data = self.get_mutation.get_mutation_data(p.d, p.t, p.m, self.uniprot_id, mutation)
                    uniprot_mutation = self.get_mutation.evaluate_mutation(mut_data)
                except (errors.MutationOutsideDomain,
                        errors.MutationOutsideInterface) as e:
                    self.log.debug('{}: {}; OK'.format(type(e), e.message))
                    continue
                except possible_mutation_errors as e:
                    uniprot_mutation.mutation_errors = '{}: {}'.format(type(e), e.message)
                    self.log.debug(uniprot_mutation.mutation_errors)
                self.log.info('Adding mutation {}\n\n\n'.format(mutation))
                p.mut.append(uniprot_mutation)
                self.db.add_uniprot_mutation(uniprot_mutation, p.d.path_to_data)
        self.log.info('Finished processing all mutations for {} {}\n'.format(self.uniprot_id, self.mutations))


    def __print_header(self, p):
        self.log.info(
            p.d.uniprot_domain_id
            if isinstance(p.d, sql_db.UniprotDomain)
            else p.d.uniprot_domain_pair_id)
        self.log.info('{}\t{}\t{}'.format(p.d, p.t, p.m,))


    def __clear_provean_temp_files(self):
        subprocess.check_call('rm -rf ' + self.provean_temp_path + '/provean*', shell=True)
        subprocess.check_call('rm -rf ' + self.provean_temp_path + '/*cdhit*', shell=True)


    def __prepare_temp_folder(self):
        # create the basic tmp directory
        # delete its content if it exists (AS: disabled so that I can continue from previous run)
        if not os.path.isdir(self.temp_path):
            subprocess.check_call('mkdir -p ' + self.temp_path, shell=True)

        # unique_temp_folder
        if not os.path.isdir(self.unique_temp_folder):
            subprocess.check_call('mkdir -p ' + self.unique_temp_folder, shell=True)

        # tcoffee
        if not os.path.isdir(self.unique_temp_folder + '/tcoffee'):
            # create tmp for tcoffee
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/tcoffee && ' + \
                            'mkdir -p ' + self.unique_temp_folder + '/tcoffee/tmp && ' + \
                            'mkdir -p ' + self.unique_temp_folder + '/tcoffee/lck && ' + \
                            'mkdir -p ' + self.unique_temp_folder + '/tcoffee/cache'
            subprocess.check_call(mkdir_command, shell=True)

        # FoldX
        if not os.path.isdir(self.unique_temp_folder + '/FoldX'):
            # make the directories
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/FoldX'
            # copy the executables
            cp_command = 'cp ' + self.bin_path + 'FoldX.linux64 ' + self.unique_temp_folder + '/FoldX/ && ' + \
                       'cp ' + self.bin_path + 'rotabase.txt ' + self.unique_temp_folder + '/FoldX/'
            subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
            # Copy dssp into the folder for modelling
            cp_command = 'cp ' + self.bin_path + 'dssp-2.0.4-linux-amd64 ' + self.unique_temp_folder + '/FoldX/'
            subprocess.check_call(cp_command, shell=True)

        # modeller
        if not os.path.isdir(self.unique_temp_folder + '/modeller'):
            # create workingfolder for modeller
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/modeller'
            subprocess.check_call(mkdir_command, shell=True)
            # Copy knot into the same folder as modeller
            cp_command = 'cp ' + self.bin_path + 'topol ' + self.unique_temp_folder + '/modeller'
            subprocess.check_call(cp_command, shell=True)

        # sequence conservation
        if not os.path.isdir(self.unique_temp_folder + '/sequence_conservation'):
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/sequence_conservation'
            subprocess.check_call(mkdir_command, shell=True)
            # provean
            cp_command = 'cp ' + self.bin_path + 'provean ' + self.unique_temp_folder + '/sequence_conservation/'
            subprocess.check_call(cp_command, shell=True)

        # analyze_structure
        if not os.path.isdir(self.unique_temp_folder + '/analyze_structure'):
            # create workingfolder for analyzing structure sasa and secondary structure
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/analyze_structure'
            subprocess.check_call(mkdir_command, shell=True)
            # Pops
            cp_command = 'cp ' + self.bin_path + 'pops ' + self.unique_temp_folder + '/analyze_structure'
            subprocess.check_call(cp_command, shell=True)
            # Dssp
            cp_command = 'cp ' + self.bin_path + 'dssp-2.0.4-linux-amd64 ' + self.unique_temp_folder + '/analyze_structure/dssp'
            subprocess.check_call(cp_command, shell=True)
#                # Naccess
#                cp_command = (
#                    'cp ' + self.bin_path + 'naccess ' + self.unique_temp_folder + '/analyze_structure/ && '
#                    'cp ' + self.bin_path + 'accall ' + self.unique_temp_folder + '/analyze_structure/ && '
#                    'cp ' + self.bin_path + 'standard.data ' + self.unique_temp_folder + '/analyze_structure/ && '
#                    'cp ' + self.bin_path + 'vdw.radii ' + self.unique_temp_folder + '/analyze_structure/')
#                subprocess.check_call(cp_command, shell=True)
            #MSMS
            cp_command = (
                'cp ' + self.bin_path +  'pdb_to_xyzrn ' +  self.unique_temp_folder + '/analyze_structure/ && '
                'cp ' + self.bin_path + 'standard.data ' +  self.unique_temp_folder + '/analyze_structure/ && '
                'cp ' + self.bin_path + 'atmtypenumbers ' +  self.unique_temp_folder + '/analyze_structure/ && '
                'cp ' + self.bin_path + 'msms.x86_64Linux2.2.6.1 ' + self.unique_temp_folder + '/analyze_structure/')
            subprocess.check_call(cp_command, shell=True)



if __name__ == '__main__':
    # read which configFile to use
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True)
    parser.add_argument('-i', '--input_file')
    parser.add_argument('-u', '--uniprot_id')
    parser.add_argument('-m', '--mutations', nargs='?', default=['',])
    parser.add_argument('-t', '--run_type', nargs='?', type=int, default=1)
    parser.add_argument('-n', '--n_cores', nargs='?', type=int, default=1)
    args = parser.parse_args()

    assert os.path.isfile(args.config_file)
    pipeline = Pipeline(args.config_file)

    if args.input_file \
    and os.path.isfile(args.input_file):
        uniprot_ids = []
        mutations = []
        with open(args.input_file, 'r') as fh:
            for line in fh:
                # Can skip lines by adding spaces or tabs before them
                if line[0][0] == ' ' or line[0][0] == '\t':
                    continue

                row = [ l.strip() for l in line.split('\t') ]

                # AS: Mutation does not necessarily have to be specified
                if len(row) > 1:
                    uniprot_id, mutation = row[0], row[1]
                elif len(row) == 1:
                    uniprot_id = row[0]
                    mutation = ''
                #
                uniprot_ids.append(uniprot_id)
                mutations.append(mutation)

    elif args.uniprot_id:
        uniprot_ids = [args.uniprot_id,]
        mutations = [args.mutations,] if args.mutations is not None else ['',]

    else:
        raise Exception('Need to supply either a list of uniprot_mutation combos '
            'or a flatfile with the same!')

    run_type = args.run_type
    n_cores = args.n_cores
    # Run jobs
    print uniprot_ids
    print mutations
    print run_type
    for uniprot_id, mutation in zip(uniprot_ids, mutations):
        print uniprot_id
        print mutation
        pipeline(uniprot_id, mutation, run_type, n_cores)



