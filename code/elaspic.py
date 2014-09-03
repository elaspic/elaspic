# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 18:37:36 2012

@author: niklas
"""
import os
import time
import subprocess
import tempfile
import argparse
import atexit
import signal

from Bio.SubsMat import MatrixInfo
from ConfigParser import SafeConfigParser

import domain_alignment
import domain_model
import domain_mutation

import helper_functions as hf
import errors
import sql_db


class Pipeline(object):

    def __init__(self, configFile):
        """
        """
        # Read the configuration file and set the variables
        configParser = SafeConfigParser(
            defaults={
                'db_type': sql_db.SQL_FLAVOUR,
                'path_to_archive': '/home/kimlab1/database_data/elaspic_v2/',
                'web_server': False,
            })
        configParser.read(configFile)

        # From [DEFAULT]
        self.global_temp_path = configParser.get('DEFAULT', 'global_temp_path')
        temp_path_suffix = configParser.get('DEFAULT', 'temp_path').strip('/') + '/' # So it never starts with '/';
        self.debug = configParser.getboolean('DEFAULT', 'DEBUG')
        self.path_to_archive = configParser.get('DEFAULT', 'path_to_archive')
        self.look_for_interactions = configParser.getboolean('DEFAULT', 'look_for_interactions')
        self.n_cores = configParser.get('DEFAULT', 'n_cores')
        self.db_type = configParser.get('DEFAULT', 'db_type')
        self.db_path = configParser.get('DEFAULT', 'db_path')
        self.db_is_immutable = True if self.db_type.lower().startswith('sqlite') else False
        self.web_server = configParser.get('DEFAULT', 'web_server')

        # From [SETTINGS]
        self.blast_db_path = configParser.get('SETTINGS', 'blast_db_path')
        self.pdb_path = configParser.get('SETTINGS', 'pdb_path')
        self.output_path = configParser.get('SETTINGS', 'output_path')
        self.bin_path = configParser.get('SETTINGS', 'bin_path')
        self.max_runtime = configParser.get('SETTINGS', 'max_runtime')

        # From [GET_MODEL]
        self.modeller_runs = configParser.getint('GET_MODEL', 'modeller_runs')

        # From [GET_MUTATION]
        self.foldx_water = configParser.get('GET_MUTATION', 'foldx_water')
        self.foldx_num_of_runs = configParser.get('GET_MUTATION', 'foldx_num_of_runs')
        self.matrix_type = configParser.get('GET_MUTATION', 'matrix_type')
        self.gap_start = configParser.getint('GET_MUTATION', 'gap_start')
        self.gap_extend = configParser.getint('GET_MUTATION', 'gap_extend')
        self.matrix = getattr(MatrixInfo, self.matrix_type)

        # Get the temporary folder, taking into account the $TMPDIR env variable
        self.temp_path = hf.get_temp_path(self.global_temp_path, temp_path_suffix)
        # Make sure that the correct temporary folder is created when running as
        # a job on a cluster
        hostname = hf.get_hostname()
        if (self.temp_path.startswith(os.path.join(self.global_temp_path, temp_path_suffix)) and
            any([(x.lower() in hostname) for x in ['node', 'behemoth', 'grendel', 'beagle']])):
            raise Exception('You should be using a temp folder that it specific to the particular job!')

        # Copy the sqlite version of the pdbfam database if using sqlite3
#        if self.db_type == 'sqlite_file':
#            self.__copy_sqlite_database()

        # Copy the blast databaes
        self.__copy_blast_database()


    def __copy_sqlite_database(self):
        db_filename = self.db_path.strip().split('/')[-1]
        if not os.path.isfile(self.temp_path + db_filename):
            system_command = 'cp -u ' + self.temp_path + ' ' + self.db_path + '/'
            print system_command
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, error_message = childProcess.communicate()
            if childProcess.returncode != 0:
                print result
                print error_message
                raise Exception('Could not copy the sqlite3 database!')


    def __copy_blast_database(self):
        #######################################################################
        # If running on the cluster copy the database to the tmp DIR for
        # speedup and to avoid killing the network. BLAST is very IO intensive
        # and you don't want that to be run over the network!
        # I can distinguish where the pipeline is running by checking the username.
        # You will have to adjust that!
        # My usernames are:
        # local: niklas
        # banting: nberliner
        # Scinet: joan
        username  = hf.get_username()
        hostname = hf.get_hostname()
        if username == 'strokach' and (('beagle' in hostname) or ('grendel' in hostname)):
            # The blast database should already have been copied to the local temp folder
            if not os.path.isdir(self.global_temp_path + 'blast/db/'):
#                raise Exception(
#                    'Could not find a blast database even though it should be '
#                    'present on these nodes.')
                print ('Creating a simlink to the blast folder even though the '
                    'blast database should be copied on these nodes')
                system_command = (
                    'mkdir -p ' + self.global_temp_path + 'blast/ && ' +
                    'ln -sf ' + self.blast_db_path + 'db ' +  self.global_temp_path + 'blast/')
            else:
                return
#            # Copy the blast database if running on beagle or Grendel
#            system_command = (
#                # 'mkdir -p ' + self.global_temp_path + 'blast/pdbaa_db/ && ' +
#                'mkdir -p ' + self.global_temp_path + 'blast/db/ && ' +
#                # 'rsync -rzu ' + self.blast_db_path + 'pdbaa_db/ ' + self.global_temp_path + 'blast/pdbaa_db/ && ' +
#                'rsync -rzu ' + self.blast_db_path + 'db/ ' + self.global_temp_path + 'blast/db/')
        elif username == 'strokach' or username == 'alexey':
            # Use a symbolic link to the blast database and use the home folder for temporary storage
            system_command = (
                #'rm -rf ' + self.global_temp_path + 'blast/ && ' +
                'mkdir -p ' + self.global_temp_path + 'blast/ && ' +
                #'ln -sf ' + self.blast_db_path + 'pdbaa_db ' + self.global_temp_path + 'blast/ && ' +
                'ln -sf ' + self.blast_db_path + 'db ' +  self.global_temp_path + 'blast/')
        elif username == 'witvliet' or username == 'kimadmin' or username == 'www-data':
            system_command = (
                'rm -rf ' + self.global_temp_path + 'blast/ && ' +
                'mkdir -p ' + self.global_temp_path + 'blast/ && ' +
                'ln -sf ' + self.blast_db_path + 'pdbaa_db ' + self.global_temp_path + 'blast/ && ' +
                'ln -sf ' + self.blast_db_path + 'db ' +  self.global_temp_path + 'blast/')
        elif username == 'joan':
            # For scinet, blast is already installed, but the database needs to be copied
            system_command = (
                'mkdir -p ' + self.global_temp_path + 'blast && ' +
                'cp -ru $HOME/niklas-pipeline/blastdb/pdbaa_db ' +
                self.global_temp_path + 'blast/')

        print system_command
        # Try running the system command several times in case it doesn't work
        # the first time
        n_tries = 0
        childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, error_message = childProcess.communicate()
        rc = childProcess.returncode
        while rc != 0 and n_tries < 10:
            print "Couldn't copy the blast database. Retrying..."
            time.sleep(15)
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if rc != 0:
            print result
            print error_message
            raise Exception(
                'Couldn\'t copy the blast database on {hostname} as {username}!'
                .format(hostname=hostname, username=username))


    def __call__(self, uniprot_id, mutations, run_type=1, n_cores=None, number_of_tries=[], logger=None):
        """ Run the main function of the program and parse errors
        """
        self.uniprot_id = uniprot_id
        self.mutations = mutations
        self.calculated_mutations = []
        self.run_type = run_type
        if n_cores is not None:
            self.n_cores = n_cores
        self.number_of_tries = number_of_tries

        self.unique = tempfile.mkdtemp(prefix='', dir=self.temp_path).split('/')[-1]
        self.unique_temp_folder = self.temp_path + self.unique + '/'

        if logger is None:
            logger = hf.get_logger(self.debug)
        self.logger = logger
        self.logger.info(self.unique_temp_folder)

        # Switch to the root of the unique tmp directory
        os.chdir(self.unique_temp_folder)

        # Create subfolders for temporary storage
        self.__prepare_temp_folder()

        # Set up a temporary folder for provean
        self.provean_temp_path = self.__prepare_provean_temp_folder()
        atexit.register(self.__clear_provean_temp_files)

        # Initialise the sql database for accessing all information
        self.logger.info("Connecting to a '{}' database...".format(self.db_type))
        if self.db_type.lower() == 'sqlite_file':
            self.logger.info('Path to the database: {}'.format(self.db_path))
        self.db = sql_db.MyDatabase(
            path_to_sqlite_db=self.db_path, sql_flavour=self.db_type,
            temp_path=self.temp_path, path_to_archive=self.path_to_archive,
            is_immutable=self.db_is_immutable, logger=self.logger)

        # Obtain all domains and interactions for a given uniprot
        self.logger.info('Obtaining protein domain information...')
        self.uniprot_domains = self.db.get_uniprot_domain(self.uniprot_id, True)
        if not self.uniprot_domains:
            self.logger.error('Uniprot {} has no pfam domains'.format(self.uniprot_id))
            return
        self._update_path_to_data(self.uniprot_domains)

        # Mutations
        self.uniprot_mutations = []

        # Find provean
        if run_type in [1, 5]:
            self.logger.info('\n\n\n' + '*' * 80)
            self.logger.info("Computing provean...")
            if self._compute_provean():
                if run_type == 1:
                    self.logger.info('Finished run_type {}'.format(run_type))
                    return
                # If provean was updated, we need to reload uniprot domain data
                # for all the other domains
                self.logger.info('\n\n\n')
                self.logger.info('Obtaining protein domain information...')
                self.uniprot_domains = self.db.get_uniprot_domain(self.uniprot_id, True)

        # Get interactions
        if self.look_for_interactions:
            self.logger.info('Obtaining protein domain pair information...')
            self.uniprot_domain_pairs = self.db.get_uniprot_domain_pair(self.uniprot_id, True)
            self._update_path_to_data(self.uniprot_domain_pairs)

        # Make models
        if run_type in [2, 4, 5]:
            self.get_model = domain_model.GetModel(
                self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
                self.db, self.logger, self.n_cores, self.modeller_runs)

            self.logger.info('\n\n\n' + '*' * 80)
            self.logger.info("Building models...")
            self._compute_models()
        
        # Analyse mutations
        if run_type in [3, 4, 5] and self.mutations:
            self.get_mutation = domain_mutation.GetMutation(
                self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
                self.db, self.logger, self.n_cores, self.bin_path, self.foldx_water,
                self.foldx_num_of_runs, self.matrix, self.gap_start, self.gap_extend,)
            
            self.logger.info('\n\n\n' + '*' * 80)
            self.logger.info("Analyzing mutations...")
            self._compute_mutations()


    def _update_path_to_data(self, d_list):
        if not isinstance(d_list, list):
            d_list = [d_list]
        for d in d_list:
            if not d.path_to_data or any([len(x) > 255 for x in d.path_to_data.split('/')]):
                d.path_to_data = hf.get_uniprot_base_path(d) + hf.get_uniprot_domain_path(d)
                self.db.merge_row(d)
                subprocess.check_call('mkdir -p {}'.format(self.temp_path + d.path_to_data), shell=True)


    def _compute_provean(self):
        """
        """
        d = self.uniprot_domains[0] # All uniprot domains will have the same uniprot_sequence info
        self.__print_header(d)

        if (d.uniprot_sequence.provean and
                d.uniprot_sequence.provean.provean_supset_filename and
                os.path.isfile(
                    self.path_to_archive + hf.get_uniprot_base_path(d) +
                    d.uniprot_sequence.provean.provean_supset_filename) ):
            self.logger.debug('The provean supset has already been calculated. Done!')
            return None
        elif (d.uniprot_sequence.provean and
                d.uniprot_sequence.provean.provean_supset_filename):
            self.logger.debug('Provean has been calculated but the file is missing!!!!!!! Recalculating...')
            provean = d.uniprot_sequence.provean
        elif d.uniprot_sequence.provean:
            provean = d.uniprot_sequence.provean
        else:
            provean = sql_db.Provean()
            provean.uniprot_id = d.uniprot_id

        try:
            provean.provean_supset_filename, provean.provean_supset_length = \
                domain_alignment.build_provean_supporting_set(
                    self, d.uniprot_id, d.uniprot_sequence.uniprot_name,
                    d.uniprot_sequence.uniprot_sequence)

        except (errors.ProveanError, errors.ProveanResourceError) as e:
            provean.provean_errors = self.__add_new_error(provean.provean_errors, e)
            self.logger.error(provean.provean_errors)
            self.logger.error(e.__str__())
            if isinstance(e, errors.ProveanResourceError):
                # Send the kill signal to the main process group, killing everything
                self.__clear_provean_temp_files() # these won't get cleaned once the process dies
                self.logger.error('Killing group...')
                os.killpg(e.child_process_group_id, signal.SIGTERM)
                self.logger.critical(
                    'Provean ran out of resources. Everything has been killed. '
                    'The user will probably not see this message...')

        self.__clear_provean_temp_files()
        self.db.merge_provean(provean, hf.get_uniprot_base_path(d))
        self.logger.info('Finished computing provean for {}'.format(d.uniprot_id))
        return provean


    def _compute_models(self):
        """
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        possible_model_errors = (
            errors.ModellerError,
            errors.PDBChainError,
            errors.ChainsNotInteractingError,
            errors.MutationOutsideDomainError,
            errors.MutationOutsideInterfaceError,
            errors.NoSequenceFound,
            errors.TcoffeeError,
            errors.PDBNotFoundError)

        # Go over all domains and domain pairs for a given protein
        for d in self.uniprot_domains + self.uniprot_domain_pairs:
            self.__print_header(d)
            
            # Check if we should skip the model
            if (isinstance(d, sql_db.UniprotDomain) and
                    not (d.template and d.template.cath_id) ):
                self.logger.error('No structural template availible for this domain. Skipping...')
                continue
            elif (isinstance(d, sql_db.UniprotDomainPair) and
                    not (d.template and d.template.cath_id_1 and d.template.cath_id_2) ):
                self.logger.error('No structural template availible for this domain pair. Skipping...')
                continue
            elif d.template.model and d.template.model.model_filename:
                self.logger.info('Model already calculated. Skipping...')
                continue
            elif (d.template.model and d.template.model.model_errors and
                    (('Giving up' in d.template.model.model_errors) or
                    (d.template.model.model_errors.count(';') > 2)) ):
                self.logger.info(
                    'Previous model had unfixable errors: "{}". Skipping...'
                    .format(d.template.model.model_errors))
                continue

            # Make a homology model using the domain and template information
            try:
                self.get_model(d)

            except possible_model_errors as e:
                # Find domains that were used as a template and eventually led to
                # the error in the model, and add the error to their `domain_errors`
                # or `domain_contact_errors` fields.
                self.logger.error(str(e))
                if isinstance(d, sql_db.UniprotDomain):
                    if d.template.model == None:
                        d.template.model = sql_db.UniprotDomainModel()
                        d.template.model.uniprot_domain_id = d.uniprot_domain_id
                    bad_domains = self.db.get_rows_by_ids(
                        sql_db.Domain, [sql_db.Domain.cath_id], [d.template.cath_id])
                    assert len(bad_domains) == 1
                    bad_domain = bad_domains[0]
                    bad_domain.domain_errors = str(d.uniprot_domain_id) + ': ' + str(type(e))
                    self.logger.error(
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
                    self.logger.error(
                        "Adding error '{0}' to the domain pair with cath_id_1 {1} "
                        "and cath_id_2 {2}..."
                        .format(bad_domain.domain_contact_errors, d.template.cath_id_1, d.template.cath_id_2))
                # Add the error type to the model_errors column
                d.template.model.model_errors = self.__add_new_error(d.template.model.model_errors, e)
                self.logger.error(d.template.model.model_errors)
                self.db.merge_row(bad_domain)
                # d.template.model = empty_model

            # Add either the empty model or the calculated model to the database
            self.logger.info('Adding model...')
            self.db.merge_model(d, d.path_to_data)
        self.logger.info('Finished processing all models for {} {}'.format(self.uniprot_id, self.mutations))


    def _compute_mutations(self):
        """
        """
        possible_mutation_errors = (
            errors.PDBError,
            errors.FoldxError,
            errors.ResourceError,)

        for d in self.uniprot_domains + self.uniprot_domain_pairs:
            self.__print_header(d)

            if not self.mutations:
                self.logger.debug('Not evaluating mutations because no mutations specified...')
                continue
            elif ( (isinstance(d, sql_db.UniprotDomain) and
                    not d.template.domain) or
                    (isinstance(d, sql_db.UniprotDomainPair) and
                    not (d.template.domain_1 and d.template.domain_2)) ):
                self.logger.debug('Skipping because no structural template is availible...')
                continue
            elif not d.template.model or not d.template.model.model_filename:
                self.logger.debug('Skipping because no model...')
                continue

            self.logger.debug('Going over all mutations: {}'.format(self.mutations))
            for mutation in self.mutations.split(','):
                self.logger.debug('-' * 80)
                self.logger.debug('Mutation: {}'.format(mutation))

                # Check to see if we have a precalculated mutation. Skip if all
                # parameters have been calculated; otherwise analyse the remaining
                # parameters. Create an empty mutation if the mutation has not
                # been precalculated.
                precalculated_mutation = self.db.get_uniprot_mutation(d, mutation, self.uniprot_id, True)
                self.logger.info('Have the following precalculated mutation: {}'.format(precalculated_mutation))
                if (precalculated_mutation and
                        (precalculated_mutation.provean_score and
                        precalculated_mutation.stability_energy_wt and
                        precalculated_mutation.ddg != None)):
                    self.calculated_mutations.append(precalculated_mutation)
                    self.logger.info('Mutation has already been completely evaluated. Skipping...')
                    continue
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
                    self.logger.debug('{}: {}; OK'.format(type(e), e))
                    continue
                except possible_mutation_errors as e:
                    uniprot_mutation.mutation_errors = '{}: {}'.format(type(e), e)
                    self.logger.debug(uniprot_mutation.mutation_errors)

                self.logger.info('Adding mutation {}'.format(mutation))
                self.calculated_mutations.append(uniprot_mutation)
                self.db.merge_mutation(uniprot_mutation, d.path_to_data)
                self.uniprot_mutations.append(uniprot_mutation)
        self.logger.info('Finished processing all mutations for {} {}'.format(self.uniprot_id, self.mutations))


    def __print_header(self, d):
        # self.logger.info('Domain or domain pair number: {}'.format(d_idx))
        self.logger.debug('=' * 80)
        if isinstance(d, sql_db.UniprotDomain):
            self.logger.debug('uniprot_domain_id: {}'.format(d.uniprot_domain_id))
        else:
            self.logger.debug('uniprot_domain_pair_id: {}'.format(d.uniprot_domain_pair_id))


    def __add_new_error(self, d_error_log, e):
        if d_error_log is None:
            return str(type(e))
        else:
            return '{}; {}: {}'.format(d_error_log, type(e), str(e))


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
            subprocess.check_call("mkdir -p '{}'".format(self.unique_temp_folder + '/FoldX'), shell=True)
            # copy the executables
            subprocess.check_call("cp '{}' '{}'".format(
                self.bin_path + 'foldx64Linux', self.unique_temp_folder + '/FoldX/'), shell=True)
            subprocess.check_call("cp '{}' '{}'".format(
                self.bin_path + 'FoldX.linux64', self.unique_temp_folder + '/FoldX/'), shell=True)
            subprocess.check_call("cp '{}' '{}'".format(
                self.bin_path + 'rotabase.txt', self.unique_temp_folder + '/FoldX/'), shell=True)
            # Copy dssp into the folder for modelling
#            cp_command = 'cp ' + self.bin_path + 'mkdssp ' + self.unique_temp_folder + '/FoldX/'
#            subprocess.check_call(cp_command, shell=True)

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
            # cp_command = 'cp ' + self.bin_path + 'dssp-2.0.4-linux-amd64 ' + self.unique_temp_folder + '/analyze_structure/dssp'
            # cp_command = 'cp ' + self.bin_path + 'mkdssp ' + self.unique_temp_folder + '/analyze_structure/dssp'
            # subprocess.check_call(cp_command, shell=True)
            cp_command = 'cp ' + self.bin_path + 'stride ' + self.unique_temp_folder + '/analyze_structure/stride'
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


    def __prepare_provean_temp_folder(self):
        username  = hf.get_username()
        hostname = hf.get_hostname()
        if username == 'strokach' and not (('beagle' in hostname) or ('grendel' in hostname)):
            # Can't use /tmp for temporary storage on banting because I run out
            # of space and crash the node
            self.logger.debug('Using a temp folder on kimstg for provean temp files. This may lead to poor performace...')
            provean_temp_path = '/home/kimlab1/strokach/tmp/elaspic/' + self.unique + '/provean_temp/'
        else:
            provean_temp_path = self.unique_temp_folder + 'provean_temp/'
        subprocess.check_call('mkdir -p ' + provean_temp_path, shell=True)
        return provean_temp_path


    def __clear_provean_temp_files(self):
        subprocess.check_call('rm -rf ' + self.provean_temp_path + '/provean*', shell=True)
        subprocess.check_call('rm -rf ' + self.provean_temp_path + '/*cdhit*', shell=True)



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
        error_message = (
            'Need to supply either a list of uniprot_mutation combos '
            'or a flatfile with the same!')
        raise Exception(error_message)

    run_type = args.run_type
    n_cores = args.n_cores

    # Run jobs
    for uniprot_id, mutation in zip(uniprot_ids, mutations):
        print uniprot_id
        print mutation
        print run_type
        uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutation, run_type, n_cores)

