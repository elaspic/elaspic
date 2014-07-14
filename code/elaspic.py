# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 18:37:36 2012

@author: niklas
"""
import os
import sys
import time
import subprocess
import tempfile
import argparse
import atexit
import signal
import cPickle as pickle
import sqlalchemy.ext.serializer

from Bio.SubsMat import MatrixInfo
from ConfigParser import SafeConfigParser

from modeller import ModellerError

import domain_alignment
import domain_model
import domain_mutation

import helper_functions as hf
import errors
import sql_db

try:
    from celery.utils.log import get_task_logger
except ImportError:
    pass



class Pipeline(object):

    def __init__(self, configFile):
        """
        """
        # Read the configuration file and set the variables
        configParser = SafeConfigParser(
            defaults={
                'db_type': sql_db.sql_flavor,
                'db_path': hf.get_path_to_current_file() + '/../db/pipeline.db',
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

        # Copy the sqlite version of the pdbfam database if using sqlite3
        if self.db_type == 'sqlite_file':
            self.__copy_sqlite_database()

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
                sys.exit('Could not copy the sqlite3 database!')


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


    def __call__(self, uniprot_id, mutations, run_type=1, n_cores=None, number_of_tries=[]):
        """ Run the main function of the program and parse errors
        """
        print uniprot_id
        print mutation

        self.uniprot_id = uniprot_id
        self.mutations = mutations
        self.run_type = run_type
        if n_cores is not None:
            self.n_cores = n_cores
        self.number_of_tries = number_of_tries

        self.logger = hf.get_logger(self.debug)
        self.unique = tempfile.mkdtemp(prefix='', dir=self.temp_path).split('/')[-1]
        self.unique_temp_folder = self.temp_path + self.unique + '/'
        if self.web_server: # Webserver logging is handled by Celery
            self.logger = get_task_logger('web_pipeline.tasks')

        # Switch to the root of the unique tmp directory
        os.chdir(self.unique_temp_folder)

        # Create subfolders for temporary storage
        self.__prepare_temp_folder()

        # Set up a temporary folder for provean
        self.provean_temp_path = self.__prepare_provean_temp_folder()
        atexit.register(self.__clear_provean_temp_files)

        # Initialise the sql database for accessing all information
        self.logger.info('Initializing the database...')
        self.db = sql_db.MyDatabase(
            path_to_sqlite_db=self.temp_path+'pipeline.db', sql_flavor=self.db_type,
            path_to_temp=self.temp_path, path_to_archive=self.path_to_archive,
            is_immutable=False, logger=self.logger)

        # Initialize external class objects
        self.get_template = domain_alignment.GetTemplate(
            self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
            self.db, self.logger, self.n_cores,
            self.provean_temp_path)

        self.get_model = domain_model.GetModel(
            self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
            self.db, self.logger, self.n_cores, self.modeller_runs)

        self.get_mutation = domain_mutation.GetMutation(
            self.global_temp_path, self.temp_path, self.unique, self.pdb_path,
            self.db, self.logger, self.n_cores, self.bin_path, self.foldx_water,
            self.foldx_num_of_runs, self.matrix, self.gap_start, self.gap_extend,)

        # Obtain all domains and interactions for a given uniprot
        self.logger.info('Obtaining protein domain information...')
        self.uniprot_domains = self.db.get_uniprot_domain(self.uniprot_id, True)

        if not self.uniprot_domains:
            self.logger.error('Uniprot {} has no pfam domains'.format(self.uniprot_id))
            sys.exit('Uniprot {} has no pfam domains'.format(self.uniprot_id))

        if self.look_for_interactions:
            self.logger.info('Obtaining protein domain pair information...')
            self.uniprot_domain_pairs = self.db.get_uniprot_domain_pair(self.uniprot_id, True)

        # Set the path_to_data variable
        for d in self.uniprot_domains + self.uniprot_domain_pairs:
            if not d.path_to_data:
                d.path_to_data = hf.get_uniprot_base_path(d) + hf.get_uniprot_domain_path(d)
                self.db.merge_row(d)

        if run_type in [1, 5]:
            self.logger.info('-' * 80)
            self.logger.info("Computing provean...")
            self._compute_provean()
        if run_type in [2, 4, 5]:
            self.logger.info("Building models...")
            self._compute_models()
        if run_type in [3, 4, 5]:
            self.logger.info("Analyzing mutations...")
            self._compute_mutations()


    def _compute_provean(self):
        """
        """
        d = self.uniprot_domains[0] # All uniprot domains will have the same uniprot_sequence info
        self.__print_header(d)

        if d.uniprot_provean and d.uniprot_provean.provean_supset_filename:
            self.logger.debug('The provean supset has already been calculated. Done!')
            return

        if d.uniprot_provean:
            uniprot_provean = d.uniprot_provean
        else:
            uniprot_provean = sql_db.UniprotProvean()
            uniprot_provean.uniprot_id = d.uniprot_id

        try:
            uniprot_provean.provean_supset_filename, uniprot_provean.provean_supset_length = \
                domain_alignment.build_provean_supporting_set(self, d.uniprot_id, d.uniprot_name, d.uniprot_sequence.seq)

        except (errors.ProveanError, errors.ProveanResourceError) as e:
            uniprot_provean.provean_errors = self.__add_new_error(uniprot_provean.provean_errors, e)
            self.logger.error(uniprot_provean.provean_errors)
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
        self.db.merge_row(uniprot_provean)
        # Propagate the provean supporting set to all the other domains and domain pairs
        for d_2 in self.uniprot_domains:
            d_2.uniprot_provean = uniprot_provean
        for d_2 in self.uniprot_domain_pairs:
            if d_2.uniprot_domain_1.uniprot_id == d.uniprot_id:
                d_2.uniprot_domain_1.uniprot_provean = uniprot_provean
            if d_2.uniprot_domain_2.uniprot_id == d.uniprot_id:
                d_2.uniprot_domain_2.uniprot_provean = uniprot_provean
        self.logger.info('Finished computing provean for {}\n\n\n'.format(d.uniprot_id))



    def _compute_models(self):
        """
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        possible_model_errors = (
            ModellerError,
            errors.ModellerError,
            errors.PdbChainError,
            errors.ChainsNotInteractingError,
            errors.MutationOutsideDomainError,
            errors.MutationOutsideInterfaceError,
            errors.NoSequenceFound,
            errors.TcoffeeError,
            errors.PdbNotFoundError)

        # Go over all domains and domain pairs for a given protein
        for d in self.uniprot_domains:
            self.__print_header(d)

            # Check if we should skip the model
            if (isinstance(d, sql_db.UniprotDomain) and not d.pdb_id):
                self.logger.error('No structural template availible for this domain. Skipping...')
                continue
            elif (isinstance(d, sql_db.UniprotDomainPair) and not d.pdb_id):
                self.logger.error('No structural template availible for this domain pair. Skipping...')
                continue
            elif d.model and d.model.model_filename:
                self.logger.info('Model already calculated. Skipping...')
                continue
            elif d.model and ('Giving up' in d.model.model_errors):
                self.logger.info(
                    'Previous model had unfixable errors: {}. Skipping...'
                    .format(d.model.model_errors))
                continue

            # Make a homology model using the domain and template information
            try:
                self.get_model(d)

            except possible_model_errors as e:
                # Find domains that were used as a template and eventually led to
                # the error in the model, and add the error to their `domain_errors`
                # or `domain_contact_errors` fields.
                if isinstance(d, sql_db.UniprotDomain):
                    empty_model = sql_db.UniprotDomainModel()
                    empty_model.uniprot_domain_id = d.uniprot_domain_id
                    bad_domains = self.db.get_rows_by_ids(
                        sql_db.Domain, [sql_db.Domain.cath_id], [d.cath_id])
                    assert len(bad_domains) == 1
                    bad_domain = bad_domains[0]
                    bad_domain.domain_errors = str(d.uniprot_domain_id) + ': ' + str(type(e))
                    self.logger.error(
                        'Adding error to the domain with cath_id {cath_id}...'
                        .format(cath_id=d.cath_id))
                elif isinstance(d, sql_db.UniprotDomainPair):
                    empty_model = sql_db.UniprotDomainPairModel()
                    empty_model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
                    bad_domains = self.db.get_rows_by_ids(
                        sql_db.DomainContact,
                        [sql_db.DomainContact.cath_id_1, sql_db.DomainContact.cath_id_2],
                        [d.cath_id_1, d.cath_id_2])
                    if len(bad_domains) == 0:
                        bad_domains = self.db.get_rows_by_ids(
                            sql_db.DomainContact,
                            [sql_db.DomainContact.cath_id_1, sql_db.DomainContact.cath_id_2],
                            [d.cath_id_2, d.cath_id_1])
                    assert len(bad_domains) == 1
                    bad_domain = bad_domains[0]
                    bad_domain.domain_contact_errors = str(d.uniprot_domain_pair_id) + ': ' + str(type(e))
                    self.logger.error(
                        'Adding error to the domain pair with cath_id_1 {cath_id_1} '
                        'and cath_id_2 {cath_id_2}...'
                        .format(cath_id_1=d.cath_id_1, cath_id_2=d.cath_id_2))
                # Add the error type to the model_errors column
                empty_model.model_errors = self.__add_new_error(empty_model.model_errors, e)
                self.logger.error(empty_model.model_errors)
                self.db.merge_row(bad_domain)
                d.model = empty_model

            # Add either the empty model or the calculated model to the database
            self.logger.info('Adding model...\n\n\n')
            self.db.merge_domain(d) # To save the alignments
            self.db.merge_model(d.model)
        self.logger.info('Finished processing all models for {} {}\n\n\n'.format(self.uniprot_id, self.mutations))



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
                self.logger.debug('Not evaluating mutations because no mutations specified...')
                continue
            if ((isinstance(p.d, sql_db.UniprotDomain) and not (p.t and p.t.domain_def)) or
                    (isinstance(p.d, sql_db.UniprotDomainPair) and not (p.t and p.t.domain_def_1 and p.t.domain_def_2))
                    ):
                self.logger.debug('Skipping because the template is missing domain definitions...')
                continue
            if not p.m or not p.m.model_filename:
                self.logger.debug('Skipping because no model...')
                continue
            self.logger.debug(self.mutations)
            for mutation in self.mutations.split(','):
                self.logger.debug(mutation)
                self.logger.debug('-' * 80)
                precalculated_mutations = self.db.get_uniprot_mutation(p.m, self.uniprot_id, mutation, p.d.path_to_data)
                self.logger.info(precalculated_mutations)
                if len(precalculated_mutations) == 0:
                    precalculated_mutation = None
                elif len(precalculated_mutations) == 1:
                    precalculated_mutation = precalculated_mutations[0]
                else:
                    raise Exception('Supposed to get only one precalculated mutation!')

                # Check to see if we have a precalculated mutation. Skip if all
                # parameters have been calculated; otherwise analyse the remaining
                # parameters. Create an empty mutation if the mutation has not
                # been precalculated.
                if (precalculated_mutation and
                        (precalculated_mutation.provean_score and
                        precalculated_mutation.Stability_energy_wt and
                        precalculated_mutation.ddG)):
                    p.mut.append(precalculated_mutations[0])
                    self.logger.info('Mutation has already been completely evaluated. Skipping...')
                    continue
                elif precalculated_mutation:
                    uniprot_mutation = precalculated_mutation
                # Construct empty models that will be used if we have errors
                elif isinstance(p.d, sql_db.UniprotDomain):
                    uniprot_mutation = sql_db.UniprotDomainMutation()
                    uniprot_mutation.uniprot_domain_id = p.d.uniprot_domain_id
                    uniprot_mutation.uniprot_id = self.uniprot_id
                    uniprot_mutation.mutation = mutation
                elif isinstance(p.d, sql_db.UniprotDomainPair):
                    uniprot_mutation = sql_db.UniprotDomainPairMutation()
                    uniprot_mutation.uniprot_domain_pair_id = p.d.uniprot_domain_pair_id
                    uniprot_mutation.uniprot_id = self.uniprot_id
                    uniprot_mutation.mutation = mutation

                try:
                    mut_data = self.get_mutation.get_mutation_data(p.d, p.t, p.m, self.uniprot_id, mutation)
                    uniprot_mutation = self.get_mutation.evaluate_mutation(mut_data, uniprot_mutation)
                except (errors.MutationOutsideDomain,
                        errors.MutationOutsideInterface) as e:
                    self.logger.debug('{}: {}; OK'.format(type(e), e.__str__()))
                    continue
                except possible_mutation_errors as e:
                    uniprot_mutation.mutation_errors = '{}: {}'.format(type(e), e.__str__())
                    self.logger.debug(uniprot_mutation.mutation_errors)
                self.logger.info('Adding mutation {}\n\n\n'.format(mutation))
                p.mut.append(uniprot_mutation)
                self.db.add_uniprot_mutation(uniprot_mutation, p.d.path_to_data)
        self.logger.info('Finished processing all mutations for {} {}\n\n\n'.format(self.uniprot_id, self.mutations))




    def __print_header(self, d):
        # self.logger.info('Domain or domain pair number: {}'.format(d_idx))
        if isinstance(d, sql_db.UniprotDomain):
            self.logger.debug('uniprot_domain_id: {}'.format(d.uniprot_domain_id))
        else:
            self.logger.debug('uniprot_domain_pair_id: {}'.format(d.uniprot_domain_pair_id))


    def __add_new_error(self, d_error_log, e):
        if d_error_log is None:
            return str(type(e))
        else:
            return '{}; {}'.format(d_error_log, type(e))


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


    def __prepare_provean_temp_folder(self):
        username  = hf.get_username()
        hostname = hf.get_hostname()
        if username == 'strokach' and not (('beagle' in hostname) or ('grendel' in hostname)):
            # Can't use /tmp for temporary storage on banting because I run out
            # of space and crash the node
            self.logger.debug('Using a temp folder on kimstg for provean temp files. This may lead to poor performace...')
            provean_temp_path = '/home/kimlab1/strokach/tmp/elaspic/' + self.unique
            subprocess.check_call('mkdir -p ' + self.provean_temp_path, shell=True)
        else:
            provean_temp_path = self.unique_temp_folder
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
        uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutation, run_type, n_cores)
        temp = sqlalchemy.ext.serializer.dumps(uniprot_domains_and_domain_pairs)
        with open('/tmp/elaspic/sa_pickles/{}_{}.pickle'.format(uniprot_id, mutation), 'wb') as ofh:
            pickle.dump(temp, ofh, pickle.HIGHEST_PROTOCOL)

        with open('/tmp/elaspic/sa_pickles/{}_{}.pickle'.format(uniprot_id, mutation), 'rb') as ifh:
            uniprot_domains_and_domain_pairs = pickle.load(ifh)

