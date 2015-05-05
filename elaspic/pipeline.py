# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import object

import os
import re
import time
import subprocess
import tempfile
import atexit
import signal
import six

from . import conf
from . import errors
from . import helper_functions as hf

sql_db = None
domain_alignment = None
domain_model = None
domain_mutation = None


class Pipeline(object):

    def __init__(self, configurations):
        """
        It should be possible to initialize one pipeline and call it in parallel using different
        mutations as input
        """
        global sql_db
        global domain_alignment
        global domain_model
        global domain_mutation
        
        # Read the configuration file and set the variables
        if isinstance(configurations, six.string_types):
            conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            conf.configs = configurations.copy()
        
        # Can imort sql_db only after the configuration file has been read
        # Otherwise classes in `sql_db` won't be properly configured
        from . import sql_db
        from . import domain_alignment
        from . import domain_model
        from . import domain_mutation

        # Make sure that the temp path was created successfully
        # In particular, a local temporary folder must be used on cluasters with limited local storage
        self.__validate_temp_path(conf.configs)

        # Initialize a logger
        self.logger = hf.get_logger(conf.configs['debug'])

        # Copy the sqlite version of the pdbfam database. # Deprecated! Not tested!
        if conf.configs['db_type'] == 'sqlite_file':
            self.__copy_sqlite_database()

        # Copy the blast databaes
        self.__copy_blast_database()


    def __validate_temp_path(self, configs):
        """TODO: remove so error message does not appear in a production release
        """
        hostname = hf.get_hostname()
        no_job_specific_folder = (
            configs['temp_path'].startswith(
                os.path.join(configs['global_temp_path'], configs['temp_path_suffix']))
        )
        on_node_with_manditory_job_specific_folder = (
            any([(x.lower() in hostname) for x in ['node', 'behemoth', 'grendel', 'beagle']])
        )
        if no_job_specific_folder and on_node_with_manditory_job_specific_folder:
            raise Exception('You should be using a temp folder that it specific to the particular job!')


    def __copy_sqlite_database(self):
        db_filename = conf.configs['sqlite_db_path'].strip().split('/')[-1]
        if not os.path.isfile(conf.configs['temp_path'] + db_filename):
            system_command = 'cp -u ' + conf.configs['temp_path'] + ' ' + conf.configs['sqlite_db_path'] + '/'
            self.logger.debug(system_command)
            result, error_message, return_code = hf.popen(system_command)
            if return_code != 0:
                self.logger.error(result)
                self.logger.error(error_message)
                raise Exception('Could not copy the sqlite3 database!')


    def __copy_blast_database(self):
        if (os.path.isdir(conf.configs['global_temp_path'] + 'blast/db/') and
            os.path.isfile(conf.configs['global_temp_path'] + 'blast/db/nr.pal') and
            os.path.isfile(conf.configs['global_temp_path'] + 'blast/db/pdbaa.pal')):
                system_command = (
                    'if [[ -L {global_temp_path}blast/db ]]; '
                    'then echo "True"; '
                    'else echo "False"; '
                    'fi'
                )
                system_command = system_command.format(**conf.configs)
                result, error_message, return_code = hf.popen(system_command)
                if result.strip() == 'False':
                    self.logger.warning(
                        'A symlink to the blast database already exists, but you can speed up '
                        'Provean by manually copying the blast databases to the temp folder!')
                elif result.strip() == 'True':
                    self.logger.info('A local copy of the blast database already exists. Good!')
                else:
                    raise Exception(
                        "Unexpected output '{result}' from system command '{system_command}`."
                        .format(result=result, system_command=system_command))
                return
        
        # Make a symbolic link from the remote blast database location to a local path
        self.logger.info(
            'Creating a simbolic link to the blast database folder...\n'
            'To speed up Provean, you may want to rsync blast databases to the temporary folder, '
            'rather than creating symbolic links.'
        )
        system_command = (
            'rm -rf ' + conf.configs['global_temp_path'] + 'blast/ && ' +
            'mkdir -p ' + conf.configs['global_temp_path'] + 'blast/ && ' +
            'ln -sf ' + conf.configs['blast_db_path'] + 'db ' +  conf.configs['global_temp_path'] + 'blast/'
        )
        
        # Try running the system command several times in case it doesn't work the first time
        n_tries = 0
        result, error_message, return_code = hf.popen(system_command)
        while return_code != 0 and n_tries < 10:
            self.logger.info("Couldn't copy the blast database. Retrying...")
            self.logger.debug(result)
            self.logger.debug(error_message)
            time.sleep(15)
            result, error_message, return_code = hf.popen(system_command)
        if return_code != 0:
            self.logger.error(result)
            self.logger.error(error_message)


    def __call__(self, uniprot_id, mutations, run_type=1, number_of_tries=[]):
        """ Run the main function of the program and parse errors
        """
        self.uniprot_id = uniprot_id
        self.mutations = mutations
        self.calculated_mutations = []
        self.run_type = run_type
        self.number_of_tries = number_of_tries

        self.unique = tempfile.mkdtemp(prefix='', dir=conf.configs['temp_path']).split('/')[-1]
        self.unique_temp_folder = os.path.join(conf.configs['temp_path'], self.unique) + '/'
        self.logger.info(self.unique_temp_folder)
        self.logger.info(conf.configs['db_schema'])

        # Switch to the root of the unique tmp directory
        os.chdir(self.unique_temp_folder)

        # Create subfolders for temporary storage
        self.__prepare_temp_folder()

        # Set up a temporary folder for provean
        self.provean_temp_path = self.__prepare_provean_temp_folder()
        atexit.register(self.__clear_provean_temp_files)

        ### TODO: We don't really need to provide lengthy conf input to ``MyDatabase``,
        ### ``GetModel``, and ``GetMutation``. They can read configurations from the config file
        ### themselves.

        # Initialise the sql database for accessing all information
        self.db = sql_db.MyDatabase(conf.configs, logger=self.logger)

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
            self.logger.info('\n\n\n' + '*' * 110)
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
        if conf.configs['look_for_interactions']:
            self.logger.info('Obtaining protein domain pair information...')
            self.uniprot_domain_pairs = self.db.get_uniprot_domain_pair(self.uniprot_id, True)
            self._update_path_to_data(self.uniprot_domain_pairs)

        # Make models
        if run_type in [2, 4, 5]:
            self.get_model = domain_model.GetModel(
                self.unique_temp_folder, self.db, self.logger, conf.configs)
            self.logger.info('\n\n\n' + '*' * 110)
            self.logger.info("Building models...")
            self._compute_models()

        # Analyse mutations
        if run_type in [3, 4, 5] and self.mutations:
            self.get_mutation = domain_mutation.GetMutation(
                self.unique_temp_folder, self.db, self.logger, self.provean_temp_path, conf.configs)
            self.logger.info('\n\n\n' + '*' * 110)
            self.logger.info("Analyzing mutations...")
            self._compute_mutations()


    def _update_path_to_data(self, d_list):
        if not isinstance(d_list, list):
            d_list = [d_list]
        for d in d_list:
            if not d.path_to_data or any([len(x) > 255 for x in d.path_to_data.split('/')]):
                d.path_to_data = sql_db.get_uniprot_base_path(d) + sql_db.get_uniprot_domain_path(d)
                self.db.merge_row(d)
            subprocess.check_call('mkdir -p {}'.format(conf.configs['temp_path'] + d.path_to_data), shell=True)


    def _compute_provean(self):
        """
        """
        d = self.uniprot_domains[0] # All uniprot domains will have the same uniprot_sequence info
        self.__print_header(d)

        if (d.uniprot_sequence.provean and
            d.uniprot_sequence.provean.provean_supset_filename):
                path_to_provean_supset = (
                    conf.configs['temp_path'] + sql_db.get_uniprot_base_path(d) +
                    d.uniprot_sequence.provean.provean_supset_filename)
#                path_to_provean_supset = (
#                    conf.configs['path_to_archive'] + sql_db.get_uniprot_base_path(d) +
#                    d.uniprot_sequence.provean.provean_supset_filename)
                if os.path.isfile(path_to_provean_supset):
                    if not conf.configs['remake_provean_supset']:
                        self.logger.debug('The provean supset has already been calculated. Done!\n')
                        return None
                    first_aa = d.uniprot_sequence.uniprot_sequence[0]
                    domain_mutation = '{0}1{0}'.format(first_aa)
                    result, error_message, return_code = \
                        domain_alignment.check_provean_supporting_set(
                            domain_mutation, d.uniprot_sequence.uniprot_sequence,
                            conf.configs, self.unique_temp_folder, self.provean_temp_path,
                            self.logger, d.uniprot_id, path_to_provean_supset,
                            save_supporting_set=False, check_mem_usage=False)
                    if return_code == 0:
                        self.logger.debug('The provean supset has already been calculated. Done!\n')
                        return None
                    self.logger.debug('Provean supporting set caused an error:\n\n{}'.format(error_message))
                    self.logger.debug('Recompiling...')
                    provean = d.uniprot_sequence.provean
                else:
                    self.logger.error(
                        'Provean has been calculated but the file is missing from:\n{}\nRecompiling...'
                        .format(path_to_provean_supset))
                    provean = d.uniprot_sequence.provean
        elif d.uniprot_sequence.provean:
            self.logger.info('Provean supporting set has not been calculated previously. Computing...')
            provean = d.uniprot_sequence.provean
        else:
            provean = sql_db.Provean()
            provean.uniprot_id = d.uniprot_id

        try:
            provean.provean_supset_filename, provean.provean_supset_length = \
                domain_alignment.build_provean_supporting_set(
                    d.uniprot_id,
                    d.uniprot_sequence.uniprot_name,
                    d.uniprot_sequence.uniprot_sequence,
                    conf.configs, self.unique_temp_folder, self.provean_temp_path, self.logger)

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
        self.db.merge_provean(provean, sql_db.get_uniprot_base_path(d))
        self.logger.info('Finished computing provean for {}\n'.format(d.uniprot_id))
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
                    (d.template.model.model_errors.count(';') > 10)) ):
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
            elif ( (isinstance(d, sql_db.UniprotDomain) and not d.template.domain) or
                    (isinstance(d, sql_db.UniprotDomainPair) and not (d.template.domain_1 and d.template.domain_2)) ):
                self.logger.debug('Skipping because no structural template is availible...')
                continue
            elif d.template.model == None or d.template.model.model_filename == None:
                self.logger.debug('Skipping because no model...')
                continue
            elif d.template.model.model_errors != None:
                self.logger.debug('Skipping because the model has errors: {}!'.format(d.template.model.model_errors))
                continue

            self.logger.debug('Going over all mutations: {}'.format(self.mutations))
            for mutation in self.mutations.split(','):
                self.logger.debug('-' * 80)
                self.logger.debug('Mutation: {}'.format(mutation))
                mutation_prototype = re.compile("^[A-z][0-9]+[A-z]$")
                if not mutation_prototype.match(mutation) or int(mutation[1:-1]) == 0:
                    self.logger.error('The mutation {} is not a supported format! Skiping!'.format(mutation))
                    continue
                # Check to see if we have a precalculated mutation. Skip if all
                # parameters have been calculated; otherwise analyse the remaining
                # parameters. Create an empty mutation if the mutation has not
                # been precalculated.
                precalculated_mutation = self.db.get_uniprot_mutation(d, mutation, self.uniprot_id, True)
                self.logger.info('Have the following precalculated mutation: {}'.format(precalculated_mutation))
                if (precalculated_mutation and
                        (precalculated_mutation.provean_score and
                        precalculated_mutation.stability_energy_wt and
                        precalculated_mutation.ddg != None and
                        precalculated_mutation.ddg != 1.0)): # TODO: Remove this line
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
                except errors.FoldXAAMismatchError as e:
                    self.logger.error('{}: {}'.format(type(e), e))
                    unique_id_string, unique_id = self.__get_unique_id(d)
                    if self.number_of_tries.count(unique_id) > 3:
                        self.logger.error(
                            'An error occured more than three times for {}: {}! It cannot be fixed!'
                            .format(unique_id_string, unique_id)
                        )
                        raise e
                    self.number_of_tries.append(unique_id)
                    self.logger.debug(
                        'Deleting the model for {}: {} and trying to run the pipeline again... number_of_tries: {}'
                        .format(unique_id_string, unique_id, self.number_of_tries))
                    self.db.remove_model(d)
                    self.__call__(
                        self.uniprot_id, self.mutations, self.run_type, conf.configs['n_cores'],
                        self.number_of_tries, self.logger
                    )
                except (errors.MutationOutsideDomainError,
                        errors.MutationOutsideInterfaceError) as e:
                    self.logger.debug('{}: {}; OK'.format(type(e), e))
                    continue
                except possible_mutation_errors as e:
                    uniprot_mutation.mutation_errors = '{}: {}'.format(type(e), e)
                    self.logger.debug(uniprot_mutation.mutation_errors)

                self.logger.info('Adding mutation {}'.format(mutation))
                self.logger.debug("Mutation attributes:")
                for attr in dir(uniprot_mutation):
                    if attr.startswith('_'):
                        continue
                    attr_field = getattr(uniprot_mutation, attr)
                    attr_type = type(attr_field)
                    self.logger.debug(attr)
                    self.logger.debug(attr_field)
                    self.logger.debug(attr_type)
                    if (six.PY2 and
                        (isinstance(attr_field, six.binary_type) or
                         isinstance(attr_field, six.text_type))):
                            self.logger.debug(
                                 'Changing attribute {} from {} to {}...'
                                 .format(attr, attr_type, str(attr_field)))
                            setattr(uniprot_mutation, attr, str(attr_field))
                    if six.PY3 and isinstance(attr_field, six.binary_type):
                        self.logger.debug(
                            'Changing attribute {} from {} to {}...'
                            .format(attr, attr_type, type(attr_field.decode())))
                        setattr(uniprot_mutation, attr, attr_field.decode())
                self.calculated_mutations.append(uniprot_mutation)
                self.db.merge_mutation(uniprot_mutation, d.path_to_data)
                self.uniprot_mutations.append(uniprot_mutation)
        self.logger.info('Finished processing all mutations for {} {}'.format(self.uniprot_id, self.mutations))


    def __get_unique_id(self, d):
        if isinstance(d, sql_db.UniprotDomain):
            return ('uniprot_domain_id', d.uniprot_domain_id)
        else:
            return('uniprot_domain_pair_id', d.uniprot_domain_pair_id)


    def __print_header(self, d):
        # self.logger.info('Domain or domain pair number: {}'.format(d_idx))
        self.logger.info('=' * 77)
        unique_id = self.__get_unique_id(d)
        self.logger.info('{}: {}'.format(*unique_id))


    def __add_new_error(self, d_error_log, e):
        if d_error_log is None:
            return str(type(e))
        else:
            return '{}; {}: {}'.format(d_error_log, type(e), str(e))


    def __prepare_temp_folder(self):
        # base temp directory
        if not os.path.isdir(conf.configs['temp_path']):
            subprocess.check_call('mkdir -p ' + conf.configs['temp_path'], shell=True)

        # unique_temp_folder
        if not os.path.isdir(self.unique_temp_folder):
            subprocess.check_call('mkdir -p ' + self.unique_temp_folder, shell=True)

        # t_coffee
        if not os.path.isdir(self.unique_temp_folder + '/tcoffee'):
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/tcoffee && ' + \
                            'mkdir -p ' + self.unique_temp_folder + '/tcoffee/tmp && ' + \
                            'mkdir -p ' + self.unique_temp_folder + '/tcoffee/lck && ' + \
                            'mkdir -p ' + self.unique_temp_folder + '/tcoffee/cache'
            subprocess.check_call(mkdir_command, shell=True)

        # folx
        if not os.path.isdir(self.unique_temp_folder + '/FoldX'):
            subprocess.check_call("mkdir -p '{}'".format(self.unique_temp_folder + '/FoldX'), shell=True)

        # modeller
        if not os.path.isdir(self.unique_temp_folder + '/modeller'):
            # create workingfolder for modeller
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/modeller'
            subprocess.check_call(mkdir_command, shell=True)

        # sequence conservation
        if not os.path.isdir(self.unique_temp_folder + '/sequence_conservation'):
            mkdir_command = 'mkdir -p ' + self.unique_temp_folder + '/sequence_conservation'
            subprocess.check_call(mkdir_command, shell=True)

        # analyze_structure
        analyze_structure_path = os.path.join(self.unique_temp_folder, 'analyze_structure') + '/'
        if not os.path.isdir(analyze_structure_path):
            mkdir_command = 'mkdir -p ' + analyze_structure_path
            subprocess.check_call(mkdir_command, shell=True)


    def __prepare_provean_temp_folder(self):
        """
        Some nodes on the cluster have a very limited amount of memory for temp storage.
        When working on those nodes, you sould use a remote location for temp storage. This is slow,
        but at least it ensures that you don't crush the nodes by filling up the hard drive.
        However, the most serious problem should be fixed with an up-to-date version of cd-hit.
        (Older versions could go into an infinite loop and generate huge temp files).
        """
        hostname = hf.get_hostname()
        if conf.configs['provean_temp_path'] and not (('beagle' in hostname) or ('grendel' in hostname)):
            # Can't use /tmp for temporary storage on banting because I run out
            # of space and crash the node
            self.logger.debug(
                "Using a temp folder '{provean_temp_path}' for provean temp files. "
                "This may lead to poor performace...".format(**conf.configs))
            provean_temp_path = os.path.join(conf.configs['provean_temp_path'], self.unique) + '/'
        else:
            provean_temp_path = self.unique_temp_folder + 'provean_temp/'
        subprocess.check_call('mkdir -p ' + provean_temp_path, shell=True)
        return provean_temp_path


    def __clear_provean_temp_files(self):
        self.logger.debug("Clearning provean temporary files from '{}'".format(self.provean_temp_path))
        try:
            subprocess.check_call('rm -rf ' + self.provean_temp_path + 'provean*', shell=True)
            subprocess.check_call('rm -rf ' + self.provean_temp_path + '*cdhit*', shell=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(
                'Removing provean temporary files failed with the following error: {}\n'
                'This is probably because the files are still being used...'
                .format(e))

