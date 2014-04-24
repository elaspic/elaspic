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

from modeller import ModellerError

import errors
import sql_db

import domain_template
import domain_model
import domain_mutation

from Bio.SubsMat import MatrixInfo
from ConfigParser import SafeConfigParser
from helper_functions import scinetCleanup




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
        # read the configuration file and set the variables
        configParser = SafeConfigParser(
            defaults={
                'global_temp_path': '/tmp/',
                'tmpPath': '/tmp/pipeline/',
                'HETATM': True,
                'DEBUG': False,
                'saveTo': '$SCRATCH/niklas-pipeline/',
                'saveScinet': False,
                'path_to_archive': '/home/kimlab1/database_data/elaspic/human/',
                'db_type': sql_db.sql_flavor,
                'db_path': errors.path_to_pipeline_code() + '/../db/pipeline.db',
                'look_for_interactions': 'True',
                'webServer': False,
                'numConsumers': '1',
                'tcoffee_parallel_runs': '1',
                'outputPath': '../results/',
                'savePDB': 'results/pdb_files/',
                'runTime': 'INFINITE',
                'matrix': 'blosum80',
                'gap_start': '-16',
                'gap_extend': '-4',
                'n_cores': '1',
                'run_provean': False})

        configParser.read(configFile)

        # from [DEFAULT]
        self.global_temp_path = configParser.get('DEFAULT', 'global_temp_path')
        tmpPath = configParser.get('DEFAULT', 'tmpPath').strip('/') + '/'
        self.DEBUG = configParser.getboolean('DEFAULT', 'DEBUG')
        self.HETATM = configParser.getboolean('DEFAULT', 'HETATM')
        self.saveTo = configParser.get('DEFAULT', 'saveTo')
        self.saveScinet = configParser.getboolean('DEFAULT', 'saveScinet')
        self.path_to_archive = configParser.get('DEFAULT', 'path_to_archive')
        self.db_type = configParser.get('DEFAULT', 'db_type')
        self.db_path = configParser.get('DEFAULT', 'db_path')
        self.look_for_interactions = configParser.getboolean('DEFAULT', 'look_for_interactions')
        self.webServer = configParser.get('DEFAULT', 'webServer')
        self.n_cores = configParser.get('DEFAULT', 'n_cores')
        self.run_provean = configParser.getboolean('DEFAULT', 'run_provean')

        # from [SETTINGS]
        self.num_consumers = configParser.getint('SETTINGS', 'numConsumers')
        self.tcoffee_parallel_runs = configParser.getint('SETTINGS', 'tcoffee_parallel_runs')
        self.outputPath = configParser.get('SETTINGS', 'outputPath')
        self.savePDB = configParser.get('SETTINGS', 'savePDB')
        self.runTime = configParser.get('SETTINGS', 'runTime')
        self.pdbPath = configParser.get('SETTINGS', 'pdbPath')
        self.executables = configParser.get('SETTINGS', 'bin')

        # from [INPUT]
        self.template_finding = configParser.getboolean('INPUT', 'mutation_uniprot')

        # from [MODELLER]
        self.modeller_runs = configParser.getint('MODELLER', 'modeller_runs')

        # from [FOLDX]
        self.foldX_WATER = configParser.get('FOLDX', 'WATER')
        self.buildModel_runs = configParser.get('FOLDX', 'buildModel_runs')


        #######################################################################
        # matrix (currently hardcoded, can be changed if needed)
        # in that case also change the gap_start and gap_extend options
        # they are used in conjuntion with the matrix to determine the
        # sequence similariy (could be used to determine the sequence similarity
        # between to interfaces. Might help to improve the modelling)

        # matrix type
        self.matrix_option = configParser.get('SETTINGS', 'matrix')
#        self.matrix_option = 'blosum80'

        # gap_start
        self.gap_start = configParser.getint('SETTINGS', 'gap_start')
#        self.gap_start = -16

        # gap_extend
        self.gap_extend = configParser.getint('SETTINGS', 'gap_extend')
#        self.gap_extend = -4

        # set the matrix for the substitution score
        if self.matrix_option == 'blosum80':
            self.matrix =  MatrixInfo.blosum80
        elif self.matrix_option == 'blosum60':
            self.matrix =  MatrixInfo.blosum60
        else:
            raise Exception('Specified matrix not found!')


        #######################################################################
        # check the TMPDIR
        # if a TMPDIR is given as environment variable the tmp directory
        # is created relative to that. This is useful when running on banting
        # (the cluster in the ccbr) and also on Scinet (I might have set the
        # environment variable on Scinet myself though..). Make sure that it
        # points to '/dev/shm/' on Scinet

        childProcess = subprocess.Popen('echo $TMPDIR', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, __ = childProcess.communicate()
        TMPDIR_CLUSTER = result.strip()
        try:
            if TMPDIR_CLUSTER[-1] == '/':
                # the tmpPath might be given as an absolute Path
                # thus the last '/' has to be removed
                self.TMPDIR_CLUSTER = TMPDIR_CLUSTER[:-1]
            else:
                self.TMPDIR_CLUSTER = TMPDIR_CLUSTER
        except IndexError:
            self.TMPDIR_CLUSTER = TMPDIR_CLUSTER

        if self.TMPDIR_CLUSTER:
            if tmpPath[0] == '/':
                # i.e. tmpPath is given as an absolute Path
                self.tmpPath = self.TMPDIR_CLUSTER + tmpPath
            else:
                self.tmpPath = self.TMPDIR_CLUSTER + '/' + tmpPath
        else:
            if tmpPath[0] == '/':
                # i.e. tmpPath is given as an absolute Path
                self.tmpPath = self.global_temp_path + tmpPath
            else:
                self.tmpPath = self.global_temp_path + '/' + tmpPath

        subprocess.check_call('mkdir -p ' + self.tmpPath, shell=True)

        #######################################################################
        # if running on the cluster copy the database to the tmp DIR for
        # speedup and to avoid killing the network. BLAST is very IO intensive
        # and you don't want that to be run over the network!
        #
        # I can distinguish where the pipeline is running by checking the username
        # you will have to adjust that!
        # my usernames are:
        # local: niklas
        # banting: nberliner
        # Scinet: joan
        childProcess = subprocess.Popen('whoami', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        whoami, __ = childProcess.communicate()
        path_to_local_blast_db = '/home/kimlab1/database_data/blast/'
#        if not os.path.isdir(self.tmpPath + 'blast/'):
##            childProcess = subprocess.Popen('hostname | cut -d. -f1', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
##            host_name, __ = childProcess.communicate()
#            if whoami.strip() == 'strokach':
#                # when running the pipeline on beagle or banting, copy the database
#                # to a local folder
#                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
#                                    'cp -r ' + path_to_local_blast_db + 'pdbaa_db/ ' + self.tmpPath + 'blast/pdbaa_db/ && ' + \
#                                    'cd ' + self.tmpPath + 'blast && ' + \
#                                    'ln -sf ' + path_to_local_blast_db + 'db'
#            elif whoami.strip() == 'alexey':
#                # when running the pipeline locally there is no need to copy the database
#                # a symlink is enough
#                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
#                                    'cd ' + self.tmpPath + 'blast && ' + \
#                                    'ln -sf ' + path_to_local_blast_db + 'pdbaa_db && ' + \
#                                    'ln -sf ' + path_to_local_blast_db + 'db'
#            elif whoami.strip() == 'witvliet':
#                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
#                                    'cd ' + self.tmpPath + 'blast && ' + \
#                                    'ln -sf /home/witvliet/working/bin/ncbi-blast-2.2.28+/pdbaa_db'
#            elif whoami.strip() == 'joan':
#                # for scinet, blast is already installed, but the database needs to be copied
#                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
#                                    'cp -ru $HOME/niklas-pipeline/blastdb/pdbaa_db ' + \
#                                    self.tmpPath + 'blast/'


        childProcess = subprocess.Popen('hostname | cut -d. -f1', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        host_name, __ = childProcess.communicate()
        if (whoami.strip() == 'strokach' and ( ('beagle' in host_name) or ('grendel' in host_name) )):
            system_command = (
                'mkdir -p ' + self.global_temp_path + 'blast/pdbaa_db/ && ' +
                'mkdir -p ' + self.global_temp_path + 'blast/db/ && ' +
                'rsync -rzu ' + path_to_local_blast_db + 'pdbaa_db/ ' + self.global_temp_path + 'blast/pdbaa_db/ && ' +
                'rsync -rzu ' + path_to_local_blast_db + 'db/ ' + self.global_temp_path + 'blast/db/')
        elif whoami.strip() == 'strokach' or whoami.strip() == 'alexey':
            system_command = (
                'mkdir -p ' + self.global_temp_path + 'blast/ && ' +
                'cd ' + self.global_temp_path + 'blast/ && ' + \
                'ln -sf ' + path_to_local_blast_db + 'pdbaa_db/ ' +
                'ln -sf ' + path_to_local_blast_db + 'db/ ')
        elif whoami.strip() == 'witvliet':
            system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                'cd ' + self.tmpPath + 'blast && ' + \
                                'ln -sf /home/witvliet/working/bin/ncbi-blast-2.2.28+/pdbaa_db'
        elif whoami.strip() == 'joan':
            # for scinet, blast is already installed, but the database needs to be copied
            system_command = 'mkdir -p ' + self.global_temp_path + 'blast && ' + \
                                'cp -ru $HOME/niklas-pipeline/blastdb/pdbaa_db ' + \
                                self.global_temp_path + 'blast/'
        print system_command
        childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, error_message = childProcess.communicate()
        if childProcess.returncode != 0:
            print result
            print error_message
            raise Exception('Couldn\'t copy the blast database!')


        #######################################################################
        # Copy the database file
        if self.db_type == 'sqlite_file' and not os.path.isfile(self.tmpPath + '/' + self.db_path.strip().split('/')[-1]):
            print self.tmpPath
            print self.db_path
            system_command = 'cp -u ' + self.db_path + ' ' + self.tmpPath + '/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, error_message = childProcess.communicate()
            if childProcess.returncode != 0:
                print result
                print error_message
                raise Exception('Could not copy the sqlite3 database!')


        #######################################################################
        # Create output folder if it doesn't exist
#        if not os.path.isdir(self.outputPath):
#            subprocess.check_call('mkdir -p ' + self.outputPath, shell=True)

        #######################################################################



    def __call__(self, uniprot_id, mutation):
        """ Run the main function of the program and parse errors
        """

        self.uniprot_id = uniprot_id
        self.mutations = mutation
        print self.uniprot_id

        self.unique = tempfile.mkdtemp(prefix='', dir=self.tmpPath).split('/')[-1]
        print self.unique

        # create temporary folders
        self.__prepare_temp_folder()

        # initialize the logger
        logger = logging.getLogger(__name__)
        if self.DEBUG:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
        if False:
            handler = logging.FileHandler(self.outputPath + self.uniprot_id + '_' +
                            self.mutation + '.log', mode='w', delay=True, maxsize=0, rotate=5)
        else:
            handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.handlers = []
        logger.addHandler(handler)
        self.log = logger

        # Initialise the sql database for accessing all information
        self.db = sql_db.MyDatabase(path_to_sqlite_db=(self.tmpPath + 'pipeline.db'),
                                 sql_flavor=self.db_type,
                                 is_immutable=False,
                                 path_to_temp=self.tmpPath,
                                 path_to_archive=self.path_to_archive)

        # go to a unique directory
        self.PWD = self.tmpPath + self.unique + '/'
        os.chdir(self.PWD)

        # Initialize the classes used for calculating templates, models and mutations
        self.get_template = domain_template.GetTemplate(
            self.global_temp_path, self.tmpPath, self.unique, self.pdbPath, self.db, self.log, self.n_cores)

        self.get_model = domain_model.GetModel(
            self.global_temp_path, self.tmpPath, self.unique, self.pdbPath, self.db, self.log,
            self.modeller_runs, self.buildModel_runs, self.PWD, self.n_cores)

        self.get_mutation = domain_mutation.GetMutation(
            self.global_temp_path, self.tmpPath, self.unique, self.pdbPath, self.db, self.log,
            self.foldX_WATER, self.buildModel_runs, self.PWD, self.matrix, self.gap_start, self.gap_extend, self.n_cores)


        #######################################################################
        # Find all domains and domain pairs for a given uniprot
        protein_domains = self.db.get_uniprot_domain(self.uniprot_id)
        for d, t, m in protein_domains:
            if not d.path_to_data:
                d.path_to_data = (d.uniprot_name.split('_')[-1].lower() + '/' +
                    d.uniprot_id[:3] + '/' + d.uniprot_id[3:5] + '/' + d.uniprot_id + '/' +
                    d.pfam_name + '*' + d.envelope_def.replace(':','-') + '/')
        if self.look_for_interactions:
            protein_domain_pairs = self.db.get_uniprot_domain_pair(self.uniprot_id)
            for d, t, m in protein_domain_pairs:
                if not d.path_to_data:
                    d.path_to_data = (d.uniprot_domain_1.uniprot_name.split('_')[-1].lower() + '/' +
                    d.uniprot_domain_1.uniprot_id[:3] + '/' + d.uniprot_domain_1.uniprot_id[3:5] + '/' + d.uniprot_domain_1.uniprot_id + '/' +
                    d.uniprot_domain_1.pfam_name + '*' + d.uniprot_domain_1.envelope_def.replace(':','-') + '/' +
                    d.uniprot_domain_2.pfam_name + '*' + d.uniprot_domain_2.envelope_def.replace(':','-') + '/' +
                    d.uniprot_domain_2.uniprot_id + '/')
            protein_domain_template_model = protein_domains + protein_domain_pairs
        else:
            protein_domain_template_model = protein_domains
        if protein_domain_template_model == []:
            self.log.error('Uniprot %s has no pfam domains' % self.uniprot_id)
        # Convert to a single object which holds domain, template, model, and mutations
        # (so that I can modify each of those parameters "in place" in a for loop)
        self.protein_definitions = []
        for d, t, m in protein_domain_template_model:
            self.protein_definitions.append(ProteinDomainData(d, t, m))


        #######################################################################
        # Do the calculations
        if (not self.run_provean):
            self.log.info("Finding templates...")
            self._compute_templates()

            self.log.info("Building models...")
            self._compute_models()

            self.log.info("Analyzing mutations...")
            self._compute_mutations()

        else:
            self.log.info("Compute provean...")
            self._compute_provean()
        # Need to add some stuff for foldX: <analyse complex>


        #######################################################################
        # save the results from ramdisk
        if self.saveScinet:
            scinetCleanup(self.outputPath, self.saveTo)

        # Leave everything the way you found it
        self.db.session.close()

        os.chdir(self.PWD)

        return self.protein_definitions

        # modeller_path = self.tmpPath + self.unique + '/modeller/'

#        if self.template_finding:

#        #######################################################################
#
#        if self.template_finding and self.mutations != '':
#            # If we are given a mutation (or a list of tab-separated mutations)
#            # calculate the dG for each of those mutations
#
#
#        #######################################################################
#
#        if not self.template_finding and self.mutations != '':
#
#            # Modell mutations using the structure from the pdb
#            template, template_mutations = self.get_pdb_and_mutations(self.mutations)
#
#            # Get the data for a particular mutation (or a set of mutations)
#            template_mutations = self.analyse_mutations(template, template_mutations)
#
#            protein_definitions = [template.update(template_mutations),]
#
#
#        # Switch back to the original directory
#        os.chdir(self.PWD)
#
#        return protein_definitions


    def _compute_templates(self):
        """
        Align uniprot domains to structural domains in pdbfam, find the best
        template, and expand domain boundaries
        """

        # Go over all domains and domain pairs for a given protein
        for p in self.protein_definitions:
            self.log.info('%s\t%s\t%s' % (p.d, p.t, p.m,) )

            if p.t:
                if False:
                    # So that I can comment out all the elifs
                    pass
#                elif d.model and d.model.model_errors:
#                    # Recalculate templates if the model has errors
#                    self.log.debug('Recalculating templates because the model had errors')
#                    pass
                elif isinstance(
                p.t, sql_db.UniprotDomainPairTemplate) and \
                p.t.domain_1 and p.t.domain_2 and \
                p.t.domain_1 == p.t.domain_2 and \
                p.t.domain_1.pdb_chain == p.t.domain_2.pdb_chain:
                    self.log.debug('Recalculating templates because both interacting partners mapped to the same pdb chain')
                    pass
#                elif not d.template.template_errors == 'no templates found':
#                    # Recalculate cases where a template was found but errors occured
#                    self.log.debug('Recalculating templates because odd errors occured the last time we tried')
#                    pass
#                elif isinstance(p.t, sql_db.UniprotDomainTemplate) \
#                and not p.t.provean_supset_filename:
#                    self.log.debug('Don\'t have a supporting sequence alignment for provean')
#                    pass
                else:
                    self.log.info('skipping template')
                    continue

            # Construct empty templates that will be used if we have errors
            if type(p.d) == sql_db.UniprotDomain:
                empty_template = sql_db.UniprotDomainTemplate()
                empty_template.uniprot_domain_id = p.d.uniprot_domain_id
            elif type(p.d) == sql_db.UniprotDomainPair:
                empty_template = sql_db.UniprotDomainPairTemplate()
                empty_template.uniprot_domain_pair_id = p.d.uniprot_domain_pair_id

            # Find templates and catch errors
            try:
                if isinstance(p.t, sql_db.UniprotDomainTemplate) \
                and p.t.alignment_filename \
                and not p.t.template_errors \
                and not p.t.provean_supset_filename:
                    template = p.t
                    (template.provean_supset_filename, template.provean_supset_length) \
                        = self.get_template.build_provean_supporting_set(p.d, template)
                else:
                    template = self.get_template(p.d)
            except (
            errors.NoStructuralTemplates,
            errors.NoSequenceFound,
            errors.NoTemplatesFound,
            errors.TcoffeeError,
            errors.ProveanError) as e:
                self.log.error(e.message)
                empty_template.template_errors = e.message
                template = empty_template

            self.log.info('adding template')
            p.t = template
            self.db.add_uniprot_template(template, p.d.path_to_data)

        self.log.info('Finished processing all templates for ' + self.uniprot_id + ' ' + self.mutations + '\n')



    def _compute_provean(self):
        """
        Align uniprot domains to structural domains in pdbfam, find the best
        template, and expand domain boundaries
        """

        # Go over all domains and domain pairs for a given protein
        for p in self.protein_definitions:
            self.log.info( '%s\t%s\t%s' % (p.d, p.t, p.m,) )

            # Make a provean supset for the template
            if (
            not p.t
            or (isinstance(p.t, sql_db.UniprotDomainTemplate) and not p.t.alignment_filename)
            or (isinstance(p.t, sql_db.UniprotDomainPairTemplate) and not (p.t.alignment_filename_1 and p.t.alignment_filename_2))):
                self.log.info('No template availible. Skipping...')
                continue

            if (isinstance(p.t, sql_db.UniprotDomainTemplate)
            and not p.t.provean_supset_filename
            and not p.t.template_errors):
                try:
                    p.t.provean_supset_filename, p.t.provean_supset_length \
                        = self.get_template.build_provean_supporting_set(p.d, p.t)
                except errors.ProveanError as e:
                    self.log.error(e.message)
                    if p.t.template_errors:
                        p.t.template_errors += ', ' + e.message
                    else:
                        p.t.template_errors = e.message
                    self.log.error('Skipping...')
                    continue
                finally:
                    rc = subprocess.check_call('rm -rf ' + self.tmpPath + self.unique + '/provean*', shell=True)
                    rc = subprocess.check_call('rm -rf ' + self.tmpPath + self.unique + '/*cdhit*', shell=True)
                    
                self.log.info('Adding template with provean info')
                self.db.add_uniprot_template(p.t, p.d.path_to_data)

            if not p.m:
                self.log.error('No model availible. Skipping...')
                continue

            # Get a provean score for each mutation
            for mutation in self.mutations.split(','):
                precalculated_mutations = self.db.get_uniprot_mutation(p.m, self.uniprot_id, mutation, p.d.path_to_data)
                self.log.info(precalculated_mutations)
                if precalculated_mutations:
                    precalculated_mutation = precalculated_mutations[0]
                    if precalculated_mutation.provean_score:
                        p.mut.append(precalculated_mutation)
                        self.log.info('skipping mutation')
                        self.log.info(precalculated_mutation)
                        continue
                    try:
                        mut_data = self.get_mutation.get_mutation_data(p.d, p.t, p.m, self.uniprot_id, mutation)
                    except errors.MutationOutsideDomain as e:
                        self.log.error(e.message)
                        continue
                    except errors.MutationOutsideInterface as e:
                        self.log.error(e.message)
                        continue
                    provean_mutation, provean_score = None, None
                    if mut_data.path_to_provean_supset:
                        provean_mutation, provean_score = \
                        self.get_mutation.get_provean_score(
                            mut_data.uniprot_domain_id,
                            mut_data.mutation_domain,
                            mut_data.domain_sequences[0],
                            mut_data.path_to_provean_supset)
                    precalculated_mutation.provean_score = provean_score
                    self.log.debug('provean mutation:')
                    self.log.debug(provean_mutation)
                    self.log.debug('provean score:')
                    self.log.debug(provean_score)

                    self.log.info('adding mutation %s' % mutation)
                    p.mut.append(precalculated_mutation)
                    self.db.add_uniprot_mutation(precalculated_mutation, p.d.path_to_data)

        self.log.info(
            'Finished adding provean info to all templates and mutations in: -u %s -m %s'
            % (self.uniprot_id, self.mutations,))



    def _compute_models(self):
        """
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """

        # Go over all domains and domain pairs for a given protein
        for p in self.protein_definitions:
            self.log.info('%s\t%s\t%s' % (p.d, p.t, p.m,) )

            # Check if we have models from a previous run
            if p.m and not p.m.model_errors:
                self.log.info('skipping model')
                continue

            # Construct empty models that will be used if we have errors
            if isinstance(p.d, sql_db.UniprotDomain):
                if not p.t.domain_def:
                    self.log.error('Domain definitions missing. Skipping...')
                    continue
                empty_model = sql_db.UniprotDomainModel()
                empty_model.uniprot_domain_id = p.t.uniprot_domain_id
            elif isinstance(p.d, sql_db.UniprotDomainPair):
                if not (p.t.domain_def_1 and p.t.domain_def_2):
                    self.log.error('Domain definitions missing. Skipping...')
                    continue
                empty_model = sql_db.UniprotDomainPairModel()
                empty_model.uniprot_domain_pair_id = p.t.uniprot_domain_pair_id

            # Make a homology model using the domain and template information
            try:
                uniprot_model = self.get_model(p.d, p.t)
            except (
            errors.PDBChainError,
            errors.ModellerFailure,
            errors.ChainsNotInteracting) as e:
                self.log.error(e.message)
                self.log.error('Skipping...')
                continue
#                empty_model.model_errors = e.message
#                uniprot_model = empty_model
            except ModellerError as e:
                self.log.error('ModellerError while trying to modellel in __getPDB for ' + self.uniprot_id + ':' + self.mutations)
                self.log.error('ModellerError args:' + '\n'.join(e.args))
                self.log.error('ModellerError message:' + e.message)
                if 'Alignment sequence not found in PDB file' in e.message:
                    empty_model.model_errors = 'alignment sequence not found in PDB file'
                else:
                    empty_model.model_errors = 'modelling error with message: %s; %s' % (e.message, ','.join(e.args),)
                uniprot_model = empty_model

            self.log.info('adding model')
            p.m = uniprot_model
            self.db.add_uniprot_model(uniprot_model, p.d.path_to_data)

        self.log.info('Finished processing all models for ' + self.uniprot_id + ' ' + self.mutations + '\n')



    def _compute_mutations(self):
        """
        """

        for p in self.protein_definitions:
            self.log.info('%s\t%s\t%s' % (p.d, p.t, p.m) )
            if not p.m \
            or not p.m.model_filename:
                continue

            for mutation in self.mutations.split(','):
                precalculated_mutation = self.db.get_uniprot_mutation(p.m, self.uniprot_id, mutation, p.d.path_to_data)
                self.log.info(precalculated_mutation)
                if precalculated_mutation:
                    p.mut.append(precalculated_mutation)
                    self.log.info('skipping mutation')
                    continue

                # Construct empty models that will be used if we have errors
                if type(p.d) == sql_db.UniprotDomain:
                    if not p.t.domain_def:
                        continue
                    empty_mutation = sql_db.UniprotDomainMutation()
                    empty_mutation.uniprot_domain_id = p.d.uniprot_domain_id
                elif type(p.d) == sql_db.UniprotDomainPair:
                    if not p.t.domain_def_1 or not p.t.domain_def_2:
                        continue
                    empty_mutation = sql_db.UniprotDomainPairMutation()
                    empty_mutation.uniprot_domain_pair_id = p.d.uniprot_domain_pair_id
                empty_mutation.uniprot_id = self.uniprot_id
                empty_mutation.mutation = mutation

                try:
                    mut_data = self.get_mutation.get_mutation_data(p.d, p.t, p.m, self.uniprot_id, mutation)
                    uniprot_mutation = self.get_mutation.evaluate_mutation(self, mut_data)

                except errors.FoldXError as e:
                    self.log.error('FoldXError while repairing the wildtype for ' + self.uniprot_id + ':' + self.mutations)
                    self.log.error('FoldXError error:' + e.message)
                    empty_mutation.mutation_errors = 'FoldXError error:' + e.message
                    uniprot_mutation = empty_mutation
                except errors.pdbError as e:
                    self.log.error(e.message)
                    empty_mutation.mutation_errors = 'pdb error:' + e.message
                    uniprot_mutation = empty_mutation
                except errors.MutationOutsideDomain as e:
                    self.log.error(e.message)
                    empty_mutation.mutation_errors = e.message
                    uniprot_mutation = empty_mutation
#                    continue
                except errors.MutationOutsideInterface as e:
                    self.log.error(e.message)
#                    empty_mutation.mutation_errors = e.message
#                    uniprot_mutation = empty_mutation
                    continue
                self.log.info('adding mutation %s' % mutation)
                p.mut.append(uniprot_mutation)
                self.db.add_uniprot_mutation(uniprot_mutation, p.d.path_to_data)

        self.log.info('Finished processing all mutations for ' + self.uniprot_id + ' ' + self.mutations + '\n')



    def __prepare_temp_folder(self):
        # create the basic tmp directory
        # delete its content if it exists (AS: disabled so that I can continue from previous run)
        if not os.path.isdir(self.tmpPath):
            subprocess.check_call('mkdir -p ' + self.tmpPath, shell=True)
#        else:
#            if not self.tmpPath[-1] == '/':
#                self.tmpPath = self.tmpPath +'/'
#            subprocess.check_call('rm -r ' + self.tmpPath + '*', shell=True)

        for i in range(1, self.num_consumers + 1):
            # the consumers
            if not os.path.isdir(self.tmpPath + self.unique):
                subprocess.check_call('mkdir -p ' + self.tmpPath + self.unique, shell=True)

            # tcoffee
            if not os.path.isdir(self.tmpPath + self.unique + '/tcoffee'):
                # create tmp for tcoffee
                mkdir_command = 'mkdir -p ' + self.tmpPath + self.unique + '/tcoffee && ' + \
                                'mkdir -p ' + self.tmpPath + self.unique + '/tcoffee/tmp && ' + \
                                'mkdir -p ' + self.tmpPath + self.unique + '/tcoffee/lck && ' + \
                                'mkdir -p ' + self.tmpPath + self.unique + '/tcoffee/cache'
                subprocess.check_call(mkdir_command, shell=True)

            # FoldX
            if not os.path.isdir(self.tmpPath + self.unique + '/FoldX'):
                # make the directories
                mkdir_command = 'mkdir -p ' + self.tmpPath + self.unique + '/FoldX'
                # copy the executables
                cp_command = 'cp ' + self.executables + 'FoldX.linux64 ' + self.tmpPath + self.unique + '/FoldX/ && ' + \
                           'cp ' + self.executables + 'rotabase.txt ' + self.tmpPath + self.unique + '/FoldX/'
                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
                # Copy dssp into the folder for modelling
                cp_command = 'cp ' + self.executables + 'dssp-2.0.4-linux-amd64 ' + self.tmpPath + self.unique + '/FoldX/'
                subprocess.check_call(cp_command, shell=True)

            # modeller
            if not os.path.isdir(self.tmpPath + self.unique + '/modeller'):
                # create workingfolder for modeller
                mkdir_command = 'mkdir -p ' + self.tmpPath + self.unique + '/modeller'
                subprocess.check_call(mkdir_command, shell=True)
                # Copy knot into the same folder as modeller
                cp_command = 'cp ' + self.executables + 'topol ' + self.tmpPath + self.unique + '/modeller'
                subprocess.check_call(cp_command, shell=True)

            # sequence conservation
            if not os.path.isdir(self.tmpPath + self.unique + '/sequence_conservation'):
                mkdir_command = 'mkdir -p ' + self.tmpPath + self.unique + '/sequence_conservation'
                subprocess.check_call(mkdir_command, shell=True)
                # provean
                cp_command = 'cp ' + self.executables + 'provean ' + self.tmpPath + self.unique + '/sequence_conservation/'
                subprocess.check_call(cp_command, shell=True)

            # analyze_structure
            if not os.path.isdir(self.tmpPath + self.unique + '/analyze_structure'):
                # create workingfolder for analyzing structure sasa and secondary structure
                mkdir_command = 'mkdir -p ' + self.tmpPath + self.unique + '/analyze_structure'
                subprocess.check_call(mkdir_command, shell=True)
                # Pops
                cp_command = 'cp ' + self.executables + 'pops ' + self.tmpPath + self.unique + '/analyze_structure'
                subprocess.check_call(cp_command, shell=True)
                # Dssp
                cp_command = 'cp ' + self.executables + 'dssp-2.0.4-linux-amd64 ' + self.tmpPath + self.unique + '/analyze_structure/dssp'
                subprocess.check_call(cp_command, shell=True)
#                # Naccess
#                cp_command = (
#                    'cp ' + self.executables + 'naccess ' + self.tmpPath + self.unique + '/analyze_structure/ && '
#                    'cp ' + self.executables + 'accall ' + self.tmpPath + self.unique + '/analyze_structure/ && '
#                    'cp ' + self.executables + 'standard.data ' + self.tmpPath + self.unique + '/analyze_structure/ && '
#                    'cp ' + self.executables + 'vdw.radii ' + self.tmpPath + self.unique + '/analyze_structure/')
#                subprocess.check_call(cp_command, shell=True)
                #MSMS
                cp_command = (
                    'cp ' + self.executables + 'pdb_to_xyzrn ' + self.tmpPath + self.unique + '/analyze_structure/ && '
                    'cp ' + self.executables + 'standard.data ' + self.tmpPath + self.unique + '/analyze_structure/ && '
                    'cp ' + self.executables + 'atmtypenumbers ' + self.tmpPath + self.unique + '/analyze_structure/ && '
                    'cp ' + self.executables + 'msms.x86_64Linux2.2.6.1 ' + self.tmpPath + self.unique + '/analyze_structure/')
                subprocess.check_call(cp_command, shell=True)

#            # create tmp for KNOT
#            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/KNOT'):
#                # make the directories
#                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/KNOT'
#                # copy the executables
#                cp_command = 'cp ' + self.executables + 'topol ' + self.tmpPath + 'Consumer-' + str(i) + '/KNOT'
#                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
#
#            # create tmp for pops
#            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/pops'):
#                # make the directories
#                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/pops'
#                # copy the executables
#                cp_command = 'cp ' + self.executables + 'pops ' + self.tmpPath + 'Consumer-' + str(i) + '/pops'
#                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
#
#            # create tmp folder for dssp
#                # make the directories
#                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/dssp'
#                # copy the executables
#                cp_command = 'cp ' + self.executables + 'dssp-2.0.4-linux-amd64 ' + self.tmpPath + 'Consumer-' + str(i) + '/dssp'
#                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)



if __name__ == '__main__':
    # read which configFile to use
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True)
    parser.add_argument('-i', '--input_file')
    parser.add_argument('-u', '--uniprot_id')
    parser.add_argument('-m', '--mutations', nargs='+')
    args = parser.parse_args()

    assert os.path.isfile(args.config_file)
    pipeline = Pipeline(args.config_file)

    uniprot_ids = []
    mutations = []
    if args.input_file \
    and os.path.isfile(args.input_file):
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
        mutations = args.mutations

    else:
        raise Exception('Need to supply either a list of uniprot_mutation combos '
            'or a flatfile with the same!')

    # Run jobs
    print args.uniprot_id
    print args.mutations
    for uniprot_id, mutation in zip(uniprot_ids, mutations):
        print uniprot_id
        print mutation
        pipeline(uniprot_id, mutation)



