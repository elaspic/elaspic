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

import class_error as error
import class_sql as sql

import class_domain_template
import class_domain_model
import class_domain_mutation

from Bio.SubsMat import MatrixInfo
from ConfigParser import SafeConfigParser
from scinetCleanup import scinetCleanup




class ProteinDomainData(object):
    
    def __init__(self, domain, template, model, mutations=[]):
        self.domain = domain
        self.template = template
        self.model = model
        self.mutations = mutations
        


class Pipeline(object):

    def __init__(self, configFile):
       
        #######################################################################
        # read the configuration file and set the variables
        configParser = SafeConfigParser(
            defaults={'tmpPath': '/tmp/pipeline/',
                      'HETATM': True,
                      'DEBUG': False,
                      'saveTo': '$SCRATCH/niklas-pipeline/',
                      'saveScinet': False,
                      'path_to_archive': '/home/kimlab1/database_data/elaspic/human/',
                      'db_type': sql.sql_flavor,
                      'db_path': error.path_to_pipeline_code() + '/../db/pipeline.db',
                      'look_for_interactions': 'True',
                      'webServer': False,
                      'numConsumers': '1',
                      'tcoffee_parallel_runs': '1',
                      'outputPath': '../results/',
                      'savePDB': 'results/pdb_files/',
                      'runTime': 'INFINITE',
                      'matrix': 'blosum80',
                      'gap_start': '-16',
                      'gap_extend': '-4'})

        configParser.read(configFile)

        # from [DEFAULT]
        tmpPath = configParser.get('DEFAULT', 'tmpPath')
        self.DEBUG = configParser.getboolean('DEFAULT', 'DEBUG')
        self.HETATM = configParser.getboolean('DEFAULT', 'HETATM')
        self.saveTo = configParser.get('DEFAULT', 'saveTo')
        self.saveScinet = configParser.getboolean('DEFAULT', 'saveScinet')
        self.path_to_archive = configParser.get('DEFAULT', 'path_to_archive')
        self.db_type = configParser.get('DEFAULT', 'db_type')
        self.db_path = configParser.get('DEFAULT', 'db_path')
        self.look_for_interactions = configParser.getboolean('DEFAULT', 'look_for_interactions')
        self.webServer = configParser.get('DEFAULT', 'webServer')
            
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
        
        if tmpPath[0] == '/':
            # i.e. tmpPath is given as an absolute Path
            self.tmpPath = self.TMPDIR_CLUSTER + tmpPath
        else:
            self.tmpPath = self.TMPDIR_CLUSTER + '/' + tmpPath
            
        
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
        path_to_local_pdbaa = '/home/kimlab1/strokach/ncbi-blast-2.2.28+/pdbaa_db'
        path_to_tmp_pdbaa = self.tmpPath + 'blast/pdbaa_db'
        if not os.path.isdir(path_to_tmp_pdbaa):
#            childProcess = subprocess.Popen('hostname | cut -d. -f1', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#            host_name, __ = childProcess.communicate()
            if whoami.strip() == 'strokach':
                # when running the pipeline on beagle or banting, copy the database
                # to a local folder
                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                    'cp -r ' + path_to_local_pdbaa + ' ' + path_to_tmp_pdbaa
            elif whoami.strip() == 'alexey':
                # when running the pipeline locally there is no need to copy the database
                # a symlink is enough
                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                    'cd ' + self.tmpPath + 'blast && ' + \
                                    'ln -sf /home/kimlab1/strokach/ncbi-blast-2.2.28+/pdbaa_db'
            elif whoami.strip() == 'witvliet':
                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                    'cd ' + self.tmpPath + 'blast && ' + \
                                    'ln -sf /home/witvliet/working/bin/ncbi-blast-2.2.28+/pdbaa_db'
            elif whoami.strip() == 'joan':
                # for scinet, blast is already installed, but the database needs to be copied
                system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                    'cp -ru $HOME/niklas-pipeline/blastdb/pdbaa_db ' + \
                                    self.tmpPath + 'blast/'
            
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, __ = childProcess.communicate()
            if childProcess.returncode != 0:
                raise Exception('Couldnt copy the blast database!')
                
        
        #######################################################################
        # Copy the database file
        if self.db_type == 'sqlite_file' and not os.path.isfile(self.tmpPath + '/' + self.db_path.strip().split('/')[-1]):
            print self.tmpPath
            print self.db_path
            system_command = 'cp -u ' + self.db_path + ' ' + self.tmpPath + '/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, __ = childProcess.communicate()
            if childProcess.returncode != 0:
                raise Exception('Could not copy the sqlite3 database!')
    
    
        #######################################################################
        # Create output folder if it doesn't exist
        if not os.path.isdir(self.outputPath):
            subprocess.check_call('mkdir -p ' + self.outputPath, shell=True)

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
        logger.addHandler(handler)
        self.log = logger

        # Initialise the sql database for accessing all information
        self.db = sql.MyDatabase(path_to_sqlite_db=(self.tmpPath + 'pipeline.db'),
                                 sql_flavor=self.db_type, 
                                 is_immutable=False,
                                 path_to_temp=self.tmpPath,
                                 path_to_archive=self.path_to_archive)
        
        # go to a unique directory
        self.PWD = self.tmpPath + self.unique + '/'
        os.chdir(self.PWD)
        
        # Initialize the classes used for calculating templates, models and mutations
        self.get_template = class_domain_template.GetTemplate(self.tmpPath, self.unique, self.pdbPath, self.db, self.log)
            
        self.get_model = class_domain_model.GetModel(self.tmpPath, self.unique, self.pdbPath, self.db, self.log, 
            self.modeller_runs, self.buildModel_runs, self.PWD)
            
        self.get_mutation = class_domain_mutation.GetMutation(self.tmpPath, self.unique, self.pdbPath, self.db, self.log, 
            self.foldX_WATER, self.buildModel_runs, self.PWD, self.matrix, self.gap_start, self.gap_extend)
        
        
        #######################################################################
        # Find all domains and domain pairs for a given uniprot
        protein_domains = self.db.get_uniprot_domain(self.uniprot_id)
        if self.look_for_interactions:
            protein_domain_pairs = self.db.get_uniprot_domain_pair(self.uniprot_id)
            protein_domain_template_model = protein_domains + protein_domain_pairs
        else:
            protein_domain_template_model = protein_domains
        if protein_domain_template_model == []:
            self.log.error('Uniprot %s has no pfam domains' % self.uniprot_id)
        # Convert to a single object which holds domain, template, model, and mutations
        # (so that I can modify each of those parameters "in place" in a for loop)
        self.protein_definitions = []
        for domain, template, model in protein_domain_template_model:
            self.protein_definitions.append(ProteinDomainData(domain, template, model))        
        
        
        #######################################################################
        # Do the calculations
        self.log.info("Finding templates...")
        self._compute_templates()
        
        self.log.info("Building models...")
        self._compute_models()

        self.log.info("Analyzing mutations...")        
        self._compute_mutations()
            
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
        for d in self.protein_definitions:
            self.log.info('%s\t%s\t%s' % (d.domain, d.template, d.model,) )
            
            if d.template:
                if False:
                    pass # So that I can comment out all the elifs
#                elif d.model and d.model.model_errors:
#                    # Recalculate templates if the model has errors
#                    self.log.debug('Recalculating templates because the model had errors')
#                    pass
                elif type(d.template) == sql.UniprotDomainPairTemplate and \
                d.template.domain_1 and d.template.domain_2 and \
                d.template.domain_1 == d.template.domain_2 and \
                d.template.domain_1.pdb_chain == d.template.domain_2.pdb_chain:
                    # Recalculate cases where the two chains are the same
                    self.log.debug('Recalculating templates because both interacting partners mapped to the same pdb chain')
                    pass
#                elif not d.template.template_errors == 'no templates found':
#                    # Recalculate cases where a template was found but errors occured
#                    self.log.debug('Recalculating templates because odd errors occured the last time we tried')
#                    pass
                else:
                    self.log.info('skipping template')
                    continue
            
            # Construct empty templates that will be used if we have errors
            if type(d.domain) == sql.UniprotDomain:
                empty_template = sql.UniprotDomainTemplate()
                empty_template.uniprot_domain_id = d.domain.uniprot_domain_id
            elif type(d.domain) == sql.UniprotDomainPair:
                empty_template = sql.UniprotDomainPairTemplate()
                empty_template.uniprot_domain_pair_id = d.domain.uniprot_domain_pair_id
            
            # Find templates and catch errors
            try:
                template = self.get_template(d.domain)
            except (error.NoStructuralTemplates,
                    error.NoSequenceFound,
                    error.NoTemplatesFound,
                    error.TcoffeeError) as e:
                self.log.error(e.error)
                empty_template.template_errors = e.error
                template = empty_template
            except:
                raise
            
            self.log.info('adding template')
            d.template = template
            self.db.add_uniprot_template(template, d.domain.path_to_data)
                
        self.log.info('Finished processing all templates for ' + self.uniprot_id + ' ' + self.mutations + '\n')



    def _compute_models(self):            
        """ 
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        
        # Go over all domains and domain pairs for a given protein
        for p in self.protein_definitions:
            self.log.info('%s\t%s\t%s' % (p.domain, p.template, p.model,) )
            
            if p.model and not p.model.model_errors:
                # Have models from a previous run
                self.log.info('skipping model')
                continue
            
            # Construct empty models that will be used if we have errors
            if type(p.domain) == sql.UniprotDomain:
                empty_model = sql.UniprotDomainModel()
                empty_model.uniprot_domain_id = p.template.uniprot_domain_id
            elif type(p.domain) == sql.UniprotDomainPair:
                empty_model = sql.UniprotDomainPairModel()
                empty_model.uniprot_domain_pair_id = p.template.uniprot_domain_pair_id
            
            
            if (type(p.template) == sql.UniprotDomainTemplate and not p.template.domain_def) \
            or (type(p.template) == sql.UniprotDomainPairTemplate and (not p.template.domain_def_1 or not p.template.domain_def_2)):
                # Template has errors and no domains, so give up on the modelling
                self.log.error('Error with the template')
                empty_model.error = 'error with the template'
                uniprot_model = empty_model
            else:
                # Make a homology model using the domain and template information
                try:
                    uniprot_model = self.get_model(p.domain, p.template)
                except error.PDBChainError as e:
                    self.log.error(e.error)
                    uniprot_model = empty_model
                    uniprot_model.model_errors = e.error
                except ModellerError as e:
                    self.log.error('ModellerError while trying to modellel in __getPDB for ' + self.uniprot_id + ':' + self.mutations)
                    self.log.error('ModellerError args:' + '\n'.join(e.args))
                    self.log.error('ModellerError message:' + e.message)
                    if 'Alignment sequence not found in PDB file' in e.message:
                        empty_model.model_errors = 'alignment sequence not found in PDB file'
                    else:
                        empty_model.model_errors = 'modelling error with message: %s' % e.message.replace('\n','; ')
                    uniprot_model = empty_model
            
            self.log.info('adding model')
            p.model = uniprot_model
            self.db.add_uniprot_model(uniprot_model, p.domain.path_to_data)
        
        self.log.info('Finished processing all models for ' + self.uniprot_id + ' ' + self.mutations + '\n')



    def _compute_mutations(self):
        """
        """
        
        for p in self.protein_definitions:
            self.log.info('%s\t%s\t%s' % (p.domain, p.template, p.model) )
            for mutation in self.mutations.split(','):
                
                precalculated_mutation = self.db.get_uniprot_mutation(p.model, mutation)
                self.log.info(precalculated_mutation)
                if precalculated_mutation:
                    p.mutations.append(precalculated_mutation)
                    self.log.info('skipping mutation')
                    continue
                
                # Construct empty models that will be used if we have errors
                if type(p.domain) == sql.UniprotDomain:
                    if not p.template.domain_def:
                        continue
                    empty_mutation = sql.UniprotDomainMutation()
                    empty_mutation.uniprot_domain_id = p.domain.uniprot_domain_id
                elif type(p.domain) == sql.UniprotDomainPair:
                    if not p.template.domain_def_1 or not p.template.domain_def_2:
                        continue
                    empty_mutation = sql.UniprotDomainPairMutation()
                    empty_mutation.uniprot_domain_pair_id = p.domain.uniprot_domain_pair_id
                
                try:
                    uniprot_mutation = self.get_mutation(p.domain, p.template, p.model, self.uniprot_id, mutation)
                except (error.MutationOutsideDomain,
                        error.MutationOutsideInterface) as e:
                    self.log.debug(e.error)
                    continue
                except error.FoldXError as e:
                    self.log.error('FoldXError while repairing the wildtype for ' + self.uniprot_id + ':' + self.mutations)
                    self.log.error('FoldXError error:' + e.error)
                    continue
                except error.pdbError as e:
                    self.log.error(e.error)
                    continue
                except error.NotInteracting as e:
                    self.log.error('Uniprot 1: %s, uniprot domain pair id: %s, Mutation %s not at interface!' 
                        % (self.uniprot_id, p.domain.uniprot_domain_pair_id, mutation))
                    continue
                
                self.log.info('adding mutation %s' % mutation)
                p.mutations.append(uniprot_mutation)
                self.db.add_uniprot_mutation(uniprot_mutation, p.domain.path_to_data)
      
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
                subprocess.check_call('mkdir ' + self.tmpPath + self.unique, shell=True)
                
            # tcoffee
            if not os.path.isdir(self.tmpPath + self.unique + '/tcoffee'):
                # create tmp for tcoffee
                mkdir_command = 'mkdir ' + self.tmpPath + self.unique + '/tcoffee && ' + \
                                'mkdir ' + self.tmpPath + self.unique + '/tcoffee/tmp && ' + \
                                'mkdir ' + self.tmpPath + self.unique + '/tcoffee/lck && ' + \
                                'mkdir ' + self.tmpPath + self.unique + '/tcoffee/cache'
                subprocess.check_call(mkdir_command, shell=True)
        
            # FoldX
            if not os.path.isdir(self.tmpPath + self.unique + '/FoldX'):
                # make the directories
                mkdir_command = 'mkdir ' + self.tmpPath + self.unique + '/FoldX'
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
                mkdir_command = 'mkdir ' + self.tmpPath + self.unique + '/modeller'
                subprocess.check_call(mkdir_command, shell=True)
                # Copy knot into the same folder as modeller
                cp_command = 'cp ' + self.executables + 'topol ' + self.tmpPath + self.unique + '/modeller'
                subprocess.check_call(cp_command, shell=True)
            
            # analyze_structure
            if not os.path.isdir(self.tmpPath + self.unique + '/analyze_structure'):
                # create workingfolder for analyzing structure sasa and secondary structure
                mkdir_command = 'mkdir ' + self.tmpPath + self.unique + '/analyze_structure'
                subprocess.check_call(mkdir_command, shell=True)                
                # Pops
                cp_command = 'cp ' + self.executables + 'pops ' + self.tmpPath + self.unique + '/analyze_structure'
                subprocess.check_call(cp_command, shell=True)
                # Dssp
                cp_command = 'cp ' + self.executables + 'dssp-2.0.4-linux-amd64 ' + self.tmpPath + self.unique + '/analyze_structure/dssp'
                subprocess.check_call(cp_command, shell=True)
                # Naccess
                cp_command = (
                    'cp ' + self.executables + 'naccess ' + self.tmpPath + self.unique + '/analyze_structure/ && '
                    'cp ' + self.executables + 'accall ' + self.tmpPath + self.unique + '/analyze_structure/ && '
                    'cp ' + self.executables + 'standard.data ' + self.tmpPath + self.unique + '/analyze_structure/ && '
                    'cp ' + self.executables + 'vdw.radii ' + self.tmpPath + self.unique + '/analyze_structure/')
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
    parser.add_argument('-u', '--uniprot_ids', nargs='+')
    arguments = parser.parse_args()
           
    assert os.path.isfile(arguments.config_file)
    pipeline = Pipeline(arguments.config_file)
    
    uniprot_ids = []
    mutations = []    
    if arguments.input_file \
    and os.path.isfile(arguments.input_file):
        with open(arguments.input_file, 'r') as fh:
            for line in fh:
                # Can skip lines by adding spaces or tabs before them
                if line[0][0] == ' ' or line[0][0] == '\t':
                    continue
                
                row = [ l.strip() for l in line.split('\t') ]
                
                # AS: Mutation does not necessarily have to be specified
                if len(line) > 1:
                    uniprot_id, mutation = line[0], line[1]
                elif len(line) == 1:
                    uniprot_id = line[0]
                    mutation = ''
                #
                uniprot_ids.append(uniprot_id)
                mutations.append(mutation)
    
    elif arguments.uniprot_ids:
        for uniprot_mut in arguments.uniprot_ids:
            uniprot_mut = uniprot_mut.split('_')
            if len(uniprot_mut) == 1:
                # No mutation was specified
                uniprot_mut.append('')
            uniprot_id, mutation = uniprot_mut
            #
            uniprot_ids.append(uniprot_id)
            mutations.append(mutation)
    
    else:
        raise Exception('Need to supply either a list of uniprot_mutation combos '
            'or a flatfile with the same!')
    
    # Run jobs
    for uniprot_id, mutation in zip(uniprot_ids, mutations):
        print uniprot_id
        print mutation
        pipeline(uniprot_id, mutation)
    
    
    