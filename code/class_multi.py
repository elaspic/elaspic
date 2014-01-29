# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 18:37:36 2012

@author: niklas
"""
import os
import logging
import multiprocessing

from class_time import runTime as rt
from modeller import ModellerError

import class_error as error
import class_sql as sql

import class_domain_template
import class_domain_model
import class_domain_mutation


class Consumer(multiprocessing.Process):
   
    def __init__(self, proc_name, task_queue, result_queue, runTime, pool, semaphore, DEBUG, outputPath, logger, webServer=False):
        multiprocessing.Process.__init__(self)
        self.proc_name = proc_name
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.DEBUG = DEBUG
        self.webServer = webServer
        
        self.pool = pool
        self.semaphore = semaphore
        self.outputPath = outputPath
        self.logger = logger
        
        # get the logger from the parent and add a handler
        self.log = logging.getLogger(proc_name)
        self.log.addHandler(self.logger)
        self.log.setLevel(logging.DEBUG)

        # to check how long the program has been running
        self.now = rt()
        self.runTime = self.setTime(runTime)
    
    
    def setTime(self, runTime):
        """
        input of the form h:m or 'INFINITE'
        """
        if runTime == 'INFINITE':
            d = 24*60*60
            return self.now() + 4*d
        else:
            h, m = runTime.split(':')
            return 60*60*int(h) + 60*int(m)
    
    # credit goes to here:
    # http://code.activestate.com/recipes/473878-timeout-function-using-threading/
    # not really clean... the child process remain running. For the purpose of
    # running them on the cluster this is not a problem since they will be
    # killed automatically
    def timeout(self, func, args=(), kwargs={}, timeout_duration=1, default='timeout'):
        """
        runs func for at most timeout_duration seconds
        returns default after timeout_duration seconds
        """
        class InterruptableThread(multiprocessing.Process):
            def __init__(self, sender):
                multiprocessing.Process.__init__(self)
                self.result = 'Null'
                self.sender = sender
    
            def run(self):
                self.result = func(*args, **kwargs)
                self.sender.send(self.result)
        
        sender, receiver = multiprocessing.Pipe()
        
        it = InterruptableThread(sender)
        it.start()
        it.join(timeout_duration)
        if it.is_alive():
            it.terminate()
            print 'returning default', default
            return default
        else:
            return receiver.recv()


    def run(self):
        i = 0
        while True:
            i += 1
            next_task = self.task_queue.get()
            # check if it was the last task
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
                        
            # check if the next calculation should be started
            # when run on a cluster the runtime is limited and to avoid that
            # the program is killed this check is made
            # the limit on scinet is 48 hours, let 8 hours for the last calculation
            # and the cleanup, i.e. 40 hours runtime, i.e. 144 000 seconds

            # do the calculations if time permits
            
            remainingTime = self.runTime - float(self.now())
            if remainingTime >= 60:
                answer = self.timeout(func=next_task, 
                                      kwargs={'unique':self.proc_name,
                                              'semaphore':self.semaphore,
                                              'pool':self.pool,
                                              'log':self.log,
                                              'DEBUG':self.DEBUG},
                                      timeout_duration=(remainingTime-60))
                self.result_queue.put(answer)
            else:
                self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'timeout'])
            self.task_queue.task_done()

        print 'exiting ', self.proc_name
        return



class Task(object):

    def __init__(self, uniprot_id, mutation, # UNIQUE FOR EACH JOB
            template_finding, savePDB, tmpPath, outputPath, pdbPath,
            matrix, gap_start, gap_extend, modeller_runs, buildModel_runs,
            foldX_WATER, path_to_archive, db_type):
        """ Initiate the variables that will be used throughout template finding
        and mutation modelling
        """
        self.uniprot_id = uniprot_id
        self.mutations = mutation
        
        self.template_finding = template_finding
        self.matrix = matrix # which matrix is used for the interface similarity
                             # calculation
        self.gap_s = gap_start  # gap penalty gap start
        self.gap_e = gap_extend # gap penalty gap extend
        
        self.HOME = os.getcwd() + '/'
        
        # savePDB is actually only used by __run as hidden input for prepareInput
        self.savePDB = self.HOME + savePDB
        self.tmpPath = tmpPath
        self.outputPath = outputPath
        self.pdbPath = pdbPath
        self.path_to_archive = path_to_archive    
        
        self.modeller_runs = modeller_runs # how many modells should be produced by modeller
        self.buildModel_runs = buildModel_runs
        
        self.foldX_WATER = foldX_WATER # water atoms in foldX
        self.db_type = db_type
        # stuff like modeller_path can't be specified here since self.unique
        # is not yet set. This happens after initialising the Task for each
        # Consumer.
        self.unique = 'proc_name'
#        self.semaphore = None
#        self.pool = None
#        self.log = None
#        self.DEBUG = None
        

    def __call__(self, unique, semaphore, pool, log, DEBUG):
        """ Run the main function of the program and parse errors
        """
        print self.uniprot_id
        
        #######################################################################
        # Usually all the self variables would be set in __init__, but we don't
        # know the unique process name and sephamor eat that point...
        self.unique = unique
        self.semaphore = semaphore
        self.pool = pool
        self.log = log
        self.DEBUG = DEBUG
        
        # Initialise the sql database for accessing all information
        self.db = sql.MyDatabase(path_to_sqlite_db = self.tmpPath + 'pipeline.db',
                                 sql_flavor = self.db_type, 
                                 is_immutable = False,
                                 path_to_temp = self.tmpPath + self.unique + '/',
                                 path_to_archive = self.path_to_archive)
        self.save_to_db = True
        
        # Set variables that require the correct "self.unique" id
        self.alignments_path = self.tmpPath + self.unique + '/tcoffee/'
        
        # go to a unique directory
        self.PWD = self.tmpPath + self.unique + '/'
        os.chdir(self.PWD)
        
        self.get_template = class_domain_template.GetTemplate(self.tmpPath, self.unique, self.pdbPath,
            self.savePDB, self.alignments_path, self.pool, self.semaphore, 
            self.db, self.log)
            
        self.get_model = class_domain_model.GetModel(self.tmpPath, self.unique, self.pdbPath,
            self.savePDB, self.alignments_path, self.pool, self.semaphore, 
            self.db, self.log,
            self.modeller_runs, self.buildModel_runs, self.PWD)
            
        self.get_mutation = class_domain_mutation.GetMutation(self.tmpPath, self.unique, self.pdbPath,
            self.savePDB, self.alignments_path, self.pool, self.semaphore, 
            self.db, self.log,
            self.foldX_WATER, self.buildModel_runs, self.PWD, self.matrix, self.gap_s, self.gap_e)
        
        # Core and interface domain definitions
        protein_domains = self.db.get_uniprot_domain(self.uniprot_id)
        
        protein_domain_pairs = self.db.get_uniprot_domain_pair(self.uniprot_id)
        # Duplicates will occur for homodimers (ie forward and reverse are the same)
        protein_domain_pairs = list(set(protein_domain_pairs[0] + protein_domain_pairs[1]))
        
        #======================================================================
        # Don't forget to uncomment
#        protein_definitions = protein_domain_pairs + protein_domains
        protein_definitions = protein_domains
        #======================================================================

        if protein_definitions == []:
            if self.DEBUG:
                raise error.ProteinDefinitionError('The protein has no pfam domains')
            else:
                return
        
        self.protein_definitions = protein_definitions
        self.protein_mutations = None
        
        #######################################################################
        # Do the calculations
        
        self.compute_templates()
        self.log.info("Finished getting core templates...\n")
        
        self.compute_models()
        self.log.info("Finished making core homology models...\n")
        
#        self.compute_mutations()
#        self.log.info("Finished evaluating mutations...\n")
            
        # Need to add some stuff for foldX: <analyse complex>
        
        #######################################################################
        # Leave everything the way you found it
        self.db.session.close()
        os.chdir(self.PWD)
        
        return [self.protein_definitions, self.protein_mutations]
        
            
        
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

    
    def compute_templates(self):
        """ Align uniprot domains to structural domains in pdbfam, find the best 
        template, and expand domain boundaries
        """
        protein_definitions_with_templates = []
        for uniprot_domain, uniprot_template, uniprot_model in self.protein_definitions:
            print uniprot_domain, uniprot_template, uniprot_model                
                
            # Check if tired expanding boundaries for this domain previously, and got erros
            if uniprot_template:
                protein_definitions_with_templates.append((uniprot_domain, uniprot_template, uniprot_model,))
                print 'skipping'
                continue
                
            # Construct empty templates that will be used if we have errors
            if type(uniprot_domain) == sql.UniprotDomain:
                empty_template = sql.UniprotDomainTemplate()
                empty_template.uniprot_domain_id = uniprot_domain.uniprot_domain_id
            elif type(uniprot_domain) == sql.UniprotDomainPair:
                empty_template = sql.UniprotDomainPairTemplate()
                empty_template.uniprot_domain_pair_id = uniprot_domain.uniprot_domain_pair_id
            
            
            # For domains that are not expanded and do not have the template cath_id:
            try:
                template = self.get_template(uniprot_domain)
                print "template:"
                
            except error.NoStructuralTemplates as e:
                if self.DEBUG:
                    raise e
                else:
                    empty_template.template_errors = e.error
                    template = empty_template
            except error.NoSequenceFound as e:
                if self.DEBUG:
                    raise e
                else:
                    empty_template.template_errors = e.error
                    template = empty_template
            except error.TemplateCoreError as e:
                self.log.error('Error while getting the core template for ' + self.uniprot_id + ':' + self.mutations)
                self.log.error('TemplateCoreError error:' + e.error)
                if self.DEBUG:
                    raise e
                else:
                    empty_template.template_errors = e.error
                    template = empty_template
            except error.TemplateInterfaceError as e:
                self.log.error('Error while getting the interface template for ' + self.uniprot_id + ':' + self.mutations)
                self.log.error('TemplateInterfaceError error:' + e.error)
                if self.DEBUG:
                    raise e
                else:
                    empty_template.template_errors = e.error
                    template = empty_template
            except error.TcoffeeError as e:
                self.log.error('t_coffee: problem occured while getting the wildtype PDB\n')
                self.log.error('t_coffee: alignment error in file: ' + e.alignInFile + '\n')
                self.log.error('t_coffee: message raised:\n')
                self.log.error('t_coffee:' + e.errors + '\n\n\n')
                if self.DEBUG:
                    raise e
                else:
                    empty_template.template_errors = e.error
                    template = empty_template
            except:
                if self.DEBUG:
                    raise
                else:
                    empty_template.template_errors = e.error
                    template = empty_template
            
            protein_definitions_with_templates.append((uniprot_domain, template, False,))
            self.db.add_uniprot_template(template, uniprot_domain.path_to_data)
                
        self.protein_definitions = protein_definitions_with_templates
        self.log.info('Finished processing all templates for ' + self.uniprot_id + ' ' + self.mutations + '\n\n\n')



    def compute_models(self):            
        """ Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        protein_definitions_with_models = []
        for uniprot_domain, uniprot_template, uniprot_model in self.protein_definitions:
            print uniprot_domain, uniprot_template, uniprot_model
            
            if uniprot_model:
                protein_definitions_with_models.append((uniprot_domain, uniprot_template, uniprot_model,))
                print 'skipping'
                continue
            
            # Construct empty models that will be used if we have errors
            if type(uniprot_domain) == sql.UniprotDomain:
                empty_model = sql.UniprotDomainModel()
                empty_model.uniprot_domain_id = uniprot_template.uniprot_domain_id
            elif type(uniprot_domain) == sql.UniprotDomainPair:
                empty_model = sql.UniprotDomainPairModel()
                empty_model.uniprot_domain_pair_id = uniprot_template.uniprot_domain_pair_id

            
#            # For domains that are not expanded and do not have the template cath_id:
#            try:
            uniprot_model = self.get_model(uniprot_domain, uniprot_template)
#
#            except error.KNOTerror as e:
#                if self.DEBUG:
#                    raise e
#                else:
#                    uniprot_model = empty_model
#            except ModellerError as e:
#                self.log.error('ModellerError while trying to modellel in __getPDB for ' + self.uniprot_id + ':' + self.mutations)
#                self.log.error('ModellerError args:' + '\n'.join(e.args))
#                self.log.error('ModellerError message:' + e.message)
#                if self.DEBUG:
#                    raise e
#                else:
#                    empty_model.modelling_errors = 'modelling error with message:' + e.message
#                    uniprot_model = empty_model                 
#            except Exception as e:
#                if self.DEBUG:
#                    raise e
#                else:
#                    # To do: have to make this part more explicit...
#                    empty_model.modelling_errors = 'unknown error'
##                    uniprot_model = empty_model
            
            protein_definitions_with_models.append((uniprot_domain, uniprot_template, uniprot_model,))
            self.db.add_uniprot_model(uniprot_model, uniprot_domain.path_to_data)
            
        self.protein_definitions = protein_definitions_with_models
        self.log.info('Finished processing all templates for ' + self.uniprot_id + ' ' + self.mutations + '\n\n\n')


    def compute_mutations(self):
                
        list_of_mutations_for_each_domain = []
        for uniprot_domain, uniprot_template, uniprot_model in self.protein_definitions:
            
            list_of_mutations = []
            for mutation in self.mutations.split(','):
                
                # Construct empty models that will be used if we have errors
                if type(uniprot_domain) == sql.UniprotDomain:
                    precalculated_mutation = self.db.get_uniprot_domain_mutation(uniprot_template.uniprot_domain_id, mutation)
                    empty_mutation = sql.UniprotDomainMutation()
                    empty_mutation.uniprot_domain_id = uniprot_domain.uniprot_domain_id
                    
                elif type(uniprot_domain) == sql.UniprotDomainPair:
                    precalculated_mutation = self.db.get_uniprot_domain_pair_mutation(uniprot_template.uniprot_domain_pair_id, mutation)
                    empty_mutation = sql.UniprotDomainPairMutation()
                    empty_mutation.uniprot_domain_pair_id = uniprot_domain.uniprot_domain_pair_id
                    
                if precalculated_mutation != []:
                    list_of_mutations.append(precalculated_mutation)
                    continue
                
                try:
                    uniprot_mutation = self.get_mutation(uniprot_domain, uniprot_template, uniprot_model, self.uniprot_id, mutation)
                    
                except error.FoldXError as e:
                    self.log.error('FoldXError while repairing the wildtype for ' + self.uniprot_id + ':' + self.mutations)
                    self.log.error('FoldXError error:' + e.error)
                    if self.DEBUG:
                        raise e
                    else:
                        return 'foldx_error'
                except error.pdbError as e:
                    self.log.error(e.error)
                    if self.DEBUG:
                        raise e
                    else:
                        return 'pdbError'
                except error.NotInteracting as e:
                    if self.DEBUG:
                        print ('Uniprot 1: %s, uniprot domain pair id: %s, Mutation %s not at interface!' 
                            % (self.uniprot_id, uniprot_domain.uniprot_domain_pair_id, mutation,))
                        continue
                    else:
                        print ('Uniprot 1: %s, uniprot domain pair id: %s, Mutation %s not at interface!' 
                            % (self.uniprot_id, uniprot_domain.uniprot_domain_pair_id, mutation,))
                        continue
                        
                list_of_mutations.append(uniprot_mutation)
                self.db.add_uniprot_mutation(uniprot_mutation, uniprot_domain.path_to_data)          
                    
            list_of_mutations_for_each_domain.append(list_of_mutations)
        
        self.protein_mutations = list_of_mutations_for_each_domain        
        self.log.info('Finished processing all templates for ' + self.uniprot_id + ' ' + self.mutations + '\n\n\n')
