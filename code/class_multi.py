# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 18:37:36 2012

@author: niklas
"""
from class_pdbTemplate import pdbTemplate
from class_time import runTime as rt
from class_interfaceSize import interfaceSize
from class_interfaceSize import getDSSP
from class_pysicoChemicalProperties import pysiChem
import class_modeller as mod
from class_foldX import foldX
import class_error as errors

import class_get_uniprot_template_core_and_interface as UniprotTemplates
from class_get_uniprot_template_core_and_interface import get_template_interface
from class_get_uniprot_template_core_and_interface import get_template_core

from HelperFunction_prepareModeller import prepareModeller

from modeller import ModellerError
import shutil
from os import getcwd, chdir
import logging
import multiprocessing
import subprocess

from Bio import SeqIO

import tarfile, os

import json

#from class_pdbTemplate import pdbClean
#import class_prepareModeller as prepMod
#from Bio.SubsMat import MatrixInfo
#import shlex
#from copy import deepcopy
#import sys
#import threading
#from Bio.Alphabet import IUPAC
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Seq import MutableSeq
#from Bio import AlignIO
#import class_callTcoffee as tc



class TimeoutException(Exception): 
    pass 

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
            # set the unique working folder
            next_task.unique = self.proc_name
            next_task.semaphore = self.semaphore
            next_task.pool = self.pool
            next_task.log = self.log
            
            # set debugging
            if self.DEBUG:
                next_task.DEBUG = True
            
            # check if the next calculation should be started
            # when run on a cluster the runtime is limited and to avoid that
            # the program is killed this check is made
            # the limit on scinet is 48 hours, let 8 hours for the last calculation
            # and the cleanup, i.e. 40 hours runtime, i.e. 144 000 seconds

            # do the calculations if time permits
            
            remainingTime = self.runTime - float(self.now())
            if remainingTime >= 60:
                answer = self.timeout(next_task, timeout_duration=(remainingTime-60))
                if answer == 'timeout':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'timeout'])
                elif answer == 'tcoffee_error':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'tcoffee_error'])
                elif answer == 'knotted':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'knotted'])
                elif answer == 'modeller_error':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'modeller_error'])
                elif answer == 'templateCoreError':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'templateCoreError'])
                elif answer == 'templateInterfaceError':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'templateInterfaceError'])
                elif answer == 'pdbError':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'pdbError'])
                elif answer == 'foldx_error':
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'foldx_error'])
                elif answer == 'unknown' or answer == None:
                    self.result_queue.put([None, next_task.uniprot_id + '_' + next_task.mutation, 'unknown'])
                else:
                    self.result_queue.put(answer)
            else:
                self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'timeout'])
            self.task_queue.task_done()

        print 'exiting ', self.proc_name
        return



class Task(object):

    def __init__(self,
                 templateFinding, # supplied as a pdb or a structure
                 uniprot_id, # UNIQUE FOR EACH JOB
                 mutation, # UNIQUE FOR EACH JOB
                 savePDB, 
                 tmpPath, 
                 outputPath,
                 pdbPath,
                 matrix, 
                 gap_start, 
                 gap_extend, 
                 modeller_runs,
                 buildModel_runs,
                 foldX_WATER,
                 DomainDefinition,
                 DomainInteraction,
                 UniprotSequence,
                 PDBResolution,
                 ProteinDefinition,
                 ProteinInteraction
                 ):
        """ Initiate the variables that will be used throughout template finding
        and mutation modelling
        """
        
        self.uniprot_id = uniprot_id
        
        self.templateFinding = templateFinding
        self.mutation = mutation
                
        self.matrix = matrix # which matrix is used for the interface similarity
                             # calculation
        self.gap_s = gap_start  # gap penalty gap start
        self.gap_e = gap_extend # gap penalty gap extend
        
        self.HOME = getcwd() + '/'
        
        # savePDB is actually only used by __run as hidden input for prepareInput
        self.savePDB = self.HOME + savePDB 
        self.tmpPath = tmpPath
        self.outputPath = outputPath
        self.pdbPath = pdbPath
        
        # stuff like modeller_path can't be specified here since self.unique
        # is not yet set. This happens after initialising the Task for each
        # Consumer.
        self.unique = 'proc_name'
        
        self.modeller_runs = modeller_runs # how many modells should be produced by modeller

        self.buildModel_runs = buildModel_runs
        self.foldX_WATER = foldX_WATER # water atoms in foldX
        
        # the database lookup functions
        self.UniprotSequence = UniprotSequence
        self.DomainDefinition = DomainDefinition
        self.DomainInteraction = DomainInteraction
        self.PDBResolution = PDBResolution
        self.ProteinDefinition = ProteinDefinition
        self.ProteinInteraction = ProteinInteraction
        
        self.pool = None
        self.semaphore = None
        
        self.log = None
        
        self.PWD = None
        
        self.DEBUG = False

        

    def __call__(self):
        """ Run the main function of the program and parse errors
        """
        
        # go to a unique directory
        chdir(self.tmpPath + self.unique)
        self.PWD = getcwd() + '/'
        
        try:
            return self.__run()
        
        # Take care of all possible errors
        except errors.KNOTerror as e:
            if self.DEBUG:
                raise e
            else:
                return 'knotted'
        except errors.TcoffeeError as e:
            self.log.error('t_coffee: problem occured while getting the wildtype PDB\n')
            self.log.error('t_coffee: alignment error in file: ' + e.alignInFile + '\n')
            self.log.error('t_coffee: message raised:\n')
            self.log.error('t_coffee:' + e.errors + '\n\n\n')
            if self.DEBUG:
                raise e
            else:
                return 'tcoffee_error'
        except ModellerError as e:
            self.log.error('ModellerError while trying to modellel in __getPDB for ' + self.uniprot_id + ':' + self.mutation)
            self.log.error('ModellerError args:' + '\n'.join(e.args))
            self.log.error('ModellerError message:' + e.message)
            if self.DEBUG:
                raise e
            else:
                return 'modeller_error'
        except errors.FoldXError as e:
            self.log.error('FoldXError while repairing the wildtype for ' + self.uniprot_id + ':' + self.mutation)
            self.log.error('FoldXError error:' + e.error)
            if self.DEBUG:
                raise e
            else:
                return 'foldx_error'
        except errors.TemplateCoreError as e:
            self.log.error('Error while getting the core template for ' + self.uniprot_id + ':' + self.mutation)
            self.log.error('TemplateCoreError error:' + e.error)
            if self.DEBUG:
                raise e
            else:
                return 'templateCoreError'
        except errors.TemplateInterfaceError as e:
            self.log.error('Error while getting the interface template for ' + self.uniprot_id + ':' + self.mutation)
            self.log.error('TemplateInterfaceError error:' + e.error)
            if self.DEBUG:
                raise e
            else:
                return 'templateInterfaceError'
        except errors.pdbError as e:
            self.log.error(e.error)
            if self.DEBUG:
                raise e
            else:
                return 'pdbError'
        except Exception as e:
            if self.DEBUG:
                raise e
            else:
                return 'unknown'
                
        finally:
            # Switch back to the original directory
            chdir(self.PWD)
    


    def __run(self):
        """ The workhorse of the class. Sequentially calls all functions required
        for precalculating structures, making mutations and collecting results.
        """
        
        if self.template_finding:
            # Find all templates for the given uniprot (and mutation)
            
            # Make a model for each domain in the query uniprot
            # (or use precalculated data if available)
                        
            # Save calculated templates to the core database folder
            
            
            
            # Make models of all possible interactions based on the uniprot,
            # uniprot domains, and pdbfam
            # (or use precalculated data if available)
                        
            # Save calculated templates to the interface database folder
                        
            templates = self.findTemplate()
    
            # if errors occured, return the errors and don't do further analysis
            if templates == [[]]:
                print 'no templates found'
                raise Exception
            elif templates[0].has_key('errors'):
                print 'no templates found'
                raise Exception
            
            if self.mutation != '':
                
                for template in templates:
                    
                    template_mutations = self.get_template_mutations(self.mutation)
                    if len(template_mutations['mutations']) == 0:
                        continue
                    
                    # Get the data for a particular mutation (or a set of mutations)
                    template_mutations = self.analyse_mutations(template, template_mutations)
    
                    template.update(template_mutations)
                
                self.log.info('Finished processing all templates for ' + self.uniprot_id + ' ' + self.mutation + '\n\n\n')
        
        
        if not self.template_finding and self.mutation != '':
            # Modell mutations using the structure from the pdb
        
            template, template_mutations = self.get_pdb_and_mutations(self.mutation)                                  
                                
            # Get the data for a particular mutation (or a set of mutations)
            template_mutations = self.analyse_mutations(template, template_mutations)
            
            templates = [template.update(template_mutations),]
        
        
        return templates      
                


    def score_pairwise(self, seq1, seq2, matrix, gap_s, gap_e):
        """ Required to get the BLOSUM score
        """
        
        def score_match(pair_match, matrix_match):
            if pair_match not in matrix_match:
                return matrix_match[(tuple(reversed(pair_match)))]
            else:
                return matrix_match[pair_match]
        
        score = 0
        gap = False
        for i in range(len(seq1)):
            pair = (seq1[i], seq2[i])
            if not gap:
                if '-' in pair:
                    gap = True
                    score += gap_s
                else:
                    score += score_match(pair, matrix)
            else:
                if '-' not in pair:
                    gap = False
                    score += score_match(pair, matrix)
                else:
                    score += gap_e
        return score
        
    
    def __getModel(self, 
                       alignments, 
                       target_ids, 
                       template_ids,
                       HETATMsInChain_SEQnumbering,
                       chains,
                       savePDB
                       ):
        outFile = self.tmpPath + self.unique + '/outFile_wildtype'
        
        # generate the input for modeller from the above generated alignment
        prepareModeller(outFile, 
                        alignments, 
                        target_ids, 
                        template_ids, 
                        HETATMsInChain_SEQnumbering,
                        chains
                        )
        
        inFile = outFile

        modeller_target_id = '_'.join(target_ids)
        
        if len(template_ids) == 1:
            modeller_template_id = template_ids[0]
        else:
            modeller_template_id = template_ids[0] + template_ids[1][-1] # if more than two chains are used this has to be adjusted

        modeller_path = self.tmpPath + self.unique + '/modeller/'
        chdir(modeller_path) # from os

        modeller = mod.modeller([inFile], 
                                modeller_target_id, 
                                modeller_template_id, 
                                savePDB, 
                                self.tmpPath + self.unique + '/', 
                                self.modeller_runs,
                                loopRefinement=True)
        normDOPE, pdbFile = modeller.run()

        chdir(self.PWD) # go back to the working directory
        
        return normDOPE, pdbFile

        
        
    def __getCrystalStructure(self, pdbCode, chains, savePDB, FULL_PATH=False):
        domains = [['Null', 'Null'] for i in range(len(chains)) ]
        pdb = pdbTemplate(self.pdbPath, pdbCode, chains, domains, savePDB)
        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()
        
        SWITCH_CHAIN = False
        if chains != chains_pdb_order:
            SWITCH_CHAIN = True
            chains = chains_pdb_order
        
        return '-', pdbCode.split('.')[0] + ''.join(chains) + '.pdb', SWITCH_CHAIN, chains

       
    def prepareMutationFoldX(self, sequence, mutations):
        """
        create the sequence snippets for the mutation with FoldX
        also, record the position of the mutation in the snippet
        """
        mutations_foldX = list()
        for mutation in mutations:
            if int(mutation[3:-1]) < 7:
                left = 0
                m_pos = int(mutation[3:-1]) # the position of the muation in the snippet
            else:
                left = int(mutation[3:-1])-7
                m_pos = 7
            if int(mutation[3:-1]) >= len(sequence) - 8:
                right = len(sequence)
            else:
                right = int(mutation[3:-1])+8

            wt_seq = sequence.seq[left:right]
            mut_seq = sequence.seq[left:int(mutation[3:-1])-1] + mutation[-1] + sequence.seq[int(mutation[3:-1]):right]

            assert( sequence.seq[int(mutation[3:-1])-1] == mutation[2] )
            
            mutations_foldX.append((m_pos, str(wt_seq) + '\n' + str(mut_seq)))

        return mutations_foldX


    def prepareInput(self, pdbCode, chains, domains, sequences, alignments, savePDB):

        pdb = pdbTemplate(self.pdbPath, pdbCode, chains, domains, savePDB)

        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract() # dict with chain ID as key and Nr. of HETATM as value

        SWITCH_CHAIN = False
        if chains != chains_pdb_order:
            SWITCH_CHAIN = True
            sequence_new_order = list()
            alignments_new_order = list()
            for item in chains_pdb_order:
                i = chains.index(item)
                sequence_new_order.append(sequences[i])
                alignments_new_order.append(alignments[i])
            sequences = sequence_new_order
            alignments = alignments_new_order
            
            chains = chains_pdb_order
        
        # to get the range of the domain make sure to call extract before getChainNumberingNOHETATMS!!
        chainNumberingDomain = dict()

        for i in range(0,len(chains)):
            chainNumberingDomain[chains[i]] = pdb.getChainNumbering(chains[i])

        HETATMsInChain_SEQnumbering = [ dict() for x in range(0,len(chains)) ]
        
        for i in range(0,len(chains)):
            HETATMsInChain_SEQnumbering[i][chains[i]] = list()
            for item in HETATMsInChain_PDBnumbering[chains[i]]:
                    try:
                        HETATMsInChain_SEQnumbering[i][chains[i]].append(chainNumberingDomain[chains[i]].index(item))
                    except ValueError:
                        HETATMsInChain_SEQnumbering[i][chains[i]].append(10001) # add one add the end that will be outside of the sequence

        return sequences, alignments, chains, SWITCH_CHAIN, HETflag, HETATMsInChain_SEQnumbering
    
    
    def setChainID(self, chains, mutations, HETflag, SWITCH_CHAIN, do_modelling):
        """
        modeller renames the chains, in order to keep track of the chain IDs
        they are renamed here to match the modeller output. Modeller labels the
        chains subsequently, i.e. A, B, C etc., thus, depending on wether HETATMs
        are present as seperate chain, the chain labels are A, B or A, C
        """
        mut = list()
        for mutation in mutations:
            mutChain = str(chains.index(mutation[0]))

            # get the chain names correct
            if do_modelling:
                if len(chains) == 1:
                    # if no HETATMs are present, the chain ID is empty
                    # otherwise the chain IDs are A and B for the atoms 
                    # and HETATMs respectively
                    if HETflag[chains[0]] == False:
                        chains_new = [' ']
                        mutChain = ' '
                    else:
                        chains_new = ['A']
                        mutChain = 'A'
                else:
                    a = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
                    index = 0
                    chains_new = list()
                    for chainID in chains:
                        if HETflag[chainID] == True:
                            # if True, then there are HETATMs associated to this chain
                            # in modeller they will be listed as two different chains.
                            # thus the chain names have to jump by two
                            chains_new.append(a[index])
                            index += 2
                        else:
                            chains_new.append(a[index])
                            index += 1
                    
                    if mutChain == '0':
                        mutChain = 'A'
                    elif mutChain == '1':
                        if chains_new[1] == 'B':
                            mutChain = 'B'
                        elif chains_new[1] == 'C':
                            mutChain = 'C'
                        else:
                            print 'could not assign the chain for mutation correctly!'
                            mutChain = None
                    
            mut.append( mutChain + '_' + mutation[2:] )
        return mut



    def findTemplate(self):
        # initialise the template selection classes
        # the unique variable is set after the __init__ is called, therefore
        # they need to be 'initialised' here.
    
        # Path where to put all the temporary alignments
        self.alignments_path = self.tmpPath + self.unique + '/tcoffee/'
        
        self.GetTemplateInterface = get_template_interface(self.tmpPath, 
                                                           self.unique, 
                                                           self.pdbPath, 
                                                           self.savePDB, 
                                                           self.alignments_path,
                                                           self.pool, 
                                                           self.semaphore, 
                                                           self.get_uniprot_sequence, 
                                                           self.interaction_database, 
                                                           self.threeDID_database, 
                                                           self.pdb_resolution_database,
                                                           self.log)
                                                           
        self.GetTemplateCore = get_template_core(self.tmpPath, 
                                                 self.unique, 
                                                 self.pdbPath, 
                                                 self.savePDB,
                                                 self.alignments_path,
                                                 self.pool, 
                                                 self.semaphore, 
                                                 self.get_uniprot_sequence, 
                                                 self.core_template_database,
                                                 self.pdb_resolution_database,
                                                 self.log)
        output_data = list()
        
        # Get templates for interface mutations
        templates, new_sequences = self.GetTemplateInterface(self.uniprot_id, self.mutation)
        self.log.info("Finished getting interface templates...\n\n")
#        self.log.debug("Interface templates:")
#        self.log.debug(templates)
        
        if not ((templates == ()) or (templates == [()]) or 
                (templates == []) or (templates == [[]])):
            # Interface mutations were found
            is_in_core = False
            
            # Collect structure information for best uniprot pairs modelling all unique pfam interactions
            for template in templates:
                output_dict = self.findTemplateHelper(template, is_in_core)
                output_data.append(output_dict)
            
        else:
            # No interface mutations were found
            is_in_core = True
            
            # Check if the mutations falls into the core
            templates, new_sequences = self.GetTemplateCore(self.uniprot_id, self.mutation)
            self.log.info("Finished getting core templates...\n\n")
#            self.log.debug("Core templates:")
#            self.log.debug(templates)
            
            if templates == 'not in core' or templates == 'no template':
                output_dict = {'errors': 'no template found: ' + str(self.uniprot_id) + '_' + str(self.mutation)}
                output_data.append(output_dict)
            elif templates == 'in gap':
                # add the logger here to report that the mutation did fall into a gap
                # in the alignment!
                output_dict = {'errors': 'no template found: ' + str(self.uniprot_id) + '_' + str(self.mutation)}
                output_data.append(output_dict)
            else:
                for template in templates:
                    output_dict = self.findTemplateHelper(template, is_in_core)
                    output_data.append(output_dict)

        output_data[0]['new_sequences'] = new_sequences
        
        self.log.info("Finished preparing template structures...\n\n")
#        self.log.debug("Output data:")
#        self.log.debug(output_data)
        
        return output_data


    def findTemplateHelper(self, template, is_in_core):
        # Unique template identifier for the template
        saveFolder = ('-'.join(template['uniprot_ids']) + '_' +
             '-'.join(template['pfam_names']) + '_' + 
             '_'.join(['-'.join([str(i) for i in x]) for x in template['domain_defs']]))
        
        # Folder for storing files for export to output
        savePath =  self.tmpPath + self.unique + '/' + saveFolder + '/'
        subprocess.check_call('mkdir -p ' + savePath, shell=True)
                    
        # Template data 
        pdbCode = template['pdb_id']
        chains_pdb = template['chains_pdb']
        mutations_pdb = [template['chains_pdb'][0] + '_' + self.mutation[0] + \
                            str(template['mutation_position_domain']) + self.mutation[-1], ]
        sequences   = template['uniprot_domain_sequences'] # is a list containing the uniprot sequences cut to the domain
        alignments  = template['alignments']
        
        # PDB domain definitons of query and partner
        domains_pdb = template['pdb_domain_defs']
        
        self.log.debug(savePath)
        self.log.debug("pdb_code:")
        self.log.debug(pdbCode)
        self.log.debug('chains_pdb:')
        self.log.debug(chains_pdb)
        self.log.debug('mutations_pdb:')
        self.log.debug(mutations_pdb)
        self.log.debug('sequences[0]:')
        self.log.debug(sequences[0])
        self.log.debug('sequences[1]:')
        self.log.debug(sequences[1])
        self.log.debug('alignments:')
        self.log.debug(alignments)
        self.log.debug('domains_pdb:')
        self.log.debug(domains_pdb)
        
        # check if all the templates have 100% sequence identity. If not,
        # set a flag so that the structure is created with modeller
        # this bit has to be implemented properly!!
        do_modelling = True   
#        for score in scores:
#        if float(scores[0]) >= 99.0:
#            do_modelling = False
        
        # savePath is where the pdb sequences and the pdb with the required chains are saved
        sequences, alignments, chains_modeller, switch_chain, HETflag, HETATMsInChain_SEQnumbering = \
            self.prepareInput(pdbCode, chains_pdb, domains_pdb, sequences, alignments, savePath)
        
        self.log.debug("sequences:")
        self.log.debug(sequences)
        self.log.debug("alignments:")
        self.log.debug(alignments)
        self.log.debug("chains_modeller:")
        self.log.debug(chains_modeller)
        self.log.debug("switch_chain:")
        self.log.debug(switch_chain)
        self.log.debug("HETflag:")
        self.log.debug(HETflag)
        self.log.debug("HETATMsInChain_SEQnumbering:")
        self.log.debug(HETATMsInChain_SEQnumbering)
        
        
        # copy the chain IDs
        # they are used afterwards for renaming and copy is needed so that
        # they are not overwritten
        target_ids   = list()
        template_ids = list()
        for i in range(0,len(alignments)):
            target_ids.append(alignments[i][0].id)
            template_ids.append(alignments[i][1].id)
        
        self.log.debug("target_ids:")
        self.log.debug(target_ids)
        self.log.debug("template_ids:")
        self.log.debug(template_ids)
        
        
        if do_modelling:
            normDOPE_wt, pdbFile_wt = self.__getModel(alignments, target_ids, template_ids, HETATMsInChain_SEQnumbering, chains_modeller, savePath)
            modeller_path = self.tmpPath + self.unique + '/modeller/'
        else:
            # add a function to take the crytal structure. Compare the __run()
            # method for a possibility. 
            pass

        # Copy Modeller pdb file
        shutil.copyfile(modeller_path + pdbFile_wt, savePath + pdbFile_wt)

        self.log.debug("normDOPE_wt:")
        self.log.debug(normDOPE_wt)
        self.log.debug("pdbFile_wt:")
        self.log.debug(pdbFile_wt + '\n\n')
        
        # Rename the wildtype pdb file
#        pdbFile_wt_renamed = self.uniprot_id + '_' + self.mutation + '.pdb'
#        system_command = 'mv ' + modeller_path + pdbFile_wt + ' ' + modeller_path + pdbFile_wt_renamed
#        subprocess.check_call(system_command, shell=True)
        
        # Copy alignments
        for alignment in template['alignments']:
            shutil.copyfile(self.alignments_path + alignment[0].id + '_' + alignment[1].id + '.aln',
                savePath + '/' + saveFolder + '_' + alignment[0].id + '_' + alignment[1].id + '.aln')
        
        # Add all calculated values to the template dictionary
        template['is_in_core'] = is_in_core
        template['saveFolder'] = saveFolder
        template['savePath'] = savePath
        template['pdbFile_wt'] = pdbFile_wt
        template['chains_modeller'] = chains_modeller
        template['mutations_pdb'] = mutations_pdb
        template['HETflag'] = HETflag
        template['normDOPE_wt'] = normDOPE_wt
        template['switch_chain'] = switch_chain

        return template



    def get_pdb_and_mutations(self, mutations):
        # mutations is of the form I_V70A (2VIR_AB_C_TC131I, 1PGA_A__VA29P)
        pdbCode, chains_pdb1, chains_pdb2, mutations = mutations.split('_')
        pdb_type, pdb_resolution = self.pdb_resolution_database(pdbCode)
        
        # Equivalent to class_get_uniprot_template_core_and_interface.py output        
        template = {'uniprot_ids': ('', '',),
                    'pfam_names': ('', '',),
                    'domain_defs': ('', '',),
                    'pdb_id': pdbCode,
                    'chains_pdb': [chains_pdb1, chains_pdb2],
                    'pdb_domain_defs': ((), (),),
                    'pdb_type': pdb_type,
                    'pdb_resolution': pdb_resolution,
                    'uniprot_domain_sequences': ('', '',),
                    'alignments': ('', '',),
                    'alignment_scores': (100, 100, 100,)}
        
        # Unique template identifier for the template
        saveFolder = (template['pdb_id'] + '_' + '-'.join(template['chains_pdb']))
        savePath =  self.tmpPath + self.unique + '/' + saveFolder + '/'
        subprocess.check_call('mkdir -p ' + savePath, shell=True)
        
        chains_get_pdb = [ item for item in chains_pdb1 ]
        chains_get_pdb.extend( [ item for item in chains_pdb2 ] )
        
        normDOPE_wt, pdbFile_wt, SWITCH_CHAIN, chains_get_pdb = self.__getCrystalStructure(pdbCode, chains_get_pdb, savePath)
        
        template_mutations = {'mutation_positions': [], 'mutations': [], 'mutations_pdb': [], 'mutations_modeller': []}
        for mutation in mutations.split(','):
    #        mutations_pdb = [mutation[1] + '_' + mutation[0] + mutation[2:], ] # So its [Chain_AAwt10AAmut, ]
            # The section below is required for making foldx_mutations
            sequence, chainNumbering = self._get_pdb_sequence(pdbCode, mutation[1])
            mutation_pos = chainNumbering.index(int(mutation[2:-1]) + 1)  # +1 to have it start with one
            mutation_pdb = mutation[1] + '_' + mutation[0] + str(mutation_pos) + mutation[-1]
            
            template_mutations['mutation_positions'].append(int(mutation[2:-1]))
            template_mutations['mutations'].append(mutation)
            template_mutations['mutations_pdb'].append(mutation_pdb)
            template_mutations['mutations_modeller'].append(mutation_pdb)
        
        # Equivalent to class_multi.py additions
        template['is_in_core'] = True
        template['normDOPE_wt'] = normDOPE_wt
        template['saveFolder'] = saveFolder
        template['savePath'] = savePath
        template['pdbFile_wt'] = pdbFile_wt
        template['sequences'] = (sequence,)
        template['chains_modeller'] = chains_get_pdb
        template['switch_chain'] = False
                
        return template, template_mutations


    
    def _get_pdb_sequence(self, pdbCode, chain):
        """
        Return the pdb file sequence (not SEQRES)
        """
        pdbPath = self.pdbPath
        domains = [['Null', 'Null'], ] # extract the full sequence

        pdb = pdbTemplate(pdbPath, pdbCode, chain, domains, self.tmpPath + self.unique + '/')
        
        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()
        chainNumberingDomain = pdb.getChainNumberingNOHETATMS(chain)
        if chainNumberingDomain == []:
            raise errors.pdbError('Could not get the pdb numbering for ' + pdbCode + '_' + chain)
        
        return next(SeqIO.parse(self.tmpPath + self.unique + '/' + pdbCode + chain + '.seq.txt', 'fasta')), chainNumberingDomain



    def get_template_mutations(self, template, mutations):
        
        template_mutations = {'mutation_positions': [], 'mutations': [], 'mutations_pdb': [], 'mutations_modeller': []}
        
        for mutation in mutations.split(','):
            # mutation numbered from the beginning of uniprot domain, rather than the beginning of uniprot
            mutation_position_domain = int(mutation[1:-1]) - template['domain_defs'][0][0] + 1 # +1 to have the numbering staring with 1
            
            mutation_pdb_no_chain = UniprotTemplates.map_to_pdb_sequence(template['alignments'][0],
                                                                template['alignment_ids'][0],
                                                                mutation_position_domain)
            if mutation_pdb_no_chain == 'in gap':
                continue
            
            mutation_pdb = template['chains_pdb'][0] + '_' + mutation[0] + str(mutation_pdb_no_chain) + mutation[-1]
            
            contacts_chain1 = UniprotTemplates.check_structure(template['pdb_id'],
                                                               template['chains_pdb'][0],
                                                               mutation_pdb)
            if not contacts_chain1[template['chains_pdb'][1]]:
                continue
            
            # Rename chains in mutation by modeller output
            mutation_modeller = self.setChainID(template['chains_modeller'], template['mutations_pdb'], template['HETflag'], template['switch_chain'], True)
            
            template_mutations['mutation_positions'].append(int(mutation[1:-1]))
            template_mutations['mutations'].append(mutation)
            template_mutations['mutations_pdb'].append(mutation_pdb)
            template_mutations['mutations_modeller'].append(mutation_modeller)
            
            self.log.debug("mutation_pdb:")
            self.log.debug(mutation_pdb)             
            self.log.debug("mutation_modeller:")
            self.log.debug(mutation_modeller)     
            
        return template_mutations



    def analyze_mutations(self, template, template_mutations):
        """
        AS: renamed the __run() function with some restructuring.
        
        input:  savePath        path to PDB files for modeller
                self.tmpPath    path to store some tmp data (not tcoffee)
                self.unique     set environment for tcoffee (for parallelization)
                pdbFile_wt      name of the wildtype structure file (either raw or from modeller)
                mutations_pdb   list of mutations in `Chain`_`AAwt``AAnumber``AAmut` format, where AAnumber is 
                                the index of the AA in the pdb
                chains_modeller     list of chains in pdbFile_wt, with the first chain in the list being mutated
                mutations_modeller  same as above, only compensated for chain switching and other changes made by modeller
                sequences       list of sequences of each chain in the pdb
                switch_chain    whether or not the order of the sequences was reversed during modelling
        """
        
        savePath, pdbFile_wt, chains, sequences, switch_chain = (
            template['savePath'], template['pdbFile_wt'], template['chains'],
            template['sequences'], template['switch_chain'], )
            
        mutations_pdb, mutations_modeller = (
            template_mutations['mutations_pdb'],
            template_mutations['mutations_modeller'], )
            
            
        #######################################################################
        # Logging
        self.log.debug("savePath:")
        self.log.debug(template['savePath'])          
        self.log.debug("pdbFile_wt:")
        self.log.debug(template['pdbFile_wt'])
        self.log.debug("chains:")
        self.log.debug(template['chains_modeller'])
        
        
        #######################################################################
        ## 1st: get a snippet of AAs around the mutation, as well as the mutation residue (no chain)
        # This is required input for FoldX
        if not switch_chain:
            mutations_foldX = self.prepareMutationFoldX(sequences[0], mutations_pdb)
        else:
            mutations_foldX = self.prepareMutationFoldX(sequences[1], mutations_pdb)

        self.log.debug("mutations_foldX:")
        self.log.debug(mutations_foldX)   

       
        #######################################################################
        ## 2nd: use the 'Repair' feature of FoldX to optimise the structure
        foldX_path = self.tmpPath + self.unique + '/FoldX/'
        try:
            chdir(foldX_path) # from os
            fX = foldX(self.tmpPath + self.unique, 
                               savePath + pdbFile_wt, 
                               chains[0], 
                               self.unique, 
                               self.buildModel_runs,
                               self.foldX_WATER
                               )
            repairedPDB_wt = fX.run('RepairPDB')
        except:
            raise
        finally:
            chdir(self.PWD)
        
        
        #######################################################################
        ## 3rd: introduce the mutation using FoldX
        
        # compile a list of mutations
        mutCodes = list()
        for mut in mutations_foldX:
            mutCodes.append(mut[1])
        
        fX_wt = foldX(self.tmpPath + self.unique, 
                      repairedPDB_wt, 
                      chains[0], 
                      self.unique, 
                      self.buildModel_runs,
                      self.foldX_WATER
                      )
        
        # do the mutation with foldX
        repairedPDB_wt_list, repairedPDB_mut_list = fX_wt.run('BuildModel', mutCodes)
        
        
        #######################################################################
        ## 4th: set up the classes for the wildtype and the mutant structures
        fX_wt_list = list()       
        for wPDB in repairedPDB_wt_list:
            fX_wt_list.append(foldX(self.tmpPath + self.unique, 
                                    wPDB, 
                                    chains[0], 
                                    self.unique, 
                                    self.buildModel_runs,
                                    self.foldX_WATER
                                    )
                               )
            # Copy the foldX wildtype pdb file                          
            shutil.copyfile(wPDB,
                 savePath + str(mutations_modeller[0]) + '-' +  wPDB.split('/')[-1])            
        
        
        fX_mut_list = list()
        for mPDB in repairedPDB_mut_list:
            fX_mut_list.append(foldX(self.tmpPath + self.unique, 
                                     mPDB, 
                                     chains[0], 
                                     self.unique, 
                                     self.buildModel_runs,
                                     self.foldX_WATER
                                     )
                               )
            # Copy the foldX mutant pdb file
            shutil.copyfile(mPDB,
                 savePath + str(mutations_modeller[0]) + '-MUT_' +  mPDB.split('/')[-1])
        
        
        #######################################################################
        ## 5th: calculate the energy for the wildtype
        FoldX_StabilityEnergy_wt = list()
        FoldX_AnalyseComplex_wt  = list()
        for fX_wt in fX_wt_list:
            # stability of the interaction
            if len([ item for item in chains if item != '' ]) == 1:
                AnalyseComplex = [ '0' for i in range(25) ]
            else:
                AnalyseComplex = fX_wt.run('AnalyseComplex')
                AnalyseComplex = [ item for item in AnalyseComplex ] # convert to list
            
            # stability
            StabilityEnergy = fX_wt.run('Stability')
            StabilityEnergy = [ item for item in StabilityEnergy ] # convert to list
            
            # append the results            
            FoldX_StabilityEnergy_wt.append(StabilityEnergy)
            FoldX_AnalyseComplex_wt.append(AnalyseComplex)
        
        
        #######################################################################
        ## 6th: calculate the energy for the mutant
        FoldX_StabilityEnergy_mut = list()
        FoldX_AnalyseComplex_mut = list()
        for fX_mut in fX_mut_list:
            # stability of the interaction
            if len([ item for item in chains if item != '' ]) == 1:
                AnalyseComplex = [ '0' for i in range(25) ]
            else:
                AnalyseComplex = fX_mut.run('AnalyseComplex')
                AnalyseComplex = [ item for item in AnalyseComplex ] # convert to list
            
            # stability
            StabilityEnergy = fX_mut.run('Stability')
            StabilityEnergy = [ item for item in StabilityEnergy ] # convert to list
            
            # append the results            
            FoldX_StabilityEnergy_mut.append(StabilityEnergy)
            FoldX_AnalyseComplex_mut.append(AnalyseComplex)
        
        
        #######################################################################
        ## 7th: combine the energy results
        # AnalyseComplex
        if len(FoldX_AnalyseComplex_wt) == 1:
            AnalyseComplex_energy_wt = FoldX_AnalyseComplex_wt[0]
            AnalyseComplex_energy_mut = FoldX_AnalyseComplex_mut[0]
        else:
            # zip the output for the wildtype
            AnalyseComplex_energy_wt = [ [] for i in range(len(FoldX_AnalyseComplex_wt[0])) ]
            for item in FoldX_AnalyseComplex_wt:
                i = 0
                for element in item:
                    AnalyseComplex_energy_wt[i].append(element)
                    i += 1
            AnalyseComplex_energy_wt = [ ','.join(item) for item in AnalyseComplex_energy_wt ]
            
            # zip the output for the mutant
            AnalyseComplex_energy_mut = [ [] for i in range(len(FoldX_AnalyseComplex_mut[0])) ]
            for item in FoldX_AnalyseComplex_mut:
                i = 0
                for element in item:
                    AnalyseComplex_energy_mut[i].append(element)
                    i += 1
            AnalyseComplex_energy_mut = [ ','.join(item) for item in AnalyseComplex_energy_mut ]

        
        # Stability
        if len(FoldX_StabilityEnergy_wt) == 1:
            Stability_energy_wt = FoldX_StabilityEnergy_wt[0]
            Stability_energy_mut = FoldX_StabilityEnergy_mut[0]
        else:
            # zip the output for the wildtype
            Stability_energy_wt = [ [] for i in range(len(FoldX_StabilityEnergy_wt[0])) ]
            for item in FoldX_StabilityEnergy_wt:
                i = 0
                for element in item:
                    Stability_energy_wt[i].append(element)
                    i += 1
            Stability_energy_wt = [ ','.join(item) for item in Stability_energy_wt ]
            
            # zip the output for the mutant
            Stability_energy_mut = [ [] for i in range(len(FoldX_StabilityEnergy_mut[0])) ]
            for item in FoldX_StabilityEnergy_mut:
                i = 0
                for element in item:
                    Stability_energy_mut[i].append(element)
                    i += 1
            Stability_energy_mut = [ ','.join(item) for item in Stability_energy_mut ]
        
        
        #######################################################################
        ## 8th: calculate the interface size
        ## don't do it if no complex is modelled
        if len([ item for item in chains if item != '' ]) == 1:
            interface_size = ['0', '0', '0']
        else:
            try:
                pops = interfaceSize('', self.tmpPath + self.unique + '/pops/')
                # here the interface between chain 0 and all the other chains is
                # calculated. If more than two chains are present, you can change
                # this behaviour here by adding the addintional chain IDs to the list
                if self.templateFinding:
                    interface_size = pops(repairedPDB_wt_list[0], ['A', ])
                else:
                    chains_complex = [ chain for chain in chains[0] ]
                    interface_size = pops(repairedPDB_wt_list[0], chains_complex)
            except:
                if self.DEBUG:
                    raise
                interface_size = ['0', '0', '0']

        
        #######################################################################
        ## 9th: calculate the pysico-chemical properties
        ## copy the pdb files to the result folder
        get_atomicContactVector = pysiChem(5.0, 4.0, self.tmpPath + self.unique + '/')
        # mutations is of the form I_V70A
        # pysiChem need is like I70 (chain that is mutated and the position)
        # self.mutations is a list which can contain more than one mutations
        mut_physChem = [ mut[0] + mut[3:-1] for mut in mutations_modeller ]
        physChem_mut = list()
        physChem_wt = list()
        physChem_mut_ownChain = list()
        physChem_wt_ownChain = list()
        chdir(self.PWD) # from os
        
        for item in repairedPDB_wt_list:
            # calculate the contact vector
            res_wt = [0, 0, 0, 0]
            res_wt_ownChain = [0, 0, 0, 0]
            i = -1
            for mut in mut_physChem:
                i += 1
                mut_snippet     = mutations_foldX[i][1].split('\n')[0]
                mut_snippet_pos = mutations_foldX[i][0]
                
                chains_complex = [ chain for chain in chains[0] ]
                atomicContactVector, atomicContactVector_ownChain, mutpos_tmp = get_atomicContactVector(item, mut[0], chains_complex, mut_snippet_pos, mut_snippet)
                for index in range(0,4):
                    res_wt[index] += atomicContactVector[index]
                for index in range(0,4):
                    res_wt_ownChain[index] += atomicContactVector_ownChain[index]
            physChem_wt.append(res_wt)
            physChem_wt_ownChain.append(res_wt_ownChain)
        
        for item in repairedPDB_mut_list:
            # calculate the contact vector
            res_mut = [0, 0, 0, 0]
            res_mut_ownChain = [0, 0, 0, 0]
            i = -1
            for mut in mut_physChem:
                i += 1
                mut_snippet     = mutations_foldX[i][1].split('\n')[1]
                mut_snippet_pos = mutations_foldX[i][0]
                
                chains_complex = [ chain for chain in chains[0] ]
                atomicContactVector, atomicContactVector_ownChain, mutpos_tmp = get_atomicContactVector(item, mut[0], chains_complex, mut_snippet_pos, mut_snippet)
                # what ever happens here:
                # atomicContactVector seems not to be a normal list..
                # thus I convert it to one in this ugly way...
                for index in range(0,4):
                    res_mut[index] += atomicContactVector[index]
                for index in range(0,4):
                    res_mut_ownChain[index] += atomicContactVector_ownChain[index]
            physChem_mut.append(res_mut)
            physChem_mut_ownChain.append(res_mut_ownChain)


        DSSP = getDSSP()
        if self.templateFinding:
            mutpos_DSSP = mutpos_tmp
        else:
            mutpos_DSSP = mutations_modeller[0][3:-1]
        try:
            solvent_accessibility_wt, secondary_structure_wt = DSSP(repairedPDB_wt_list[0], mutations_modeller[0][0], mutations_modeller[0][2], mutpos_DSSP)
            solvent_accessibility_mut, secondary_structure_mut = DSSP(repairedPDB_mut_list[0], mutations_modeller[0][0], mutations_modeller[0][-1], mutpos_DSSP)
        except:
            solvent_accessibility_wt, secondary_structure_wt   = '-1', '-1'
            solvent_accessibility_mut, secondary_structure_mut = '-1', '-1'
#            if self.DEBUG:
#                raise
        
        #######################################################################
        ## 10th: cobine the results
        # convert the elements to string and join them
        # they are lists in a list of the form: physChem_mut = [[0, 0, 0, 0], [0, 0, 0, 0]]
        # a little bit confusing but in the following are two list comprehensions
        # it is [item, item] with each item = [element, element, element, element]
        # first convert every element to str
        physChem_mut          = [ [ str(element) for element in item ] for item in physChem_mut ]
        physChem_mut_ownChain = [ [ str(element) for element in item ] for item in physChem_mut_ownChain ]
        physChem_wt           = [ [ str(element) for element in item ] for item in physChem_wt ]
        physChem_wt_ownChain  = [ [ str(element) for element in item ] for item in physChem_wt_ownChain ]
        # then join every 'element' via ':'
        physChem_mut          = [ ':'.join(item) for item in physChem_mut ]
        physChem_mut_ownChain = [ ':'.join(item) for item in physChem_mut_ownChain ]
        physChem_wt           = [ ':'.join(item) for item in physChem_wt ]
        physChem_wt_ownChain  = [ ':'.join(item) for item in physChem_wt_ownChain ]
        
        
        #######################################################################
        ## 11th: get the BLOSUM (or what ever matrix is given) score
        # self.mutations is a list containing the mutations of the form A_H56T
        matrix_score = 0
        for mutation in mutations_modeller:
            fromAA = mutation.split('_')[1][0]
            toAA   = mutation.split('_')[1][-1]
            matrix_score += self.score_pairwise(fromAA, toAA, self.matrix, self.gap_s, self.gap_e)
        
        
        #######################################################################
        ## Compile the data into a dictionary object and return
                
        template_mutations['AnalyseComplex_energy_wt'] = AnalyseComplex_energy_wt
        template_mutations['Stability_energy_wt'] = Stability_energy_wt            
        template_mutations['AnalyseComplex_energy_mut'] = AnalyseComplex_energy_mut
        template_mutations['Stability_energy_mut'] = Stability_energy_mut
        
        template_mutations['interface_size'] = interface_size
        
        template_mutations['physChem_wt'] = physChem_wt
        template_mutations['physChem_wt_ownChain'] = physChem_wt_ownChain            
        template_mutations['physChem_mut'] = physChem_mut
        template_mutations['physChem_mut_ownChain'] = physChem_mut_ownChain

        template_mutations['matrix_score'] = matrix_score
        
        template_mutations['secondary_structure_wt'] = secondary_structure_wt
        template_mutations['solvent_accessibility_wt'] = solvent_accessibility_wt
        template_mutations['secondary_structure_mut'] = secondary_structure_mut
        template_mutations['solvent_accessibility_mut'] = solvent_accessibility_mut
        
        
        #######################################################################
        # Save that data as a json file
        export_data_name = template['saveFolder'] + '_' + ','.join(template_mutations['mutations'])
        with open(template['savePath'] + export_data_name + '.json', 'wb') as fh:
            json.dump(template_mutations, fh)
            
        # Save alignments and modeller models to output database for storage
        # Move template files to the output folder as a tar archive
        self.make_tarfile(self.HOME + self.outputPath + export_data_name + '.tar.bz2', 
                          savePath[:-1])
        
        
        #######################################################################
        self.log.info('Finished processing template:')
        self.log.info(savePath.split('/')[-2])
        self.log.info('\n\n')
                
        return template_mutations

    
    def make_tarfile(self, output_filename, source_dir):
        with tarfile.open(output_filename, "w:bz2") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))
