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
import cPickle as pickle

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
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'timeout'])
                elif answer == 'tcoffee_error':
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'tcoffee_error'])
                elif answer == 'knotted':
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'knotted'])
                elif answer == 'modeller_error':
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'modeller_error'])
                elif answer == 'templateCoreError':
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'templateCoreError'])
                elif answer == 'templateInterfaceError':
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'templateInterfaceError'])
                elif answer == 'pdbError':
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'pdbError'])
                elif answer == 'foldx_error':
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'foldx_error'])
                elif answer == 'unknown' or answer == None:
                    self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'unknown'])
                else:
                    self.result_queue.put(answer)
            else:
                self.result_queue.put([None, next_task.uniprotKB + '_' + next_task.mutation_uniprot, 'timeout'])
            self.task_queue.task_done()
            

        print 'exiting ', self.proc_name
        return



class Task(object):

    def __init__(self,
                 templateFinding,
                 uniprotKB, 
                 mutation,
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
                 interaction_database,
                 threeDID_database,
                 get_uniprot_sequence,
                 pdb_resolution_database,
                 core_template_database,
                 include_all_pfam_interactions
                 ):
        
        self.templateFinding = templateFinding
        self.uniprotKB = uniprotKB
        self.mutation = mutation
        
        self.include_all_pfam_interactions = include_all_pfam_interactions
        
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
        self.interaction_database    = interaction_database
        self.threeDID_database       = threeDID_database
        self.get_uniprot_sequence    = get_uniprot_sequence
        self.pdb_resolution_database = pdb_resolution_database
        self.core_template_database  = core_template_database
        
        self.pool = None
        self.semaphore = None
        
        self.log = None
        
        self.PWD = None
        
        self.DEBUG = False
        
        

        

    def __call__(self):
        # go to a unique directory
        chdir(self.tmpPath + self.unique)
        self.PWD = getcwd() + '/'
        try:
            return self.__run()
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
            self.log.error('ModellerError while trying to modellel in __getPDB for ' + self.uniprotKB + ':' + self.mutation)
            self.log.error('ModellerError args:' + '\n'.join(e.args))
            self.log.error('ModellerError message:' + e.message)
            if self.DEBUG:
                raise e
            else:
                return 'modeller_error'
        except errors.FoldXError as e:
            self.log.error('FoldXError while repairing the wildtype for ' + self.uniprotKB + ':' + self.mutation)
            self.log.error('FoldXError error:' + e.error)
            if self.DEBUG:
                raise e
            else:
                return 'foldx_error'
        except errors.TemplateCoreError as e:
            self.log.error('Error while getting the core template for ' + self.uniprotKB + ':' + self.mutation)
            self.log.error('TemplateCoreError error:' + e.error)
            if self.DEBUG:
                raise e
            else:
                return 'templateCoreError'
        except errors.TemplateInterfaceError as e:
            self.log.error('Error while getting the interface template for ' + self.uniprotKB + ':' + self.mutation)
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
            chdir(self.PWD)
    
    
    def score_pairwise(self, seq1, seq2, matrix, gap_s, gap_e):
        
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
                       targetIDs, 
                       templateIDs,
                       HETATMsInChain_SEQnumbering,
                       chains,
                       savePDB
                       ):
        outFile = self.tmpPath + self.unique + '/outFile_wildtype'
        
        # generate the input for modeller from the above generated alignment
        prepareModeller(outFile, 
                        alignments, 
                        targetIDs, 
                        templateIDs, 
                        HETATMsInChain_SEQnumbering,
                        chains
                        )
        
        inFile = outFile

        modellerTargetID = '_'.join(targetIDs)
        
        if len(templateIDs) == 1:
            modellerTemplateID = templateIDs[0]
        else:
            modellerTemplateID = templateIDs[0] + templateIDs[1][-1] # if more than two chains are used this has to be adjusted

        modeller_path = self.tmpPath + self.unique + '/modeller/'
        chdir(modeller_path) # from os

        modeller = mod.modeller([inFile], 
                                modellerTargetID, 
                                modellerTemplateID, 
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
        
        self.getTemplateInterface = get_template_interface(self.tmpPath, 
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
                                                           self.include_all_pfam_interactions,
                                                           self.log)
                                                           
        self.getTemplateCore = get_template_core(self.tmpPath, 
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
        templates, new_sequences = self.getTemplateInterface(self.uniprotKB, self.mutation)
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
            templates, new_sequences = self.getTemplateCore(self.uniprotKB, self.mutation)
            self.log.info("Finished getting core templates...\n\n")
#            self.log.debug("Core templates:")
#            self.log.debug(templates)
            
            if templates == 'not in core' or templates == 'no template':
                output_dict = {'errors': 'no template found: ' + str(self.uniprotKB) + '_' + str(self.mutation)}
                output_data.append(output_dict)
            elif templates == 'in gap':
                # add the logger here to report that the mutation did fall into a gap
                # in the alignment!
                output_dict = {'errors': 'no template found: ' + str(self.uniprotKB) + '_' + str(self.mutation)}
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
        template_folder_name = ('_'.join(template['uniprotIDs']) + '_' +
             '_'.join(template['pfamIDs']) + '_' + 
             '_'.join(['-'.join([str(i) for i in x]) for x in template['domain_defs']]) + '_' + 
             self.mutation)
        
        # Folder for storing files for export to output
        savePath =  self.tmpPath + self.unique + '/' + template_folder_name + '/'
        subprocess.check_call('mkdir -p ' + savePath, shell=True)
                    
        # Template data 
        pdbCode = template['pdbID']
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
        sequences, alignments, chains_modeller, SWITCH_CHAIN, HETflag, HETATMsInChain_SEQnumbering = \
            self.prepareInput(pdbCode, chains_pdb, domains_pdb, sequences, alignments, savePath)
        
        self.log.debug("sequences:")
        self.log.debug(sequences)
        self.log.debug("alignments:")
        self.log.debug(alignments)
        self.log.debug("chains_modeller:")
        self.log.debug(chains_modeller)
        self.log.debug("SWITCH_CHAIN:")
        self.log.debug(SWITCH_CHAIN)
        self.log.debug("HETflag:")
        self.log.debug(HETflag)
        self.log.debug("HETATMsInChain_SEQnumbering:")
        self.log.debug(HETATMsInChain_SEQnumbering)
        
        
        if SWITCH_CHAIN:
            mutations_foldX = self.prepareMutationFoldX(sequences[1], mutations_pdb)
        else:
            mutations_foldX = self.prepareMutationFoldX(sequences[0], mutations_pdb)

        self.log.debug("mutations_foldX:")
        self.log.debug(mutations_foldX)
        
        
        # Rename chains in mutation by modeller output
        mutations_modeller = self.setChainID(chains_modeller, mutations_pdb, HETflag, SWITCH_CHAIN, do_modelling)
        
        self.log.debug("mutations_modeller:")
        self.log.debug(mutations_modeller)
        
        
        # copy the chain IDs
        # they are used afterwards for renaming and copy is needed so that
        # they are not overwritten
        targetIDs   = list()
        templateIDs = list()
        for i in range(0,len(alignments)):
            targetIDs.append(alignments[i][0].id)
            templateIDs.append(alignments[i][1].id)
        
        self.log.debug("targetIDs:")
        self.log.debug(targetIDs)
        self.log.debug("templateIDs:")
        self.log.debug(templateIDs)
        
        
        if do_modelling:
            normDOPE_wt, pdbFile_wt = self.__getModel(alignments, targetIDs, templateIDs, HETATMsInChain_SEQnumbering, chains_modeller, savePath)
            modeller_path = self.tmpPath + self.unique + '/modeller/'
        else:
            # add a function to take the crytal structure. Compare the __run()
            # method for a possibility. 
            pass

        self.log.debug("normDOPE_wt:")
        self.log.debug(normDOPE_wt)
        self.log.debug("pdbFile_wt:")
        self.log.debug(pdbFile_wt + '\n\n')
        
        
        # Copy Modeller pdb file
        shutil.copyfile(modeller_path + pdbFile_wt, savePath + pdbFile_wt)
        
        # Rename the wildtype pdb file
#        pdbFile_wt_renamed = self.uniprotKB + '_' + self.mutation + '.pdb'
#        system_command = 'mv ' + modeller_path + pdbFile_wt + ' ' + modeller_path + pdbFile_wt_renamed
#        subprocess.check_call(system_command, shell=True)
        
        # Copy alignments
        for alignment in template['alignments']:
            shutil.copyfile(self.alignments_path + alignment[0].id + '_' + alignment[1].id + '.aln',
                savePath + '/' + template_folder_name + '_' + alignment[0].id + '_' + alignment[1].id + '.aln')
        
        
        # Add all calculated values to the template dictionary
        template['is_in_core'] = is_in_core
        template['template_folder_name'] = template_folder_name
        template['savePath'] = savePath
        template['pdbFile_wt'] = pdbFile_wt
        template['chains_modeller'] = chains_modeller
        template['mutations_pdb'] = mutations_pdb
        template['mutations_modeller'] = mutations_modeller
        template['mutations_foldX'] = mutations_foldX
        template['modeller_path'] = modeller_path
        template['normDOPE_wt'] = normDOPE_wt

        return template


    def findStructure(self):
        # mutations is of the form I_V70A
        pdbCode, chains_pdb1, chains_pdb2, mutation = self.mutation.split('_')
        pdb_type, pdb_resolution = self.pdb_resolution_database(pdbCode)
        
        # Equivalent to class_get_uniprot_template_core_and_interface.py output        
        template = {'uniprotIDs': (self.uniprotKB, ''),
                    'pfamIDs'   : ('', ''),
                    'domain_defs':([], []),
                    'mutation'  : '',
                    'pdbID'     : pdbCode,
                    'chains_pdb' : [chains_pdb1, chains_pdb2],
                    'pdb_domain_defs' : ((), ()),
                    'mutation_position_domain': '',
                    'mutation_pdb': mutation,
                    'pdb_type': pdb_type,
                    'pdb_resolution' : pdb_resolution,
                    'uniprot_domain_sequences' : ('', ''), 
                    'alignments' : ('', ''),
                    'alignment_scores' : (100, 100, 100)}
        
        # Unique template identifier for the template
        template_folder_name = ('_'.join(template['uniprotIDs']) + '_' +
             '_'.join(template['pfamIDs']) + '_' + 
             '_'.join(['-'.join([str(i) for i in x]) for x in template['domain_defs']]) + '_' + 
             self.mutation)
             
        savePath =  self.tmpPath + self.unique + '/' + template_folder_name + '/'
        subprocess.check_call('mkdir -p ' + savePath, shell=True)
        
        modeller_path = self.tmpPath + self.unique + '/'

        chains_get_pdb = [ item for item in chains_pdb1 ]
        chains_get_pdb.extend( [ item for item in chains_pdb2 ] )
        
        mutations_pdb = [mutation[1] + '_' + mutation[0] + mutation[2:], ]
        
        sequence, chainNumbering = self.get_pdb_sequence(pdbCode, mutation[1])
        mutation_pos = chainNumbering.index(int(mutation[2:-1]) + 1)  # +1 to have it start with one
        mutation_foldX = mutation[1] + '_' + mutation[0] + str(mutation_pos) + mutation[-1]
        mutations_foldX = self.prepareMutationFoldX(sequence, [mutation_foldX, ])
        
        normDOPE_wt, pdbFile_wt, SWITCH_CHAIN, chains_get_pdb = self.__getCrystalStructure(pdbCode, chains_get_pdb, modeller_path)

        # Copy Modeller pdb file (could've just replaced the modeller path with savepath above)
        shutil.copyfile(modeller_path + pdbFile_wt, savePath + pdbFile_wt)
        
        # Equivalent to class_multi.py additions        
        template['is_in_core'] = True
        template['template_folder_name'] = template_folder_name
        template['chains_modeller'] = chains_get_pdb
        template['savePath'] = savePath
        template['pdbFile_wt'] = pdbFile_wt
        template['mutations_pdb'] = mutations_pdb
        template['mutations_modeller'] = mutations_pdb
        template['mutations_foldX'] = mutations_foldX
        template['modeller_path'] = modeller_path
        template['normDOPE_wt'] = normDOPE_wt
        
        return [template, ]

    
    def get_pdb_sequence(self, pdbCode, chain):
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


    def __run(self):
        """
        run the pipeline
        
        input:  filePath    path to PDB files for modeller
                tmpPath     path to store some tmp data (not tcoffee)
                unique      set environment for tcoffee (for parallelization)
        """        
        ########################################
        ## 1st: get the pdb structure file
        ## either by modelling or by selecting the pdb crystal structure                     
        if self.templateFinding: # find a template
            templates = self.findTemplate()
        else: # get the structure
            templates = self.findStructure()
            
        # if errors occured, return the errors and don't do further analysis
        if templates == [[]]:
            templates = [{'errors', 'no templates found'}]
            return templates
        elif templates[0].has_key('errors'):
            return templates

        # run analysis for each template
        for template in templates:
            # Unique template identifier for the template
            savePath = template['savePath']     
            pdbFile_wt = template['pdbFile_wt']
            chains = template['chains_modeller']
            mutations_modeller = template['mutations_modeller'] # changed from mutations pdb
            mutations_foldX = template['mutations_foldX']
            
            self.log.debug("savePath:")
            self.log.debug(template['savePath'])            
            self.log.debug("pdbFile_wt:")
            self.log.debug(template['pdbFile_wt'])
            self.log.debug("chains:")
            self.log.debug(template['chains_modeller'])
            self.log.debug("mutations_pdb:")
            self.log.debug(template['mutations_pdb'])
            self.log.debug("mutations_modeller:")
            self.log.debug(template['mutations_modeller'])
            self.log.debug("modeller_path:")
            self.log.debug(template['modeller_path'])
            self.log.debug("mutations_foldx:")
            self.log.debug(template['mutations_foldX'])
                            
            ########################################
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
                repairedPDB_wt_list = [ repairedPDB_wt ]
            except:
                raise
            finally:
                chdir(self.PWD)
            
            
            ########################################
            ## 3rd: introduce the mutation using FoldX
            fX_wt = foldX(self.tmpPath + self.unique, 
                          repairedPDB_wt, 
                          chains[0], 
                          self.unique, 
                          self.buildModel_runs,
                          self.foldX_WATER
                          )
            
            fX_wt_list = [ fX_wt ]
            fX_mut_list = list()
            repairedPDB_mut_list = list()
            
    
            # do the mutation with foldX
            mutCodes = list()
            fX_wt_list = [] # flush it
            repairedPDB_wt_list = [] # flush it
        
            for mut in mutations_foldX:
                mutCodes.append(mut[1])
    
            repairedPDB_wt_list, repairedPDB_mut_list = fX_wt.run('BuildModel', mutCodes)

    
            ########################################
            ## 4th: set up the classes for the wildtype and the mutant structures
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
                
            
            ########################################
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
            
            
            ######################################
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
            
            
            #############################
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
            
            
            #############################
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

            
            #############################
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
            
            #############################
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
            
            
            #############################
            ## 11th: get the BLOSUM (or what ever matrix is given) score
            # self.mutations is a list containing the mutations of the form A_H56T
            matrix_score = 0
            for mutation in mutations_modeller:
                fromAA = mutation.split('_')[1][0]
                toAA   = mutation.split('_')[1][-1]
                matrix_score += self.score_pairwise(fromAA, toAA, self.matrix, self.gap_s, self.gap_e)
            
            template['normDOPE_mut'] = '-'
            template['AnalyseComplex_energy_mut'] = AnalyseComplex_energy_mut
            template['Stability_energy_mut'] = Stability_energy_mut
            template['AnalyseComplex_energy_wt'] = AnalyseComplex_energy_wt
            template['Stability_energy_wt'] = Stability_energy_wt
            template['interface_size'] = interface_size
            template['physChem_mut'] = physChem_mut
            template['physChem_mut_ownChain'] = physChem_mut_ownChain
            template['physChem_wt'] = physChem_wt
            template['physChem_wt_ownChain'] = physChem_wt_ownChain
            template['matrix_score'] = matrix_score
            template['secondary_structure_wt'] = secondary_structure_wt
            template['solvent_accessibility_wt'] = solvent_accessibility_wt
            template['secondary_structure_mut'] = secondary_structure_mut
            template['solvent_accessibility_mut'] = solvent_accessibility_mut

            self.log.info('Finished processing template:')
            self.log.info(savePath.split('/')[-2])
            self.log.info('\n\n')
            
            # Save template dictionary as a pickled object
            filename = (savePath + savePath.split('/')[-2] + '.pickle')
            f = open(filename, 'wb')
            pickle.dump(template, f)
            f.close()
            
            # Move template files to the output folder as a tar archive
            self.make_tarfile(self.HOME + self.outputPath + savePath.split('/')[-2] + '.tar.gz', 
                              savePath[:-1])
            
        self.log.info('Finished processing all templates for ' + self.uniprotKB + ' ' + self.mutation + '\n\n\n')
        return templates
    
    def make_tarfile(self, output_filename, source_dir):
        with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))   
