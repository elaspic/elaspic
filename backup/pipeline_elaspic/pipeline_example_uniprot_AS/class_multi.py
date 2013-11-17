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
   
    def __init__(self, proc_name, task_queue, result_queue, runTime, pool, semaphore, DEBUG, outputPath, logger):
        multiprocessing.Process.__init__(self)
        self.proc_name = proc_name
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.DEBUG = DEBUG
        
        self.pool = pool
        self.semaphore = semaphore
        self.outputPath = outputPath
        self.logger = logger
        
        # get the logger from the parent and add a handler
        self.log = logging.getLogger('')
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
                 core_template_database
                 ):
        
        self.templateFinding = templateFinding
        self.uniprotKB = uniprotKB
        self.mutation_uniprot = mutation
        self.mutation_pdb = mutation
        
        
        self.matrix = matrix # which matrix is used for the interface similarity
                             # calculation
        self.gap_s = gap_start  # gap penalty gap start
        self.gap_e = gap_extend # gap penalty gap extend
        
        self.HOME = getcwd() + '/'
        
        self.savePDB = self.HOME + savePDB
        self.tmpPath = tmpPath
        self.outputPath = outputPath
        self.pdbPath = pdbPath
        self.saveAlignments = self.HOME + outputPath + 'alignments/'
        
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
        except errors.KNOTerror:
            if self.DEBUG:
                raise
            else:
                return 'knotted'
        except errors.TcoffeeError, e:
            self.log.error('t_coffee: problem occured while getting the wildtype PDB\n')
            self.log.error('t_coffee: alignment error in file: ' + e.alignInFile + '\n')
            self.log.error('t_coffee: message raised:\n')
            self.log.error('t_coffee:' + e.errors + '\n\n\n')
            if self.DEBUG:
                raise
            else:
                return 'tcoffee_error'
        except ModellerError as e:
            self.log.error('ModellerError while trying to modellel in __getPDB for ' + self.uniprotKB + ':' + self.mutation_uniprot)
            self.log.error('ModellerError args:' + '\n'.join(e.args))
            self.log.error('ModellerError message:' + e.message)
            if self.DEBUG:
                raise
            else:
                return 'modeller_error'
        except errors.FoldXError as e:
            self.log.error('FoldXError while repairing the wildtype for ' + self.uniprotKB + ':' + self.mutation_uniprot)
            self.log.error('FoldXError error:' + e.error)
            if self.DEBUG:
                raise
            else:
                return 'foldx_error'
        except errors.TemplateCoreError as e:
            self.log.error('Error while getting the core template for ' + self.uniprotKB + ':' + self.mutation_uniprot)
            self.log.error('TemplateCoreError error:' + e.error)
            if self.DEBUG:
                raise
            else:
                return 'templateCoreError'
        except errors.TemplateInterfaceError as e:
            self.log.error('Error while getting the interface template for ' + self.uniprotKB + ':' + self.mutation_uniprot)
            self.log.error('TemplateInterfaceError error:' + e.error)
            if self.DEBUG:
                raise
            else:
                return 'templateInterfaceError'
        except errors.pdbError as e:
            self.log.error(e.error)
            if self.DEBUG:
                raise
            else:
                return 'pdbError'
        except:
            if self.DEBUG:
                raise
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
                       chains
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
                                self.savePDB, 
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


    def prepareInput(self, pdbCode, chains, domains, sequences, alignments):

        pdb = pdbTemplate(self.pdbPath, pdbCode, chains, domains, self.savePDB)

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
        self.getTemplateInterface = get_template_interface(self.tmpPath, 
                                                           self.unique, 
                                                           self.pdbPath, 
                                                           self.savePDB, 
                                                           self.saveAlignments,
                                                           self.pool, 
                                                           self.semaphore, 
                                                           self.get_uniprot_sequence, 
                                                           self.interaction_database, 
                                                           self.threeDID_database, 
                                                           self.pdb_resolution_database
                                                           )
        self.getTemplateCore = get_template_core(self.tmpPath, 
                                                 self.unique, 
                                                 self.pdbPath, 
                                                 self.savePDB,
                                                 self.saveAlignments,
                                                 self.pool, 
                                                 self.semaphore, 
                                                 self.get_uniprot_sequence, 
                                                 self.core_template_database
                                                 )
        

                    

        
        # if not, check if it is in the interface and find a template    
        is_in_core = True
        template, new_sequences = self.getTemplateInterface(self.uniprotKB, self.mutation_uniprot)
        if template == []:
            pass
        else:
            is_in_core = False
            pdbCode     = template[0]
            chains      = [template[1], template[2]]
            scores      = [template[3], template[4], template[5]]
            sequences   = [template[9], template[10]] # is a list containing the uniprot sequences cut to the domain
            mutation_position_domain_uniprot = template[13]
            alignments  = [template[11], template[12]]
            mutation    = self.mutation_uniprot[0] + str(template[13]) + self.mutation_uniprot[-1]
            mutations   = [ chains[0] + '_' + mutation, ]
            mutations_uniprot = [chains[0] + '_' + self.mutation_uniprot[0] + str(mutation_position_domain_uniprot) + self.mutation_uniprot[-1], ]
            domains_pdb = [ [int(template[14].split('-')[0]), int(template[14].split('-')[1])], 
                            [int(template[15].split('-')[0]), int(template[15].split('-')[1])]
                          ]

        # check if the mutations falls into the core
        if is_in_core == True:
            template, new_sequences = self.getTemplateCore(self.uniprotKB, self.mutation_uniprot)

            if template == 'not in core' or template == 'no template':
                is_in_core = False
                return 'no template found', self.uniprotKB + '_' + self.mutation_uniprot
            elif template == 'in gap':
                # add the logger here to report that the mutation did fall into a gap
                # in the alignment!
                return 'no template found', self.uniprotKB + '_' + self.mutation_uniprot
            else:
                pdbCode     = template[0]
                chains      = [template[1], ]
                domains_pdb = [template[2], ]
                scores      = [template[3], template[3], 0]
                alignments  = [template[4], ]
                mutation    = self.mutation_uniprot[0] + str(template[7]) + self.mutation_uniprot[-1]
                mutations   = [ chains[0] + '_' + mutation, ]
                sequences   = [template[6], ]
                mutation_position_domain_uniprot = template[7]
                mutations_uniprot = [chains[0] + '_' + self.mutation_uniprot[0] + str(mutation_position_domain_uniprot) + self.mutation_uniprot[-1], ]

        # check if all the templates have 100% sequence identity. If not,
        # set a flag so that the structure is created with modeller
        # this bit has to be implemented properly!!
        do_modelling = True   
#        for score in scores:
#        if float(scores[0]) >= 99.0:
#            do_modelling = False
        
        sequences, alignments, chains, SWITCH_CHAIN, HETflag, HETATMsInChain_SEQnumbering = self.prepareInput(pdbCode, chains, domains_pdb, sequences, alignments)
        
        if SWITCH_CHAIN:
            mutations_foldX = self.prepareMutationFoldX(sequences[1], mutations_uniprot)
        else:
            mutations_foldX = self.prepareMutationFoldX(sequences[0], mutations_uniprot)

        mutations = self.setChainID(chains, mutations, HETflag, SWITCH_CHAIN, do_modelling)
        
       
        # copy the chain IDs
        # they are used afterwards for renaming and copy is needed so that
        # they are not overwritten
        targetIDs   = list()
        templateIDs = list()

        for i in range(0,len(alignments)):
            targetIDs.append(alignments[i][0].id)
            templateIDs.append(alignments[i][1].id)
        

        if do_modelling:
            normDOPE_wt, pdbFile_wt = self.__getModel(alignments, targetIDs, templateIDs, HETATMsInChain_SEQnumbering, chains)
            modeller_path = self.tmpPath + self.unique + '/modeller/'
        else:
            # add a function to take the crytal structure. Compare the __run()
            # method for a possibility. 
            pass
        
        # rename the wildtype pdb file
        pdbFile_wt_renamed = self.uniprotKB + '_' + self.mutation_uniprot + '.pdb'
        system_command = 'mv ' + modeller_path + pdbFile_wt + ' ' + modeller_path + pdbFile_wt_renamed
        subprocess.check_call(system_command, shell=True)
        pdbFile_wt = self.uniprotKB + '_' + self.mutation_uniprot + '.pdb'
        
        return normDOPE_wt, pdbFile_wt_renamed, chains, mutations, modeller_path, mutations_foldX, is_in_core, new_sequences, scores
    
    
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
            result_findTemplate = self.findTemplate()
            if result_findTemplate[0] == 'no template found':
                return result_findTemplate
            else:
                normDOPE_wt, pdbFile_wt, chains, mutations, modeller_path, mutations_foldX, is_in_core, new_sequences, scores = result_findTemplate
                mutation = mutations[0]
                
        else: # get the crystal structure
            modeller_path = self.tmpPath + self.unique + '/'
            is_in_core = True
            new_sequences = []
            scores = ['100', '100', '100']
            # mutations is of the form I_V70A
            pdbCode, chains_complex1, chains_complex2, mutation = self.mutation_pdb.split('_')
            chains = [chains_complex1, chains_complex2]
            chains_get_pdb = [ item for item in chains_complex1 ]
            chains_get_pdb.extend( [ item for item in chains_complex2 ] )
            
            mutations = [mutation[1] + '_' + mutation[0] + mutation[2:], ]
            
            
            sequence, chainNumbering = self.get_pdb_sequence(pdbCode, mutation[1])
            mutation_pos = chainNumbering.index(int(mutation[2:-1]) + 1)  # +1 to have it start with one
            mutation_foldX = mutation[1] + '_' + mutation[0] + str(mutation_pos) + mutation[-1]
            mutations_foldX = self.prepareMutationFoldX(sequence, [mutation_foldX, ])

            
            normDOPE_wt, pdbFile_wt, SWITCH_CHAIN, chains_get_pdb = self.__getCrystalStructure(pdbCode, chains_get_pdb, self.tmpPath + self.unique + '/')
        
        ########################################
        ## 2nd: use the 'Repair' feature of FoldX to optimise the structure
        foldX_path = self.tmpPath + self.unique + '/FoldX/'
        try:
            chdir(foldX_path) # from os
            fX = foldX(self.tmpPath + self.unique, 
                               modeller_path + pdbFile_wt, 
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
        normDOPE_mut = '-'
        mutCodes = list()
        fX_wt_list = [] # flush it
        repairedPDB_wt_list = [] # flush it
    
        for mut in mutations_foldX:
            mutCodes.append(mut[1])

        referenceWT, mutatedPDB = fX_wt.run('BuildModel', mutCodes)

        ########################################
        ## 4th: set up the classes for the wildtype and the mutant structures
        for mPDB in mutatedPDB:
            fX_mut_list.append(foldX(self.tmpPath + self.unique, 
                                     mPDB, 
                                     chains[0], 
                                     self.unique, 
                                     self.buildModel_runs,
                                     self.foldX_WATER
                                     )
                               )
        for wPDB in referenceWT:
            fX_wt_list.append(foldX(self.tmpPath + self.unique, 
                                    wPDB, 
                                    chains[0], 
                                    self.unique, 
                                    self.buildModel_runs,
                                    self.foldX_WATER
                                    )
                               )
                               
        repairedPDB_mut_list = mutatedPDB
        repairedPDB_wt_list = referenceWT
            

            
        
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
        mut_physChem = [ mut[0] + mut[3:-1] for mut in mutations ]
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
            
            # copy the pdb file
#            shutil.copyfile(item, self.HOME + self.outputPath + 'bestModels/Mut_R_' + item.split('/')[-1])
            shutil.copyfile(item, self.HOME + self.outputPath + 'bestModels/Mut_R_' + item.split('/')[-1] + '_' + str(mutation))
        
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
            

            # copy the pdb file
#            shutil.copyfile(item, self.HOME + self.outputPath + 'bestModels/WT_R_' + item.split('/')[-1])
            shutil.copyfile(item, self.HOME + self.outputPath + 'bestModels/WT_R_' + item.split('/')[-1] + '_' + str(mutation))
        

        DSSP = getDSSP()
        if self.templateFinding:
            mutpos_DSSP = mutpos_tmp
        else:
            mutpos_DSSP = mutations[0][3:-1]
        try:
            solvent_accessibility_wt, secondary_structure_wt = DSSP(repairedPDB_wt_list[0], mutations[0][0], mutations[0][2], mutpos_DSSP)
            solvent_accessibility_mut, secondary_structure_mut = DSSP(repairedPDB_mut_list[0], mutations[0][0], mutations[0][-1], mutpos_DSSP)
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
        for mutation in mutations:
            fromAA = mutation.split('_')[1][0]
            toAA   = mutation.split('_')[1][-1]
            matrix_score += self.score_pairwise(fromAA, toAA, self.matrix, self.gap_s, self.gap_e)
            
        if is_in_core:
            is_in_core = 'core'
        else:
            is_in_core = 'interface'
        
        scores = [ str(s) for s in scores ]
        

        return [normDOPE_wt,
                normDOPE_mut,
                AnalyseComplex_energy_mut,
                Stability_energy_mut, 
                AnalyseComplex_energy_wt,
                Stability_energy_wt, 
                interface_size,
                physChem_mut,
                physChem_mut_ownChain,
                physChem_wt,
                physChem_wt_ownChain,
                matrix_score,
                scores,
                new_sequences,
                self.mutation_uniprot,
                self.uniprotKB,
                is_in_core,
                secondary_structure_wt,
                solvent_accessibility_wt,
                secondary_structure_mut,
                solvent_accessibility_mut
                ]
    

    
    
