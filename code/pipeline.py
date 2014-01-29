# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:09:11 2013

@author: niklas
"""
import os
import subprocess
import logging
import optparse
import multiprocessing

from Bio.SubsMat import MatrixInfo
from ConfigParser import SafeConfigParser

import class_error as error

from class_multi import Consumer, Task
from class_logging import MultiProcessingLog
from scinetCleanup import scinetCleanup

# credit goes to here:
# http://pymotw.com/2/multiprocessing/communication.html#controlling-access-to-resources
class ActivePool(object):
    """
    Used to control how many parallel T_Coffee calls can be run
    Originally implemented because I had problems with running T_Coffee in
    parallel it can now be used if there are memory limitations
    
    The problem I had with T_Coffee is that I made the programm call from the
    same directory. Switching to a unique directory for every call solved
    the issue.
    """
    def __init__(self):
        super(ActivePool, self).__init__()
        self.mgr = multiprocessing.Manager()
        self.active = self.mgr.list()
        self.lock = multiprocessing.Lock()
    def makeActive(self, name):
        with self.lock:
            self.active.append(name)
    def makeInactive(self, name):
        with self.lock:
            self.active.remove(name)
    def __str__(self):
        with self.lock:
            return str(self.active)


class pipeline():
    
    def __init__(self, configFile):
        #######################################################################  
        # read the configuration file and set the variables
        configParser = SafeConfigParser(
            defaults={'tmpPath':'/tmp/pipeline/',
                      'HETATM':True,
                      'DEBUG':False,
                      'saveTo':'$SCRATCH/niklas-pipeline/',
                      'saveScinet':False,
                      'path_to_archive': '/home/kimlab1/database_data/elaspic/human/',
                      'db_type': 'sqlite_file',
                      'db_path': error.path_to_pipeline_code() + '/../db/pipeline.db',
                      'webServer': False,
                      'name': 'pipeline',
                      'numConsumers': multiprocessing.cpu_count(),
                      'tcoffee_parallel_runs': multiprocessing.cpu_count(),
                      'outputPath':'results/',
                      'savePDB': 'results/pdb_files/',
                      'runTime': 'INFINITE'})

        configParser.read(configFile)

        ## from [DEFAULT]
        tmpPath = configParser.get('DEFAULT', 'tmpPath')
        self.DEBUG = configParser.getboolean('DEFAULT', 'DEBUG')
        self.HETATM = configParser.getboolean('DEFAULT', 'HETATM')
        self.saveTo = configParser.get('DEFAULT', 'saveTo')
        self.saveScinet = configParser.getboolean('DEFAULT', 'saveScinet')
        self.path_to_archive = configParser.get('DEFAULT', 'path_to_archive')
        self.db_type = configParser.get('DEFAULT', 'db_type')
        self.db_path = configParser.get('DEFAULT', 'db_path')
        self.webServer = configParser.get('DEFAULT', 'webServer')
            
        ## from [SETTINGS]
        self.name = configParser.get('SETTINGS', 'name')
        self.num_consumers = configParser.getint('SETTINGS', 'numConsumers')
        self.tcoffee_parallel_runs = configParser.getint('SETTINGS', 'tcoffee_parallel_runs')
        self.outputPath = configParser.get('SETTINGS', 'outputPath')
        self.savePDB = configParser.get('SETTINGS', 'savePDB')
        self.runTime = configParser.get('SETTINGS', 'runTime')
        
        # pdbPath
        if configParser.has_option('SETTINGS', 'pdbPath'):
            self.pdbPath = configParser.get('SETTINGS', 'pdbPath')
        else:
            raise error.ConfigError('pdbPath')
        # bin
        if configParser.has_option('SETTINGS', 'bin'):
            self.executables = configParser.get('SETTINGS', 'bin')
        else:
            raise error.ConfigError('bin')
            
        ## from [INPUT]
        # file
        if configParser.has_option('INPUT', 'file'):
            self.inputFile = configParser.get('INPUT', 'file')
            if not os.path.isfile(self.inputFile):
                raise error.DataError(self.inputFile)
        else:
            raise error.ConfigError('file')
        # mutation_uniprot ### should be renamed template_finding
        if configParser.has_option('INPUT', 'mutation_uniprot'):
            self.mutation_uniprot = configParser.getboolean('INPUT', 'mutation_uniprot')
        else:
            raise error.ConfigError('mutation_uniprot')
        
        ## from [MODELLER]
        # modeller_runs
        if configParser.has_option('MODELLER', 'modeller_runs'):
            self.modeller_runs = configParser.getint('MODELLER', 'modeller_runs')
        else:
            raise error.ConfigError('modeller_runs')
        
        ## from [FOLDX]
        # WATER
        if configParser.has_option('FOLDX', 'WATER'):
            self.foldX_WATER = configParser.get('FOLDX', 'WATER')
        else:
            raise error.ConfigError('WATER')
        # buildModel_runs
        if configParser.has_option('FOLDX', 'buildModel_runs'):
            self.buildModel_runs = configParser.get('FOLDX', 'buildModel_runs')
        else:
            raise error.ConfigError('buildModel_runs')
        
        #######################################################################
        # matrix (currently hardcoded, can be changed if needed)
        # in that case also change the gap_start and gap_extend options
        # they are used in conjuntion with the matrix to determine the
        # sequence similariy (could be used to determine the sequence similarity
        # between to interfaces. Might help to improve the modelling)
#        if configParser.has_option('SETTINGS', 'matrix'):
#            self.matrix_option = configParser.get('SETTINGS', 'matrix')
        self.matrix_option = 'blosum80'
        # gap_start
#        if configParser.has_option('SETTINGS', 'gap_start'):
#            self.gap_start = configParser.getint('SETTINGS', 'gap_start')
#        else:
#            raise error.ConfigError('gap_start')
        self.gap_start = -16
        # gap_extend
#        if configParser.has_option('SETTINGS', 'gap_extend'):
#            self.gap_extend = configParser.getint('SETTINGS', 'gap_extend')
#        else:
#            raise error.ConfigError('gap_extend')
        self.gap_extend = -4


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

        # create the tmp directories and copy the binary files
        self.__prepareTMP()
        self.__prepareOutputPaths()

        # set the matrix for the substitution score
        self.matrix = self.__selectMatrix(self.matrix_option)
        
        
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
        if whoami.strip() == 'strokach':
            # when running the pipeline on beagle or banting
            system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                'cp -r /home/kimlab1/strokach/ncbi-blast-2.2.28+/pdbaa_db ' + \
                                self.tmpPath + 'blast/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, __ = childProcess.communicate()
            assert childProcess.returncode == 0
            if childProcess.returncode != 0:
                print 'couldnt copy the blast database!!\n\n\n\n'
        if whoami.strip() == 'alexey':
            # when running the pipeline locally there is no need to copy the database
            # a symlink is enough
            system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                'cd ' + self.tmpPath + 'blast && ' + \
                                'ln -sf /home/kimlab1/strokach/ncbi-blast-2.2.28+/pdbaa_db'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, __ = childProcess.communicate()
            assert childProcess.returncode == 0
        if whoami.strip() == 'joan':
            # for scinet, blast is already installed, but the database needs to be copied
            system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                'cp -ru $HOME/niklas-pipeline/blastdb/pdbaa_db ' + \
                                self.tmpPath + 'blast/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, __ = childProcess.communicate()
            assert childProcess.returncode == 0
            
        
        #######################################################################
        # Copy the database file
        if self.db_type == 'sqlite_file':
            print self.db_path
            print self.tmpPath
            system_command = 'cp -u ' + self.db_path + ' ' + self.tmpPath + '/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, __ = childProcess.communicate()
            assert childProcess.returncode == 0
    
    
    
    def __prepareTMP(self):
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
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i)):
                subprocess.check_call('mkdir ' + self.tmpPath + 'Consumer-' + str(i), shell=True)
                
            # tcoffee
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/tcoffee'):
                # create tmp for tcoffee
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee/tmp && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee/lck && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee/cache'
                subprocess.check_call(mkdir_command, shell=True)
        
            # FoldX
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/FoldX'):
                # make the directories
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/FoldX'
                # copy the executables
                cp_command = 'cp ' + self.executables + 'FoldX.linux64 ' + self.tmpPath + 'Consumer-' + str(i) + '/FoldX/ && ' + \
                           'cp ' + self.executables + 'rotabase.txt ' + self.tmpPath + 'Consumer-' + str(i) + '/FoldX/'
                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
                # Copy dssp into the folder for modelling
                cp_command = 'cp ' + self.executables + 'dssp-2.0.4-linux-amd64 ' + self.tmpPath + 'Consumer-' + str(i) + '/FoldX/'
                subprocess.check_call(cp_command, shell=True)
            
            # modeller
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/modeller'):
                # create workingfolder for modeller
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/modeller'
                subprocess.check_call(mkdir_command, shell=True)
                # Copy knot into the same folder as modeller
                cp_command = 'cp ' + self.executables + 'topol ' + self.tmpPath + 'Consumer-' + str(i) + '/modeller'
                subprocess.check_call(cp_command, shell=True)
                # Copy pops into the folder for modelling
                cp_command = 'cp ' + self.executables + 'pops ' + self.tmpPath + 'Consumer-' + str(i) + '/modeller'
                subprocess.check_call(cp_command, shell=True)
                # Copy dssp into the folder for modelling
                cp_command = 'cp ' + self.executables + 'dssp-2.0.4-linux-amd64 ' + self.tmpPath + 'Consumer-' + str(i) + '/modeller'
                subprocess.check_call(cp_command, shell=True)
                
#            
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
#




    def __prepareOutputPaths(self):
        if not os.path.isdir(self.outputPath):
            subprocess.check_call('mkdir ' + self.outputPath, shell=True)
        # Files are compressed into individual tar archives, do not need these paths anymore
#        if not os.path.isdir(self.outputPath + 'alignments/'):
#            subprocess.check_call('mkdir ' + self.outputPath + 'alignments/', shell=True)
#        if not os.path.isdir(self.outputPath + 'bestModels/'):
#            subprocess.check_call('mkdir ' + self.outputPath + 'bestModels/', shell=True)
#        if not os.path.isdir(self.outputPath + 'pdbFiles/'):
#            subprocess.check_call('mkdir ' +  self.outputPath + 'pdbFiles/', shell=True)

        
    def __selectMatrix(self, matrix):
        if matrix == 'blosum80':
            return MatrixInfo.blosum80
        if matrix == 'blosum60':
            return MatrixInfo.blosum60
        else:
            print 'specified matrix not found!'
        
        

    def run(self):
                    
        ## Part for multiprocessing ##
        # see: http://doughellmann.com/2009/04/pymotw-multiprocessing-part-1.html
        # Establish communication queues
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        
        # create a logger instance
        # I started to play with the logger a little bit but didn't have the time
        # to fully implement it to be really usefull. It works, one just has
        # to set the desired logging with the information where needed
        logger = MultiProcessingLog(self.outputPath + self.name + '.log', mode='w', maxsize=0, rotate=1)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        logger.setFormatter(formatter)
        
        log = logging.getLogger(__name__)
        log.addHandler(logger)
        if self.DEBUG:
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.INFO)
        
        # create the pool to control the number of t_coffee instances
        # that are allowed to run in parallel
        pool = ActivePool()
        s = multiprocessing.Semaphore(self.tcoffee_parallel_runs)
        
        # Start consumers
        print 'Creating %d consumers' % self.num_consumers
        proc_name = [ 'Consumer-' + str(i) for i in range(1, self.num_consumers + 1) ]
        consumers = [ Consumer(proc_name[i-1], tasks, results, self.runTime, pool, s, self.DEBUG, self.outputPath, logger, webServer=self.webServer)
                      for i in range(1, self.num_consumers + 1) ]

        for w in consumers:
            w.start()
    
        # check if a skip file is present
        # I normally didn't use it but the idea is that if you run the pipeline
        # and it reached the maximum runtime, you could easily re-run and skip
        # the already calculated entries.
        if os.path.isfile('processed.log'):
            skipFile = open('processed.log', 'r')
            skip = list()
            for line in skipFile:
                skip.append(line.strip())
            skipFile.close()
            SKIP = True
        else:
            SKIP = False
        
        num_jobs = 0
        with open(self.inputFile, 'r') as f:
            for l in f:
                if l[0][0] == ' ' or l[0][0] == '\t':
                    continue
                
                line = [ ll.strip() for ll in l.split('\t') ]
                
                # AS: Mutation does not necessarily have to be specified
                if len(line) > 1:
                    uniprotKB, mutation = line[0], line[1]
                elif len(line) == 1:
                    uniprotKB = line[0]
                    mutation = ''
                    
                # check if some entries should be skipped
                if SKIP:
                    if uniprotKB + mutation in skip:
                        print 'skipping ', uniprotKB, mutation
                        continue
                
                log.debug('Added to queue uniprot %s with mutation %s' % (uniprotKB, mutation) )
                
                # Enqueue jobs                
                num_jobs += 1
                tasks.put( Task(uniprotKB, 
                                mutation,
                                self.mutation_uniprot,
                                self.savePDB,
                                self.tmpPath,
                                self.outputPath,
                                self.pdbPath,
                                self.matrix,
                                self.gap_start,
                                self.gap_extend,
                                self.modeller_runs,
                                self.buildModel_runs,
                                self.foldX_WATER,
                                self.path_to_archive,
                                self.db_type))
        
        # Add a poison pill for each consumer
        for i in range(1, self.num_consumers+1):
            tasks.put( None )
       
        # Wait for all of the tasks to finish
        tasks.join()

        # process the result
        res_wt = open(self.outputPath + 'result_wt.log', 'w')
        res_mut = open(self.outputPath + 'result_mut.log', 'w')
        res_if = open(self.outputPath + 'result_additional_information.log', 'w')
        skiplog = open(self.outputPath + 'not_processed.log', 'w')
        
        # Write header
        id_labels = 'uniprotIDs\t' + 'pfamIDs\t' + 'domain_defs\t' + 'mutation\t' + 'wt_or_mut\t'
                    
        value_labels = 'normDOPE\t' + \
                    'intraclashes_energy1\t' + 'intraclashes_energy2\t' + \
                    'interaction_energy\t' + 'backbone_hbond\t' + \
                    'sidechain_hbond\t' + 'van_der_waals\t' + \
                    'electrostatics\t' + 'solvation_polar\t' + \
                    'solvation_hydrophobic\t' + 'Van_der_Waals_clashes\t' + \
                    'entropy_sidechain\t' + 'entropy_mainchain\t' + \
                    'sloop_entropy\t' + 'mloop_entropy\t' + 'cis_bond\t' + \
                    'torsional_clash\t' + 'backbone_clash\t' + 'helix_dipole\t' + \
                    'water_bridge\t' + 'disulfide\t' + 'electrostatic_kon\t' + \
                    'partial_covalent_bonds\t' + 'energy_ionisation\t' + \
                    'entropy_complex\t' + 'number_of_residues\t' + \
                    'stability_energy\t' + 'stability_backbone_hbond\t' + \
                    'stability_sidechain_hbond\t' + 'stability_Van_der_Waals\t' + \
                    'stability_electrostatics\t' + 'stability_solvation_polar\t' + \
                    'stability_solvation_hydrophobic\t' + 'stability_Van_der_Waals_clashes\t' + \
                    'stability_entropy_sidechain\t' + 'stability_entropy_mainchain\t' + \
                    'stability_sloop_entropy\t' + 'stability_mloop_entropy\t' + \
                    'stability_cis_bond\t' + 'stability_torsional_clash\t' + \
                    'stability_backbone_clash\t' + 'stability_helix_dipole\t' + \
                    'stability_water_bridge\t' + 'stability_disulfide\t' + \
                    'stability_electrostatic_kon\t' + 'stability_partial_covalent_bonds\t' + \
                    'stability_energy_ionisation\t' + 'stability_entropy_complex\t' + \
                    'stability_number_of_residues\n'
                    
        value_labels_extra = 'core_or_interface\t' + 'seq_id_avg\t' + \
                    'seq_id_chain1\t' + 'seq_id_chain2\t' + \
                    'matrix_score\t' + 'if_hydrophobic\t' + \
                    'if_hydrophilic\t' + 'if_total\t' + \
                    'contactVector_wt_ownChain\t' + 'contactVector_wt\t' + \
                    'contactVector_mut_ownChain\t' + 'contactVector_mut\t' + \
                    'secondary_structure_wt\t' + 'solvent_accessibility_wt\t' + \
                    'secondary_structure_mut\t' + 'solvent_accessibility_mut\n'
        
        res_wt.write(id_labels + value_labels)
        res_mut.write(id_labels + value_labels)
        res_if.write(id_labels + value_labels_extra)

        # Start printing results
        while num_jobs:
            num_jobs -= 1
            # if the timeout was reached and the calculation stopped before
            # every task was calculated the queue is empty before num_jobs is 0
            
            # check if the que is empty
            if results.empty():
                break
            
            output_data = results.get()
            log.debug('output_data: \n%s' % output_data)            
            for output_dict in output_data:
                if isinstance(output_dict, list) and (output_dict[0] == None):
                    skiplog.write(output_dict[1] + '\t' + output_dict[2] + '\n')
                    continue
                elif isinstance(output_dict, list) and (output_dict[0] == 'no template found'):
                    skiplog.write(output_dict[1] + '\t' + 'no template found' + '\n')
                    continue
                else:
                    # Add sequences that were not present already in the database
                    if output_dict.has_key('new_sequences'):
                        self.get_uniprot_sequence.add(output_dict['new_sequences'])
                        
                    # Unique identifier for each line
                    id_data = ('-'.join(output_dict['uniprotIDs']) + '\t' + 
                             '-'.join(output_dict['pfamIDs']) + '\t' + 
                             '_'.join(['-'.join([str(i) for i in x]) for x in output_dict['domain_defs']]) + '\t' + 
                             output_dict['mutation'] + '\t')
#
#                    # Unique identifier for each line
#                    id_data =   output_dict['uniprotIDs'][0] + '\t' + \
#                                output_dict['uniprotIDs'][1] + '\t' + \
#                                output_dict['pfamIDs'][0] + '\t' + \
#                                output_dict['pfamIDs'][1] + '\t' + \
#                                str(output_dict['domain_defs'][0][0]) + '-' + str(output_dict['domain_defs'][0][1]) + '\t' + \
#                                str(output_dict['domain_def2'][1][0]) + '-' + str(output_dict['domain_def2'][1][1]) + '\t' + \
#                                output_dict['mutation'] + '\t'
                                          
                    # Make line for wildtype file
                    resForWrite_wt = (id_data + 'wt\t' +
                                        str(output_dict['normDOPE_wt'])[:6] + '\t' +
                                        '\t'.join(output_dict['AnalyseComplex_energy_wt']) + '\t' +
                                        '\n')
                    
                    # Make line for mutant file
                    resForWrite_mut = (id_data + 'mut\t' + 
                                        '-' + '\t' + # mutant structure has no normDOPE score (same as wildtype)
                                        '\t'.join(output_dict['AnalyseComplex_energy_mut']) + '\t' +
                                        '\t'.join(output_dict['Stability_energy_mut']) +
                                        '\n')
                    
                    # Make line for additional information file
                    resForWrite_if = (id_data + '-\t' +
                                str(output_dict['is_in_core']) + '\t' +
                                '\t'.join([str(i) for i in output_dict['alignment_scores']]) + '\t' +
                                 str(output_dict['matrix_score']) + '\t' +
                                 '\t'.join(output_dict['interface_size']) + '\t' +
                                 ','.join(output_dict['physChem_wt_ownChain']) + '\t' +
                                 ','.join(output_dict['physChem_wt']) + '\t' +
                                 ','.join(output_dict['physChem_mut_ownChain']) + '\t' +
                                 ','.join(output_dict['physChem_mut']) + '\t' +
                                 str(output_dict['secondary_structure_wt'])  + '\t' +
                                 str(output_dict['solvent_accessibility_wt']) + '\t' +
                                 str(output_dict['secondary_structure_mut'])  + '\t' +
                                 str(output_dict['solvent_accessibility_mut']) +
                                 '\n')

                    # Make another file to keep precalculated values?
                    # resForWrite_precalc = []
                    # output_dict['interactingAA']
                    # output_dict['surfaceAA']
                                 
                    # Write output lines             
                    res_wt.write(str(resForWrite_wt))
                    res_mut.write(str(resForWrite_mut))
                    res_if.write(str(resForWrite_if))                 
        
        # save the database (needed if new sequences were added)
        self.get_uniprot_sequence.close()
        
        # close and flush the output files        
        res_wt.close()
        res_mut.close()
        res_if.close()
        skiplog.close()
        logger.close()
        
        # save the results from ramdisk
        if self.saveScinet:
            scinetCleanup(self.outputPath, self.saveTo, self.name)



if __name__ == '__main__':
    # read which configFile to use    
    optParser = optparse.OptionParser()
    optParser.add_option('-c', '--config', action="store", dest='configFile')
    options, args = optParser.parse_args()
       
    configFile = options.configFile
    
    if not os.path.isfile(configFile):
        print 'Error: configFile not found!'
        print 'exiting'
    else:
        try:
            p = pipeline(configFile)
            p.run()
        except error.DataError, e:
            print 'Error: input file',  e.inputFile, ' not found!'
            print 'exiting...'
        except error.ConfigError, e:
            print 'Error: option', e.option, ' not found!'
            print 'exiting...'