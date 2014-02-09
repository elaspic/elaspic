# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:09:11 2013

@author: niklas
"""
import os
import logging
import optparse
import multiprocessing

import class_error as error

from class_multi import Consumer, Task
from class_time import runTime as rt


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


class Consumer(multiprocessing.Process):
   
    def __init__(self, proc_name, task_queue, result_queue, runTime, DEBUG, outputPath, webServer=False):
        multiprocessing.Process.__init__(self)
        self.proc_name = proc_name
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.DEBUG = DEBUG
        self.webServer = webServer
        
        self.outputPath = outputPath
        
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



class pipeline_multi():
    
    def run(self):
                    
        ## Part for multiprocessing ##
        # see: http://doughellmann.com/2009/04/pymotw-multiprocessing-part-1.html
        # Establish communication queues
        tasks = multiprocessing.JoinableQueue()
        
        # The main process (pipeline.py) can be collecting results from all subprocesses
        # Not used anymore but may be used in the future
        results = multiprocessing.Queue()
        
        # Start consumers
        print 'Creating %d consumers' % self.num_consumers
        proc_name = [ 'Consumer-' + str(i) for i in range(1, self.num_consumers + 1) ]
        consumers = [ Consumer(proc_name[i-1], tasks, results, self.runTime, self.DEBUG, self.outputPath, webServer=self.webServer)
                      for i in range(1, self.num_consumers + 1) ]

        for w in consumers:
            w.start()
        
        num_jobs = 0
        with open(self.inputFile, 'r') as f:
            for l in f:
                # Can skip lines by adding spaces or tabs before them
                if l[0][0] == ' ' or l[0][0] == '\t':
                    continue
                
                line = [ ll.strip() for ll in l.split('\t') ]
                
                # AS: Mutation does not necessarily have to be specified
                if len(line) > 1:
                    uniprotKB, mutation = line[0], line[1]
                elif len(line) == 1:
                    uniprotKB = line[0]
                    mutation = ''

                print uniprotKB, mutation
                
                # Enqueue jobs                
                num_jobs += 1
                tasks.put( Pipeline(uniprotKB, 
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
                
                log.info('Added to queue uniprot %s with mutation %s' % (uniprotKB, mutation) )
        
        # Add a poison pill for each consumer
        for i in range(1, self.num_consumers+1):
            tasks.put( None )
       
        # Wait for all of the tasks to finish
        tasks.join()


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