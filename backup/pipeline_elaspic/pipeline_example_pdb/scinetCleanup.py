# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:04:40 2013

@author: niklas
"""

import subprocess
from os import chdir
from time import strftime

def scinetCleanup(folder, destination, name=None):
    """
    zip and copy the results from the ramdisk to /scratch
    """
    print 'saving the result in', folder
    chdir(folder)
    if name == None:
        output_name = 'result_' + strftime("%Y_%m_%d_at_%Hh_%Mm") + '.tar.bz2'
    else:
        output_name = name + '_' + strftime("%Y_%m_%d_at_%Hh_%Mm") + '.tar.bz2'
    copy = 'cp ' + output_name + ' ' + destination + output_name
#    copy = 'cp ' + output_name + ' $SCRATCH/niklas-pipeline/' + output_name
#    copy = 'cp ' + output_name + ' /home/niklas/tmp/' + output_name
    system_command = 'tar -acf ' + output_name + ' * && ' + copy
    
    childProcess = subprocess.Popen(system_command, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    shell=True, 
                                    )
    result, error = childProcess.communicate()
    if childProcess.returncode != 0:
        print 'error', error
    return