# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:04:40 2013

@author: niklas
"""
import os
import sys
import shlex
import subprocess
from os import chdir
from time import strftime
from contextlib import contextmanager
import atexit
import signal

###############################################################################
# Used to find the location of the files being executed
def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")


def path_to_pipeline_code():
    encoding = sys.getfilesystemencoding()
    if we_are_frozen():
        return os.path.dirname(unicode(sys.executable, encoding))
    return os.path.dirname(unicode(__file__, encoding))

###############################################################################


def get_username():
    child_process = subprocess.Popen('whoami', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    username, __ = child_process.communicate()
    username  = username.strip()
    return username


def get_hostname():
    child_process = subprocess.Popen('hostname | cut -d. -f1', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    hostname, __ = child_process.communicate()
    hostname = hostname.strip()
    return hostname


def get_echo(system_constant):
    child_process = subprocess.Popen('echo ' + system_constant, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    system_constant_value, __ = child_process.communicate()
    system_constant_value = system_constant_value.strip()
    return system_constant_value


def get_which(bin_name):
    child_process = subprocess.Popen('which ' + bin_name, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    bin_filename, __ = child_process.communicate()
    bin_filename = bin_filename.strip()
    return bin_filename


@contextmanager
def switch_paths(working_path):
    """
    """
    current_path = os.getcwd()
    try:
        os.chdir(working_path)
        yield
    except:
        raise
    finally:
        os.chdir(current_path)


def kill_child_process(child_process):
    if child_process.poll() is not None:
        print ('Child process with pid {} already terminated with return code {}'
                .format(child_process.pid, child_process.returncode))
        return
    try:
        print 'Trying to terminate gracefully child process with pid: {}'.child_process.pid
        os.killpg(child_process.pid, signal.SIGTERM)
#        child_process.terminate()
    except Exception as e:
        print "Didn't work because of error: {}".format(e.message)
        try:
            print 'Trying to kill child process...'
            os.killpg(child_process.pid, signal.SIGKILL)
#            child_process.kill()
        except:
            print "Didn't work because of error: {}".format(e.message)
            print "Letting it go..."
            pass
    print 'OK'


def run_subprocess_locally(working_path, system_command, **popen_argvars):
        with switch_paths(working_path):
            if isinstance(system_command, unicode):
                system_command = system_command.encode('utf8')
            args = shlex.split(system_command)
            child_process = subprocess.Popen(args, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, preexec_fn=os.setpgrp, **popen_argvars)
            atexit(kill_child_process, child_process)
            return child_process


class RunSubprocessLocally(object):

    def __init__(self, working_path, system_command, subprocess_ids, log=None, **popen_argvars):
        self.log = log
        self.subprocess_ids = subprocess_ids
        with switch_paths(working_path):
            if isinstance(system_command, unicode):
                system_command = system_command.encode('utf8')
            args = shlex.split(system_command)
            self.child_process = subprocess.Popen(args, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, preexec_fn=os.setpgrp, **popen_argvars)
        if self.log is not None:
            self.log.debug(
                'Adding subprocess {} with pid {} to the subprocess monitoring list...'
                .format(self.child_process, self.child_process.pid))
        self.subprocess_ids.append(self.child_process.pid)

    def communicate(self, monitor_function=lambda x: None):
        monitor_function(self.child_process)
        result, error_message = self.child_process.communicate()
        returncode = self.child_process.returncode
        if self.log is not None:
            self.log.debug('Removing subprocess {} with pid {} to the subprocess monitoring list...'
            .format(self.child_process, self.child_process.pid))
        self.subprocess_ids.remove(self.child_process.pid)
        return result, error_message, returncode



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


