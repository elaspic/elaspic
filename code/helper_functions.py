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
        print "Didn't work because of error: {}".format(e.__str__())
        try:
            print 'Trying to kill child process...'
            os.killpg(child_process.pid, signal.SIGKILL)
#            child_process.kill()
        except:
            print "Didn't work because of error: {}".format(e.__str__())
            print "Letting it go..."
            pass
    print 'OK'





###############################################################################
# The two functions below can be used to set the subproces group id to the same
# value as the parent process group id. This is a simple way of ensuring that
# all the child processes are terminated when the parent quits, but it makes
# it impossible to terminate the child process group while keeping the parent
# running....
def set_process_group(parent_process_group_id):
    """ This function is used to set the group id of the child process to be
    the same as the group id of the parent process. This way when you delete the
    parent process you also delete all the children
    """
    child_process_id = os.getpid() #
    os.setpgid(child_process_id, parent_process_group_id)


def run_subprocess_locally(working_path, system_command, **popen_argvars):
    with switch_paths(working_path):
        if isinstance(system_command, unicode):
            system_command = system_command.encode('utf8')
        args = shlex.split(system_command)
        child_process = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            preexec_fn=lambda: set_process_group(os.getpgrp()),
            **popen_argvars)
        return child_process


###############################################################################
#def run_subprocess_locally(
#        working_path, system_command, subprocess_ids, log=None,
#        monitor_function=lambda x: None, **popen_argvars):
#    """
#    """
#    if isinstance(system_command, unicode):
#        system_command = system_command.encode('utf8')
#    args = shlex.split(system_command)
#    with switch_paths(working_path):
#        child_process = subprocess.Popen(
#            args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
#            preexec_fn=os.setpgrp, **popen_argvars)
#    child_process_group_id = os.getpgid(child_process.pid)
#    subprocess_ids.append(child_process_group_id)
#    if log is not None:
#        log.debug(
#            'Adding child process group with id {} to the subprocess monitoring list...'
#            .format(child_process_group_id))
#    try:
#        monitor_function(child_process)
#    except:
#        raise
#    finally:
#        os.killpg(child_process_group_id)
#    result, error_message = child_process.communicate()
#    returncode = child_process.returncode
#    if log is not None:
#        log.debug('Removing child process group with id {} from the subprocess monitoring list...'
#        .format(child_process_group_id))
#    subprocess_ids.remove(child_process.pid)
#    return result, error_message, returncode



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


