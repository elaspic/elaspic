# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import range
from builtins import object

import os
import sys
import shlex
import subprocess
import signal
import tarfile
import datetime
import logging
import six

from contextlib import contextmanager

from Bio.PDB.PDBParser import PDBParser

if six.PY3:
    from importlib import reload


#%%
canonical_amino_acids = 'ARNDCEQGHILKMFPSTWYV'
uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
logger = None
 

#%%
class WritableObject(object):
    """ A class for collecting all the print statements from modeller in order
    to redirect them to the logger later on
    """
    def __init__(self, logger):
        self.logger = logger
    def write(self, string):
        self.logger.debug(string.strip())


@contextmanager
def log_print_statements(logger):
    """ Channel print statements to the debug logger
    """
    original_stdout = sys.stdout
    original_formatters = []
    for i in range(len(logger.handlers)):
        original_formatters.append(logger.handlers[0].formatter)
        logger.handlers[i].formatter = logging.Formatter('%(message)s')
    wo = WritableObject(logger)
    try:
        sys.stdout = wo
        yield
    except:
        raise
    finally:
        sys.stdout = original_stdout
        for i in range(len(logger.handlers)):
            logger.handlers[i].formatter = original_formatters[i]


#%%
def get_path_to_current_file():
    """ Find the location of the file that is being executed
    """
    encoding = sys.getfilesystemencoding()
    if hasattr(sys, "frozen"):
        # All of the modules are built-in to the interpreter, e.g., by py2exe
        return os.path.dirname(str(sys.executable, encoding))
    else:
        return os.path.dirname(str(__file__, encoding))


def get_logger(do_debug=True, logger_filename=None):
    global logger
    if logger is not None:
        return logger 
    import logging
    reload(logging)
    # Initialize logger
    logger = logging.getLogger(__name__)
    if do_debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.handlers = []
    # Initialize formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # Initialize streamhandler
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    # Initialize filehandler
    if logger_filename is not None:
        fh = logging.FileHandler(logger_filename, mode='w')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    return logger


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:bz2") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


#%% Working with pdb structures
def get_pdb_structure(path_to_pdb_file):
    parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
    structure = parser.get_structure('ID', path_to_pdb_file)
    return structure


#%% Helper functions for sql objects
def decode_domain_def(domains, merge=True, return_string=False):
    """ Unlike split_domain(), this function returns a tuple of tuples of strings,
    preserving letter numbering (e.g. 10B)
    """
    if not domains:
        return None

    if domains[-1] == ',':
        domains = domains[:-1]
    x = domains
    if return_string:
        domain_fragments = [ [r.strip() for r in ro.split(':')] for ro in x.split(',') ]
    else:
        domain_fragments = [ [int(r.strip()) for r in ro.split(':')] for ro in x.split(',') ]
    domain_merged = [ domain_fragments[0][0], domain_fragments[-1][-1] ]
    if merge:
        return domain_merged
    else:
        return domain_fragments


def encode_domain(domains, merged=True):
    x = domains
    if merged:
        return ':'.join([str(r) for r in x])
    else:
        return ','.join([':'.join([str(r) for r in ro]) for ro in x])


def decode_aa_list(interface_aa):
    """
    """
    if interface_aa and (interface_aa != '') and (interface_aa != 'NULL'):
        if interface_aa[-1] == ',':
            interface_aa = interface_aa[:-1]
        x  = interface_aa
        return_tuple = tuple([int(r.strip()) for r in x.split(',')])
    else:
        return_tuple = []
    return return_tuple


def row2dict(row):
    d = {}
    for column in row.__table__.columns:
        d[column.name] = getattr(row, column.name)
        if type(d[column.name]) == datetime.datetime:
            d[column.name] = d[column.name].strftime('%Y-%m-%d %H-%M-%S-%f')
    return d


#%% Helper functions for different subprocess commands

def popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True):
    child_process = subprocess.Popen(system_command, stdout=stdout, stderr=stderr, shell=shell)
    result, error_message = child_process.communicate()
    return_code = child_process.returncode
    if six.PY3:
        result = result.decode()
        error_message = error_message.decode()
    return result, error_message, return_code


def get_username():
    username, __, __ = popen('whoami')
    return username.strip()


def get_hostname():
    hostname, __, __ = popen('hostname | cut -d. -f1')
    return hostname.strip()


def get_echo(system_constant):
    system_constant_value, __, __ = popen('echo ' + system_constant)
    return system_constant_value.strip()


def get_which(bin_name):
    bin_filename, __, __ = popen('which ' + bin_name)
    return bin_filename.strip()


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
        print('Trying to terminate gracefully child process with pid: {}'.child_process.pid)
        os.killpg(child_process.pid, signal.SIGTERM)
#        child_process.terminate()
    except Exception as e:
        print("Didn't work because of error: {}".format(e.__str__()))
        try:
            print('Trying to kill child process...')
            os.killpg(child_process.pid, signal.SIGKILL)
#            child_process.kill()
        except:
            print("Didn't work because of error: {}".format(e.__str__()))
            print("Letting it go...")
            pass
    print('OK')


###############################################################################
# The two functions below can be used to set the subproces group id to the same
# value as the parent process group id. This is a simple way of ensuring that
# all the child processes are terminated when the parent quits, but it makes
# it impossible to terminate the child process group while keeping the parent
# running....
def set_process_group(parent_process_group_id):
    """ This function is used to set the group id of the child process to be
    the same as the group id of the parent process. This way when you delete the
    parent process you also delete all the children.
    """
    child_process_id = os.getpid() #
    os.setpgid(child_process_id, parent_process_group_id)


def run_subprocess_locally(working_path, system_command, **popen_argvars):
    with switch_paths(working_path):
        args = shlex.split(system_command)
        child_process = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            preexec_fn=lambda: set_process_group(os.getpgrp()),
            **popen_argvars)
        return child_process


def run_subprocess_locally_full(working_path, system_command, **popen_argvars):
    child_process = run_subprocess_locally(working_path, system_command, **popen_argvars)
    result, error_message = child_process.communicate()
    if six.PY3:
        result = str(result, encoding='utf-8')
        error_message = str(error_message, encoding='utf-8')
    return_code = child_process.returncode
    return result, error_message, return_code
    
    
    
    
    