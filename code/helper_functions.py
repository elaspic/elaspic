# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:04:40 2013

@author: niklas
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from builtins import object
#%%
import os
import sys
import shlex
import subprocess
import signal
import tarfile
import datetime
import logging
import time

from contextlib import contextmanager

from Bio.PDB.PDBParser import PDBParser

from . import sql_db



#%% Elaspic-specific helper functions
def get_uniprot_base_path(d):
    """ The uniprot id is cut into several chunks to create folders that will
    hold a manageable number of pdbs.
    """
    if isinstance(d, sql_db.UniprotDomain):
        uniprot_id = d.uniprot_id
        uniprot_name = d.uniprot_sequence.uniprot_name
    elif isinstance(d, sql_db.UniprotDomainPair):
        uniprot_id = d.uniprot_domain_1.uniprot_id
        uniprot_name = d.uniprot_domain_1.uniprot_sequence.uniprot_name
    elif isinstance(d, dict):
        uniprot_id = d['uniprot_id']
        uniprot_name = d['uniprot_name']
    else:
        raise Exception('Input parameter type is not supported!')

    uniprot_base_path = (
        '{organism_name}/{uniprot_id_part_1}/{uniprot_id_part_2}/{uniprot_id_part_3}/'
        .format(
            organism_name=uniprot_name.split('_')[-1].lower(),
            uniprot_id_part_1=uniprot_id[:3],
            uniprot_id_part_2=uniprot_id[3:5],
            uniprot_id_part_3=uniprot_id,))
    return uniprot_base_path



def get_uniprot_domain_path(d):
    """ Return the path to individual domains or domain pairs.
    """
    if isinstance(d, sql_db.UniprotDomain):
        uniprot_domain_path = (
            '{pfam_clan:.36}.{alignment_def}/'
            .format(
                pfam_clan=d.pfam_clan,
                alignment_def=d.alignment_def.replace(':','-'),))
    elif isinstance(d, sql_db.UniprotDomainPair):
        uniprot_domain_path = (
            '{pfam_clan_1:.36}.{alignment_def_1}/{pfam_clan_2:.36}.{alignment_def_2}/{uniprot_id_2}/'
            .format(
                pfam_clan_1 = d.uniprot_domain_1.pfam_clan,
                alignment_def_1 = d.uniprot_domain_1.alignment_def.replace(':','-'),
                pfam_clan_2 = d.uniprot_domain_2.pfam_clan,
                alignment_def_2 = d.uniprot_domain_2.alignment_def.replace(':','-'),
                uniprot_id_2 = d.uniprot_domain_2.uniprot_id,))
    return uniprot_domain_path



def scinetCleanup(folder, destination, name=None):
    """
    zip and copy the results from the ramdisk to /scratch
    """
    print('saving the result in', folder)
    os.chdir(folder)
    if name == None:
        output_name = 'result_' + time.strftime("%Y_%m_%d_at_%Hh_%Mm") + '.tar.bz2'
    else:
        output_name = name + '_' + time.strftime("%Y_%m_%d_at_%Hh_%Mm") + '.tar.bz2'
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
        print('error', error)
    return


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


def get_logger(do_debug=True):
    import logging
    reload(logging)
    logger = logging.getLogger(__name__)
    if do_debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.handlers = [handler]
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


def get_temp_path(global_temp_path='/tmp', temp_path_suffix=''):
    """ If a TMPDIR is given as an environment variable, the tmp directory
    is created relative to that. This is useful when running on banting
    (the cluster in the ccbr) and also on Scinet. Make sure that it
    points to '/dev/shm/' on Scinet.
    """
    temp_path = os.path.join(os.environ.get('TMPDIR', global_temp_path), temp_path_suffix)
    subprocess.check_call('mkdir -p ' + temp_path, shell=True)
    return temp_path


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
    parent process you also delete all the children
    """
    child_process_id = os.getpid() #
    os.setpgid(child_process_id, parent_process_group_id)


def run_subprocess_locally(working_path, system_command, **popen_argvars):
    with switch_paths(working_path):
        if isinstance(system_command, str):
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
