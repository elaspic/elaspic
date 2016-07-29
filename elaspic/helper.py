import os
import sys
import shlex
import shutil
import subprocess
import logging
import json
import string
import functools
from contextlib import contextmanager

logger = logging.getLogger(__name__)


# Logging
class WritableObject:
    """A writable object which writes everything to the logger."""

    def __init__(self, logger):
        self.logger = logger

    def write(self, string):
        self.logger.debug(string.strip())


def slugify(filename_string):
    valid_chars = "-_.()" + string.ascii_letters + string.digits
    return ''.join(c if c in valid_chars else '_' for c in filename_string)


@contextmanager
def log_print_statements(logger):
    """Channel print statements to the debug logger.

    Useful for modules that default to printing things instead of using a logger (Modeller...).
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


# Subprocess
def _set_process_group(parent_process_group_id):
    """Set group_id of the child process to the group_id of the parent process.

    This way when you delete the parent process you also delete all the children.
    """
    child_process_id = os.getpid()
    os.setpgid(child_process_id, parent_process_group_id)


@functools.wraps(subprocess.run)
def run(system_command, **vargs):
    if not isinstance(system_command, (list, tuple)):
        system_command = shlex.split(system_command)
    p = subprocess.run(
        system_command, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        preexec_fn=lambda: _set_process_group(os.getpgrp()),
        **vargs)
    p.stdout = p.stdout.strip()
    p.stderr = p.stderr.strip()
    return p


def get_hostname():
    return run('hostname | cut -d. -f1').stdout


def get_which(bin_name):
    return run('which ' + bin_name).stdout


# Retry
def _check_exception(exc, valid_exc):
    logger.error('The following exception occured:\n{}'.format(exc))
    to_retry = isinstance(exc, valid_exc)
    if to_retry:
        logger.error('Retrying...')
    return to_retry


def retry_database(fn):
    """Decorator to keep probing the database untill you succeed."""
    from retrying import retry
    import sqlalchemy as sa
    r = retry(
        retry_on_exception=lambda exc:
            _check_exception(exc, valid_exc=sa.exc.OperationalError),
        wait_exponential_multiplier=1000,
        wait_exponential_max=60000,
        stop_max_attempt_number=7)
    return r(fn)


def retry_archive(fn):
    """Decorator to keep probing the database untill you succeed."""
    from retrying import retry
    from elaspic import errors
    r = retry(
        retry_on_exception=lambda exc:
            _check_exception(exc, valid_exc=errors.Archive7zipError),
        wait_fixed=2000,
        stop_max_attempt_number=2)
    return r(fn)


# chmod
def copyfile(infile, outfile, mode=None):
    shutil.copyfile(infile, outfile)
    if mode is not None:
        os.chmod(outfile, mode)


def makedirs(path, mode=None, exist_ok=True):
    if mode is None:
        os.makedirs(path, exist_ok=exist_ok)
    else:
        original_umask = os.umask(0)
        try:
            os.makedirs(path, mode=mode, exist_ok=exist_ok)
        finally:
            os.umask(original_umask)


# Locks
def lock(fn):
    """Allow only a single instance of function `fn`, and save results to a lock file."""
    @functools.wraps(fn)
    def locked_fn(self, *args, **kwargs):
        """.

        Returns
        -------
        lock_filename : str
            Lock file that contains function output in json format.

        """
        # Get the lock filename
        if fn.__name__ == 'calculate_provean':
            lock_filename = '{}{}_provean.json'.format(self.pdb_id, args[0])
        elif fn.__name__ == 'calculate_model':
            lock_filename = '{}_modeller.json'.format(self.pdb_id)
        elif fn.__name__ == 'calculate_mutation':
            lock_filename = '{}{}_mutation_{}.json'.format(self.pdb_id, args[0], args[1])
        else:
            raise Exception("Function {} is not supported!".format(fn))

        # Make sure that we can get exclusive rights on the lock
        try:
            lock = open(lock_filename, 'x')
        except FileExistsError:
            try:
                results = json.load(open(lock_filename, 'r'))
                info_message = (
                    "Results have already been calculated and are in file: '{}'.\n"
                    .format(lock_filename, results)
                )
                logger.info(info_message)
                return lock_filename, results
            except ValueError:
                info_message = (
                    "Another process is currently running this function.\n"
                    "If you believe this is an error, delete lock file '{}' and try again."
                    .format(lock_filename)
                )
                logger.info(info_message)
                return lock_filename, None

        # Run the function and write results
        try:
            results = fn(self, *args, **kwargs)
            json.dump(results, lock)
            lock.close()
            return lock_filename, results
        except:
            lock.close()
            os.remove(lock.name)
            raise
    return locked_fn
