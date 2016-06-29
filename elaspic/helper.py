import os
import sys
import shlex
import subprocess
import logging
import json
import string
import re
import fcntl
import functools
from contextlib import contextmanager

logger = logging.getLogger(__name__)


@contextmanager
def decompress(file):
    """Temporarly decompress a file."""
    try:
        print("Gunzipping file '{}'...".format(file))
        subprocess.check_call("gunzip '{}'".format(file), shell=True)
    except Exception as e:
        print('Unzipping the file failed with an error: {}'.format(e))
        raise e
    else:
        yield
    finally:
        print("Gzipping the file back again...")
        subprocess.check_call("gzip '{}'".format(file.rstrip('.gz')), shell=True)


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


def decode_domain_def(domains, merge=True, return_string=False):
    """Return a tuple of tuples of strings, preserving letter numbering (e.g. 10B)."""
    if not domains:
        return None, None

    if domains[-1] == ',':
        domains = domains[:-1]
    x = domains
    if return_string:
        domain_fragments = [[r.strip() for r in ro.split(':')] for ro in x.split(',')]
    else:
        domain_fragments = [[int(r.strip()) for r in ro.split(':')] for ro in x.split(',')]
    domain_merged = [domain_fragments[0][0], domain_fragments[-1][-1]]
    if merge:
        return domain_merged
    else:
        return domain_fragments


# Database
def parse_connection_string(connection_string):
    """Split `connection_string` into a dictionary of connection properties.

    Examples
    --------
    >>> from pprint import pprint
    >>> pprint(parse_connection_string('mysql://user:@localhost'))
    {'db_password': '',
     'db_port': '',
     'db_schema': '',
     'db_socket': '',
     'db_type': 'mysql',
     'db_url': 'localhost',
     'db_username': 'user'}
    >>> pprint(parse_connection_string('mysql://user:pass@192.168.0.1:3306/test'))
    {'db_password': 'pass',
     'db_port': '3306',
     'db_schema': 'test',
     'db_socket': '',
     'db_type': 'mysql',
     'db_url': '192.168.0.1',
     'db_username': 'user'}
    >>> pprint(parse_connection_string('sqlite:////absolute/path/to/foo.db'))
    {'db_password': '',
     'db_port': '',
     'db_schema': '/absolute/path/to/foo.db',
     'db_socket': '',
     'db_type': 'sqlite',
     'db_url': '',
     'db_username': ''}
    """
    db_params = {}
    (db_params['db_type'], db_params['db_username'], db_params['db_password'],
     db_params['db_url'], db_params['db_port'], db_params['db_schema'],
     db_params['db_socket']) = (
        re.match(
            '^(\w*)://(|\w*:)(|\w*)(|@localhost|@[0-9\.]*)(|:[0-9]*)(|\/.*)(|\?unix_socket=.*)$',
            connection_string)
        .groups()
    )
    db_params['db_username'] = db_params['db_username'].rstrip(':')
    db_params['db_url'] = db_params['db_url'].lstrip('@')
    db_params['db_port'] = db_params['db_port'].lstrip(':')
    db_params['db_schema'] = (
        db_params['db_schema'][1:]
        if db_params['db_schema'].startswith('/')
        else db_params['db_schema'])
    db_params['db_socket'] = db_params['db_socket'].partition('?unix_socket=')[-1]
    return db_params


def make_connection_string(**vargs):
    """Join a dictionary of connection properties (`vargs`) into a connection string.

    Examples
    --------
    >>> make_connection_string(**{ \
        'db_password': '', \
        'db_port': '', \
        'db_schema': '', \
        'db_socket': '', \
        'db_type': 'mysql', \
        'db_url': 'localhost', \
        'db_username': 'user'})
    'mysql://user:@localhost'
    >>> make_connection_string(**{ \
        'db_password': 'pass', \
        'db_port': '3306', \
        'db_schema': 'test', \
        'db_socket': '', \
        'db_type': 'mysql', \
        'db_url': '192.168.0.1', \
        'db_username': 'user'})
    'mysql://user:pass@192.168.0.1:3306/test'
    >>> make_connection_string(**{ \
        'db_password': '', \
        'db_port': '', \
        'db_schema': '/absolute/path/to/foo.db', \
        'db_socket': '', \
        'db_type': 'sqlite', \
        'db_url': '', \
        'db_username': ''})
    'sqlite:////absolute/path/to/foo.db'
    """
    if vargs['db_username']:
        vargs['db_username'] = vargs['db_username'] + ':'
    if vargs['db_url']:
        vargs['db_url'] = '@' + vargs['db_url']
    if vargs['db_port']:
        assert vargs['db_url']
        vargs['db_port'] = ':' + vargs['db_port']
    if vargs['db_schema']:
        vargs['db_schema'] = '/' + vargs['db_schema']
    if vargs['db_socket']:
        vargs['db_socket'] = '?unix_socket=' + vargs['db_socket']
    connection_string = (
        '{db_type}://{db_username}{db_password}{db_url}{db_port}{db_schema}{db_socket}'
        .format(**vargs)
    )
    return connection_string


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


# Lock
@contextmanager
def open_exclusively(filename, mode='a'):
    fd = os.open(filename, os.O_CREAT | os.O_RDWR)
    fcntl.lockf(fd, fcntl.LOCK_EX)
    try:
        f = os.fdopen(fd, mode)
        yield f
    except:
        raise
    finally:
        f.close()


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
