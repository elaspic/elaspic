# -*- coding: utf-8 -*-
import os
import json
import six
import logging
from functools import wraps

from . import conf, helper

logger = logging.getLogger(__name__)
configs = conf.Configs()

ELASPIC_LOGO = """

8888888888 888             d8888  .d8888b.  8888888b. 8888888 .d8888b.
888        888            d88888 d88P  Y88b 888   Y88b  888  d88P  Y88b
888        888           d88P888 Y88b.      888    888  888  888    888
8888888    888          d88P 888  "Y888b.   888   d88P  888  888
888        888         d88P  888     "Y88b. 8888888P"   888  888
888        888        d88P   888       "888 888         888  888    888
888        888       d8888888888 Y88b  d88P 888         888  Y88b  d88P
8888888888 88888888 d88P     888  "Y8888P"  888       8888888 "Y8888P"

"""


# %%
class Pipeline:

    def __init__(self, configurations):
        """
        It should be possible to initialize one pipeline and call it in parallel using different
        mutations as input.
        """
        # Read the configuration file and set the variables
        if isinstance(configurations, six.string_types):
            conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            configs.update(**configurations)

        self._validate_temp_path()

        # Initialize a logger
        for line in ELASPIC_LOGO.split('\n'):
            logger.info(line)

        self.PWD = os.getcwd()

        # Each one leads to the next...
        self.seqrecords = []
        self.sequences = {}
        self.models = {}
        self.predictions = {}

    def _validate_temp_path(self):
        """
        Make sure that we are using a job specific temporary folder if we are on a cluster.

        TODO: Remove so error message does not appear in a production release.
        """
        hostname = helper.get_hostname()
        no_job_specific_folder = configs['temp_dir'].startswith(configs['global_temp_dir'])
        on_node_with_manditory_job_specific_folder = (
            any([(x.lower() in hostname) for x in ['node', 'behemoth', 'grendel', 'beagle']])
        )
        if no_job_specific_folder and on_node_with_manditory_job_specific_folder:
            raise Exception(
                'You should be using a temp folder that it specific to the particular job!')


# %%
_instances = {}


def execute_and_remember(f):
    """ Some basic memoizer.
    """
    def f_new(*args, **kwargs):
        key = tuple([f] + list(args))
        if key in _instances:
            return _instances[key].result

        else:
            instance = f(*args, **kwargs)
            if instance:
                with instance:
                    instance.run()
            _instances[key] = instance
            return _instances[key].result

    return f_new


# %%
class Foo:
    def __init__(self):
        self.info = ["I as so great!"]
        print(self.info)

    def __enter__(self):
        self.info.append('I entered')
        print(self.info)

    def run(self):
        self.info.append('And I ran')
        print(self.info)
        raise Exception('Oh oh,,,')

    def __exit__(self, *exc):
        self.info.append('And I exited')
        print(self.info)
        print(exc)
        return True

    def __bool__(self):
        print(self.info)
        return True


# %%
def lock(fn):
    """
    Allow only a single instance of function `fn`,
    and save results to a lock file.
    """
    @wraps(fn)
    def locked_fn(self, *args, **kwargs):
        """

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
