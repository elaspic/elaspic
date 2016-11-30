import abc
import functools
import json
import logging
import os

import Bio.Seq
import Bio.SeqRecord
import six

import elaspic

logger = logging.getLogger(__name__)


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


class _Pipeline(abc.ABC):

    _valid_run_types = {
        # Sequence
        'sequence': 'sequence',
        # Model
        'model': 'model',
        # Mutation
        'mutation': 'mutation',
        # Sequence - model
        'sequence,model': 'sequence.model',
        'sequence.model': 'sequence.model',
        # Model - mutation
        'model,mutation': 'model.mutation',
        'model.mutation': 'model.mutation',
        # Sequence - model - mutation
        'sequence,model,mutation': 'sequence.model.mutation',
        'sequence.model.mutation': 'sequence.model.mutation',
        # All
        'all': 'sequence.model.mutation',
    }

    def __init__(self, configurations):
        """.

        It should be possible to initialize one pipeline and call it in parallel using different
        mutations as input.
        """
        # Read the configuration file and set the variables
        if isinstance(configurations, six.string_types):
            elaspic.conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            elaspic.CONFIGS.update(**configurations)

        # Initialize a logger
        for line in ELASPIC_LOGO.split('\n'):
            logger.info(line)

        self.PWD = os.getcwd()

        # Each one leads to the next...
        self.seqrecords = []
        self.sequences = {}
        self.models = {}
        self.predictions = {}

    @abc.abstractmethod
    def run(self, *args, **kwargs):
        raise NotImplementedError

    @staticmethod
    def _split_mutations(mutations):
        if mutations is None:
            return []
        elif not isinstance(mutations, str):
            return mutations
        for sep in [',', ':']:
            if sep in mutations:
                return mutations.split(sep)
        return [mutations]

    @classmethod
    def _validate_run_type(cls, run_type):
        try:
            return cls._valid_run_types[run_type]
        except KeyError:
            raise elaspic.exc.ParameterError("Wrong run_type: '{}'".format(run_type))

    @staticmethod
    def execute_and_remember(f, _instances={}):
        """A basic memoizer.

        .. warning::

            Does not consider ``kwargs``!!!
        """
        @functools.wraps(f)
        def f_new(*args, **kwargs):
            key = tuple([f, *args])
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


# Make Bio objects hashable (hack!)
# TODO: Think more closely about which fields should be used to construct the hash.
Bio.Seq.Seq.__hash__ = lambda self: hash(self.__repr__())
Bio.SeqRecord.SeqRecord.__hash__ = lambda self: hash(self.__repr__())
Bio.SeqRecord.SeqRecord.__eq__ = lambda self, other: self.__hash__() == other.__hash__()


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
