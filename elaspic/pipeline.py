import os
import six
import logging
import functools
import Bio.Seq
import Bio.SeqRecord

from . import conf, errors

logger = logging.getLogger(__name__)
configs = conf.CONFIGS

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


class Pipeline:

    _valid_run_types = {
        # Sequence
        '1': 'sequence',
        'sequence': 'sequence',
        # Model
        '2': 'model',
        'model': 'model',
        # Mutation
        '3': 'mutation',
        'mutation': 'mutation',
        # Sequence - model
        '6': 'sequence.model',
        'sequence,model': 'sequence.model',
        'sequence.model': 'sequence.model',
        # Model - mutation
        '4': 'model.mutation',
        'model,mutation': 'model.mutation',
        'model.mutation': 'model.mutation',
        # Sequence - model - mutation
        '5': 'sequence.model.mutation',
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
            conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            configs.update(**configurations)

        # Initialize a logger
        for line in ELASPIC_LOGO.split('\n'):
            logger.info(line)

        self.PWD = os.getcwd()

        # Each one leads to the next...
        self.seqrecords = []
        self.sequences = {}
        self.models = {}
        self.predictions = {}

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
            raise errors.ParameterError("Wrong run_type: '{}'".format(run_type))


# Make Bio objects hashable (hack!)
# TODO: Think more closely about which fields should be used to construct the hash.
Bio.Seq.Seq.__hash__ = lambda self: hash(self.__repr__())
Bio.SeqRecord.SeqRecord.__hash__ = lambda self: hash(self.__repr__())
Bio.SeqRecord.SeqRecord.__eq__ = lambda self, other: self.__hash__() == other.__hash__()


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
