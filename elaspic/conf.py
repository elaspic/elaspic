"""Configure ELASPIC pipeline.

.. todo::

    This is kind of a mess right now...
    We have a function for each section of the config file,
    values are not allowed to be integers, etc.

    Deprecate configuration files and do everything from the command line?
"""
import os
import os.path as op
import tempfile
import logging
import logging.config
import re
import configparser
from kmtools.db_tools import parse_connection_string, make_connection_string
from elaspic import DATA_DIR, helper

logger = logging.getLogger(__name__)

CONFIGS = {}
DEFAULT = {
    'debug': 'False',
    'look_for_interactions': 'True',
    'remake_provean_supset': 'False',
    'n_cores': '1',
    'copy_data': 'True',
    'allow_internet': 'False',
    'testing': 'False',
}


def read_configuration_file(config_file=None, **kwargs):
    if 'DEFAULT' in kwargs:
        DEFAULT.update(kwargs.pop('DEFAULT'))
    config = configparser.ConfigParser(defaults=DEFAULT)
    config_parser('DEFAULT')(config['DEFAULT'])

    if config_file is not None:
        config.read(config_file)

    for category in ['EXTERNAL_DIRS', 'DATABASE', 'MODEL', 'LOGGER']:
        if category not in config:
            config[category] = {}
        if category in kwargs:
            opts = {k: str(v) for k, v in kwargs.pop(category).items()}
            config[category].update(opts)
        config_parser(category)(config[category])

    assert not kwargs

    _prepare_temp_folders(CONFIGS)


def config_parser(category):
    if category == 'DEFAULT':
        return read_default_configs
    elif category == 'EXTERNAL_DIRS':
        return read_sequence_configs
    elif category == 'DATABASE':
        return read_database_configs
    elif category == 'MODEL':
        return read_model_configs
    elif category == 'LOGGER':
        return read_logger_configs
    else:
        raise Exception("Unknown category: '{}'".format(category))


def read_default_configs(config):
    """[DEFAULT]."""
    # These settings won't change most of the time.
    for key in DEFAULT.keys():
        CONFIGS[key] = config.get(key)
    CONFIGS['look_for_interactions'] = (
        _parse_look_for_interactions(CONFIGS['look_for_interactions'])
    )
    CONFIGS['temp_dir'] = get_temp_dir('elaspic')

    # Temporary directories
    CONFIGS['unique_temp_dir'] = config.get(
        'unique_temp_dir',
        fallback=tempfile.mkdtemp(prefix='', dir=CONFIGS['temp_dir'])
    )
    CONFIGS['unique'] = op.basename(CONFIGS['unique_temp_dir'])
    CONFIGS['data_dir'] = config.get('data_dir', fallback=DATA_DIR)


def read_sequence_configs(config):
    """[EXTERNAL_DIRS]."""
    CONFIGS['sequence_dir'] = config.get(
        'sequence_dir',
        fallback=op.join(CONFIGS['unique_temp_dir'], 'sequence')
    )
    CONFIGS['provean_temp_dir'] = op.join(CONFIGS['sequence_dir'], 'provean_temp')
    _validate_provean_temp_dir(config, CONFIGS)

    CONFIGS['pdb_dir'] = config.get('pdb_dir')
    CONFIGS['blast_db_dir'] = config.get('blast_db_dir')
    CONFIGS['blast_db_dir_fallback'] = (
        config.get('blast_db_dir_fallback', fallback=''))
    _validate_blast_db_dir(CONFIGS)

    CONFIGS['archive_dir'] = config.get('archive_dir')
    # Supported archive types are 'directory' and '7zip'
    if CONFIGS['archive_dir'] is None:
        CONFIGS['archive_type'] = None
    elif op.splitext(CONFIGS['archive_dir'])[-1] in ['.7z', '.7zip']:
        assert op.isfile(CONFIGS['archive_dir'])
        CONFIGS['archive_type'] = '7zip'
    else:
        assert op.isdir(CONFIGS['archive_dir'])
        CONFIGS['archive_type'] = 'directory'
    CONFIGS['archive_temp_dir'] = op.join(CONFIGS['temp_dir'], 'archive')


def _validate_provean_temp_dir(config, configs):
    """Some nodes on the cluster have a very limited amount of memory for temp storage.

    When working on those nodes, you sould use a remote location for temp storage. This is slow,
    but at least it ensures that you don't crush the nodes by filling up the hard drive.
    However, the most serious problem should be fixed with an up-to-date version of cd-hit.
    (Older versions could go into an infinite loop and generate huge temp files).
    """
    hostname = helper.get_hostname()
    if (('node' in hostname) or ('grendel' in hostname) or ('behemoth' in hostname)):
        try:
            configs['provean_temp_dir'] = config.get('provean_temp_dir')
        except config.NoOptionError:
            message = (
                "The 'provean_temp_dir' option is required "
                "if you are running on one of hte bc nodes!"
            )
            logger.error(message)
            raise


def read_database_configs(config):
    """[DATABASE]."""
    if config.get('connection_string'):
        CONFIGS['connection_string'] = config.get('connection_string')
        CONFIGS.update(parse_connection_string(CONFIGS['connection_string']))
    elif config.get('db_type'):
        CONFIGS['db_type'] = config.get('db_type')
        CONFIGS['db_schema'] = config.get('db_schema')
        CONFIGS['db_database'] = config.get('db_database', fallback='')
        CONFIGS['db_username'] = config.get('db_username')
        CONFIGS['db_password'] = config.get('db_password')
        CONFIGS['db_url'] = config.get('db_url')
        CONFIGS['db_port'] = config.get('db_port')
        CONFIGS['db_socket'] = _get_db_socket(
            config, CONFIGS['db_type'], CONFIGS['db_url'])
        CONFIGS['connection_string'] = make_connection_string(**CONFIGS)
    CONFIGS['db_is_immutable'] = config.get('db_is_immutable', fallback=False)


def _get_db_socket(config, db_type, db_url):
    """.

    MySQL: ?unix_socket=/usr/local/mysql5/mysqld.sock
    PostgreSQL: ?host=/var/lib/postgresql
    """
    socket_prefix = {
        'mysql': '?unix_socket=',
        'postgresql': '?host=',
    }

    if db_url == 'localhost':
        try:
            socket_file = config.get('db_socket')
        except configparser.NoOptionError:
            db_socket = ''
        else:
            db_socket = socket_prefix[db_type] + socket_file
    else:
        db_socket = ''
    return db_socket


def read_model_configs(config):
    """[MODEL]."""
    CONFIGS['model_dir'] = config.get(
        'model_dir',
        fallback=op.join(CONFIGS['unique_temp_dir'], 'model')
    )
    CONFIGS['tcoffee_dir'] = op.join(CONFIGS['model_dir'], 'tcoffee')

    # Modeller
    CONFIGS['modeller_dir'] = op.join(CONFIGS['model_dir'], 'modeller')
    CONFIGS['modeller_runs'] = config.getint('modeller_runs', 1)

    # FoldX
    CONFIGS['foldx_water'] = config.get('foldx_water', '-IGNORE')
    CONFIGS['foldx_num_of_runs'] = config.getint('foldx_num_of_runs', 1)
    CONFIGS['matrix_type'] = config.get('matrix_type', 'blosum80')
    CONFIGS['gap_start'] = config.getint('gap_start', -16)
    CONFIGS['gap_extend'] = config.getint('gap_extend', -4)


def read_logger_configs(config):
    """Standard logger configuration, with optional tee to a file.

    .. todo:: This needs a cleanup.
    """
    # default_format = '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    default_format = '%(message)s'
    LOGGING_CONFIGS = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'default': {
                'format':
                    config.get('format', fallback=default_format)
            },
        },
        'handlers': {
            'default': {
                'level': config.get('level', 'ERROR'),
                'class': 'logging.StreamHandler',
                'formatter': 'default',
            },
        },
        'loggers': {
            '': {
                'handlers': ['default'],
                'level': 'DEBUG',
                'propagate': True
            }
        }
    }
    logging.config.dictConfig(LOGGING_CONFIGS)


def _prepare_temp_folders(configs):
    for key, value in configs.items():
        if value is None:
            logger.warning("No value provided for key: '{}'".format(key))
            continue
        if key.endswith('_dir') and not re.match('{.*}', value):
            logger.debug("Creating '{}' folder: {}...".format(key, value))
            os.makedirs(value, exist_ok=True)


def _validate_blast_db_dir(configs):
    """Make sure that configs['blast_db_path'] exists and contains a blast database.

    .. todo:: Get rid of 'blast_db_dir_fallback'; it just complicates things.
    """
    def blast_db_dir_isvalid(blast_db_dir):
        return op.isdir(blast_db_dir) and op.isfile(op.join(blast_db_dir, 'nr.pal'))

    if configs['blast_db_dir'] is None or blast_db_dir_isvalid(configs['blast_db_dir']):
        pass
    elif blast_db_dir_isvalid(configs['blast_db_dir_fallback']):
        message = (
            "Using 'blast_db_dir_fallback' because 'blast_db_dir' is not valid!\n"
            "blast_db_dir: {}\n"
            "blast_db_dir_fallback: {}"
            .format(configs['blast_db_dir'], configs['blast_db_dir_fallback'])
        )
        logger.info(message)
        configs['blast_db_dir'] = configs['blast_db_dir_fallback']
    else:
        message = (
            "Both 'blast_db_dir' and 'blast_db_dir_fallback' are not valid!"
            "blast_db_dir: {}\n"
            "blast_db_dir_fallback: {}"
            .format(configs['blast_db_dir'], configs['blast_db_dir_fallback'])
        )
        logger.error(message)


def _parse_look_for_interactions(look_for_interactions):
    if look_for_interactions in ['True', 'False']:
        return int(look_for_interactions == 'True')
    elif look_for_interactions.isnumeric():
        return int(look_for_interactions)
    else:
        raise Exception()


def get_temp_dir(elaspic_temp_dir='elaspic'):
    tempdir = tempfile.gettempdir()
    tempdir = op.join(tempdir, elaspic_temp_dir)
    os.makedirs(tempdir, exist_ok=True)
    tempfile.tempdir = tempdir
    return tempdir
