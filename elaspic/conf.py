# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import os.path as op
import subprocess
import tempfile
import logging

from configparser import SafeConfigParser, NoOptionError
from Bio.SubsMat import MatrixInfo

from . import DATA_DIR, helper

logger = logging.getLogger(__name__)


# %%
# Dictionary of configuration values used throughout the ELASPIC pipeline

class Singleton(type):
    instance = None

    def __call__(cls, *args, **kw):
        if not cls.instance:
            cls.instance = super(Singleton, cls).__call__(*args, **kw)
        return cls.instance


# %%
class Configs:
    """A singleton class that keeps track of ELASPIC configuration settings.
    """
    class _Configs:
        pass

    _configs = None

    def __init__(self):
        if Configs._configs is None:
            Configs._configs = Configs._Configs()

    def __getitem__(self, key):
        return getattr(Configs._configs, key)

    def __getattr__(self, key):
        return getattr(Configs._configs, key)

    def __setitem__(self, key, value):
        setattr(Configs._configs, key, value)

    def __setattr__(self, key, value):
        setattr(Configs._configs, key, value)

    def __dir__(self):
        return sorted(set(Configs.__dict__.keys()) | set(Configs._configs.__dict__.keys()))

    def update(self, **kwargs):
        logger.debug(
            'The following configurations will be overwritten: {}'
            .format(set(Configs._configs.__dict__) & set(kwargs))
        )
        Configs._configs.__dict__.update(kwargs)

    def keys(self):
        return Configs._configs.__dict__.keys()

    def values(self):
        return Configs._configs.__dict__.values()

    def items(self):
        return Configs._configs.__dict__.items()

    def clear(self):
        Configs._configs = Configs._Configs

    def get(self, key, fallback=None):
        try:
            return getattr(Configs._configs, key)
        except AttributeError:
            return fallback

    def copy(self):
        return Configs._configs.__dict__.copy()


# %%
configs = Configs()


def read_configuration_file(config_file, unique_temp_dir=None):
    if not os.path.isfile(config_file):
        raise Exception("The configuration file '{}' does not exist!".format(config_file))

    configParser = SafeConfigParser(
        defaults={
            'global_temp_dir': '/tmp',
            'elaspic_foldername': 'elaspic',
            'archive_foldername': 'archive',
            'debug': 'False',
            'look_for_interactions': 'True',
            'remake_provean_supset': 'False',
            'n_cores': '1',
            'web_server': 'False',
            'copy_data': 'True',
            'allow_internet': 'False',
            'testing': 'False',
        })
    configParser.read(config_file)

    # ### [DEFAULT] ###

    # These settings won't change most of the time.
    configs['global_temp_dir'] = configParser.get('DEFAULT', 'global_temp_dir')
    configs['elaspic_foldername'] = configParser.get('DEFAULT', 'elaspic_foldername')
    configs['debug'] = configParser.getboolean('DEFAULT', 'debug')
    configs['look_for_interactions'] = _parse_look_for_interactions(
        configParser.get('DEFAULT', 'look_for_interactions'))
    configs['remake_provean_supset'] = configParser.getboolean('DEFAULT', 'remake_provean_supset')
    configs['n_cores'] = configParser.getint('DEFAULT', 'n_cores')
    configs['web_server'] = configParser.get('DEFAULT', 'web_server')
    configs['copy_data'] = configParser.getboolean('DEFAULT', 'copy_data')
    configs['allow_internet'] = configParser.getboolean(
        'DEFAULT', 'allow_internet')
    configs['testing'] = configParser.getboolean('DEFAULT', 'testing')

    # Temporary directories
    configs['temp_dir'] = get_temp_dir(configs['global_temp_dir'], 'elaspic')
    os.makedirs(configs['temp_dir'], exist_ok=True)
    if unique_temp_dir is not None:
        configs['unique_temp_dir'] = unique_temp_dir
    else:
        configs['unique_temp_dir'] = configParser.get(
            'SETTINGS', 'unique_temp_dir',
            fallback=tempfile.mkdtemp(prefix='', dir=configs['temp_dir'])
        )
    configs['unique'] = op.basename(configs['unique_temp_dir'])
    configs['data_dir'] = configParser.get('SETTINGS', 'data_dir', fallback=DATA_DIR)

    # ### [DATABASE]
    if configParser.has_section('DATABASE'):
        read_database_configs(configParser)

    # ### [SEQUENCE]
    if configParser.has_section('SEQUENCE'):
        read_sequence_configs(configParser)

    # ### [MODEL]
    if configParser.has_section('MODEL'):
        read_model_configs(configParser)

    # TODO: Update `unique_temp_dir` and dependencies.
    _prepare_temp_folders(configs)


def read_database_configs(configParser):
    """[DATABASE]
    """
    # Database
    configs['db_type'] = configParser.get('DATABASE', 'db_type')
    configs['db_is_immutable'] = configParser.get('DATABASE', 'db_is_immutable', fallback=False)
    if configs['db_type'] == 'sqlite':
        configs['sqlite_db_path'] = configParser.get('DATABASE', 'sqlite_db_path')
    elif configs['db_type'] in ['mysql', 'postgresql']:
        configs['db_schema'] = configParser.get('DATABASE', 'db_schema')
        configs['db_schema_uniprot'] = (
            configParser.get('DATABASE', 'db_schema_uniprot', fallback=configs['db_schema']))
        configs['db_database'] = configParser.get('DATABASE', 'db_database', fallback='')
        configs['db_username'] = configParser.get('DATABASE', 'db_username')
        configs['db_password'] = configParser.get('DATABASE', 'db_password')
        configs['db_url'] = configParser.get('DATABASE', 'db_url')
        configs['db_port'] = configParser.get('DATABASE', 'db_port')
        configs['db_socket'] = _get_db_socket(configParser, configs)
    else:
        raise Exception("Only `MySQL`, `PostgreSQL`, and `SQLite` databases are supported!")

    # Archive folder
    configs['archive_type'] = configParser.get('DATABASE', 'archive_type', fallback='directory')
    configs['archive_dir'] = configParser.get('DATABASE', 'archive_dir')
    # supported archive types are 'directory' and 'archive'
    configs['archive_temp_dir'] = op.join(configs['temp_dir'], 'archive')


def read_sequence_configs(configParser):
    """[SEQUENCE]
    """
    configs['sequence_dir'] = configParser.get(
        'SEQUENCE', 'sequence_dir',
        fallback=op.join(configs['unique_temp_dir'], 'sequence')
    )
    configs['provean_temp_dir'] = op.join(configs['sequence_dir'], 'provean_temp')
    _validate_provean_temp_dir(configParser, configs)

    configs['pdb_dir'] = configParser.get('SEQUENCE', 'pdb_dir')
    configs['blast_db_dir'] = configParser.get('SEQUENCE', 'blast_db_dir')
    configs['blast_db_dir_fallback'] = (
        configParser.get('SEQUENCE', 'blast_db_dir_fallback', fallback=''))
    _validate_blast_db_dir(configs)


def read_model_configs(configParser):
    """[MODEL]
    """
    configs['model_dir'] = configParser.get(
        'MODEL', 'model_dir',
        fallback=op.join(configs['unique_temp_dir'], 'model')
    )
    configs['tcoffee_dir'] = op.join(configs['model_dir'], 'tcoffee')

    # Modeller
    configs['modeller_dir'] = op.join(configs['model_dir'], 'modeller')
    configs['modeller_runs'] = configParser.getint('MODEL', 'modeller_runs')

    # FoldX
    configs['foldx_water'] = configParser.get('MODEL', 'foldx_water')
    configs['foldx_num_of_runs'] = configParser.getint('MODEL', 'foldx_num_of_runs')
    configs['matrix_type'] = configParser.get('MODEL', 'matrix_type')
    configs['gap_start'] = configParser.getint('MODEL', 'gap_start')
    configs['gap_extend'] = configParser.getint('MODEL', 'gap_extend')
    configs['matrix'] = getattr(MatrixInfo, configs['matrix_type'])


def _validate_provean_temp_dir(configParser, configs):
    """
    Some nodes on the cluster have a very limited amount of memory for temp storage.
    When working on those nodes, you sould use a remote location for temp storage. This is slow,
    but at least it ensures that you don't crush the nodes by filling up the hard drive.
    However, the most serious problem should be fixed with an up-to-date version of cd-hit.
    (Older versions could go into an infinite loop and generate huge temp files).
    """
    hostname = helper.get_hostname()
    if (('node' in hostname) or ('grendel' in hostname) or ('behemoth' in hostname)):
        try:
            configs['provean_temp_dir'] = configParser.get('SETTINGS', 'provean_temp_dir')
        except configParser.NoOptionError:
            message = (
                "The 'provean_temp_dir' option is required "
                "if you are running on one of hte bc nodes!"
            )
            logger.error(message)
            raise


def _prepare_temp_folders(configs):
    for key, value in configs.items():
        if key.endswith('_dir'):
            logger.debug("Creating '{}' folder: {}...".format(key, value))
            os.makedirs(value, exist_ok=True)


def _validate_blast_db_dir(configs):
    """
    Make sure that configs['blast_db_path'] exists and contains a blast database.
    """
    def blast_db_dir_isvalid(blast_db_dir):
        return op.isdir(blast_db_dir) and op.isfile(op.join(blast_db_dir, 'nr.pal'))

    if blast_db_dir_isvalid(configs['blast_db_dir']):
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


def _get_db_socket(configParser, configs):
    """
    # MySQL: ?unix_socket=/usr/local/mysql5/mysqld.sock
    # PostgreSQL: ?host=/var/lib/postgresql
    """
    socket_prefix = {
        'mysql': '?unix_socket=',
        'postgresql': '?host=',
    }

    if configs['db_url'] == 'localhost':
        try:
            socket_file = configParser.get('DATABASE', 'db_socket')
        except NoOptionError:
            db_socket = ''
        else:
            db_socket = socket_prefix[configs['db_type']] + socket_file
    else:
        db_socket = ''
    return db_socket


def get_temp_dir(global_temp_dir='/tmp', elaspic_foldername=''):
    """ If a :envvar:`TMPDIR` is given as an environment variable, the tmp directory
    is created relative to that. This is useful when running on banting
    (the cluster in the ccbr) and also on Scinet. Make sure that it
    points to '/dev/shm/' on Scinet.
    """
    temp_dir = os.path.join(os.environ.get('TMPDIR', global_temp_dir), elaspic_foldername)
    subprocess.check_call("mkdir -p '{}'".format(temp_dir), shell=True)
    return temp_dir
