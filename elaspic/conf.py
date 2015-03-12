# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:40:02 2015

@author: Alexey Strokach
"""
from __future__ import unicode_literals

import os
import subprocess
from configparser import SafeConfigParser
from Bio.SubsMat import MatrixInfo

try:
    code_path = os.path.dirname(os.path.abspath(__file__))
except:
    code_path = os.path.dirname(os.path.getcwd())

#SQL_FLAVOUR = 'sqlite_file'
SQL_FLAVOUR = 'mysql'


#%%
configs = dict()

def read_configuration_file(config_file):

    configParser = SafeConfigParser(
        defaults={
            'global_temp_path': '/tmp/',
            'temp_path': 'elaspic/',
            'debug': 'True',
            'look_for_interactions': 'True',
            'remake_provean_supset': 'False',
            'n_cores': '1',
            'schema_version': 'elaspic',
            'web_server': 'False',
            'db_type': SQL_FLAVOUR,
        })
    configParser.read(config_file)

    # From [DEFAULT]
    configs['global_temp_path'] = configParser.get('DEFAULT', 'global_temp_path')
    configs['temp_path_suffix'] = configParser.get('DEFAULT', 'temp_path').strip('/') + '/'
    configs['debug'] = configParser.getboolean('DEFAULT', 'debug')
    configs['look_for_interactions'] = configParser.getboolean('DEFAULT', 'look_for_interactions')
    configs['remake_provean_supset'] = configParser.get('DEFAULT', 'remake_provean_supset')
    configs['n_cores'] = configParser.getint('DEFAULT', 'n_cores')
    configs['schema_version'] = configParser.get('DEFAULT', 'schema_version')
    configs['web_server'] = configParser.get('DEFAULT', 'web_server')
    configs['db_type'] = configParser.get('DEFAULT', 'db_type')
    configs['db_is_immutable'] = True if configs['db_type'].lower().startswith('sqlite') else False

    # From [SETTINGS]
    configs['path_to_archive'] = configParser.get('SETTINGS', 'path_to_archive')
    configs['blast_db_path'] = configParser.get('SETTINGS', 'blast_db_path')
    configs['sqlite_db_path'] = configParser.get('SETTINGS', 'sqlite_db_path')
    configs['pdb_path'] = configParser.get('SETTINGS', 'pdb_path')
    configs['bin_path'] = configParser.get('SETTINGS', 'bin_path')

    # From [GET_MODEL]
    configs['modeller_runs'] = configParser.getint('GET_MODEL', 'modeller_runs')

    # From [GET_MUTATION]
    configs['foldx_water'] = configParser.get('GET_MUTATION', 'foldx_water')
    configs['foldx_num_of_runs'] = configParser.get('GET_MUTATION', 'foldx_num_of_runs')
    configs['matrix_type'] = configParser.get('GET_MUTATION', 'matrix_type')
    configs['gap_start'] = configParser.getint('GET_MUTATION', 'gap_start')
    configs['gap_extend'] = configParser.getint('GET_MUTATION', 'gap_extend')
    configs['matrix'] = getattr(MatrixInfo, configs['matrix_type'])

    configs['temp_path'] = get_temp_path(configs['global_temp_path'], configs['temp_path_suffix'])


def get_temp_path(global_temp_path='/tmp', temp_path_suffix=''):
    """ If a TMPDIR is given as an environment variable, the tmp directory
    is created relative to that. This is useful when running on banting
    (the cluster in the ccbr) and also on Scinet. Make sure that it
    points to '/dev/shm/' on Scinet.
    """
    temp_path = os.path.join(os.environ.get('TMPDIR', global_temp_path), temp_path_suffix)
    subprocess.check_call('mkdir -p ' + temp_path, shell=True)
    return temp_path