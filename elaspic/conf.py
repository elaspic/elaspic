# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:40:02 2015

@author: Alexey Strokach
"""
from __future__ import unicode_literals

import os
import subprocess
from configparser import SafeConfigParser, NoOptionError
from Bio.SubsMat import MatrixInfo


#%%
try:
    code_path = os.path.dirname(os.path.abspath(__file__))
except:
    code_path = os.path.dirname(os.getcwd())


#%%
configs = dict()

def read_configuration_file(config_file):
    
    configParser = SafeConfigParser(
        defaults={
            'global_temp_path': '/tmp/',
            'temp_path_suffix': 'elaspic/',
            'debug': 'False',
            'look_for_interactions': 'True',
            'remake_provean_supset': 'False',
            'n_cores': '1',
            'web_server': 'False',
            'provean_temp_path': '',
        })
    configParser.read(config_file)

    # From [DEFAULT]
    configs['global_temp_path'] = configParser.get('DEFAULT', 'global_temp_path')
    configs['temp_path_suffix'] = configParser.get('DEFAULT', 'temp_path_suffix').strip('/') + '/'
    configs['debug'] = configParser.getboolean('DEFAULT', 'debug')
    configs['look_for_interactions'] = configParser.getboolean('DEFAULT', 'look_for_interactions')
    configs['remake_provean_supset'] = configParser.getboolean('DEFAULT', 'remake_provean_supset')
    configs['n_cores'] = configParser.getint('DEFAULT', 'n_cores')
    configs['web_server'] = configParser.get('DEFAULT', 'web_server')
    configs['provean_temp_path'] = configParser.get('DEFAULT', 'provean_temp_path')

    # From [DATABASE]
    configs['db_type'] = configParser.get('DATABASE', 'db_type')
    configs['db_schema'] = configParser.get('DATABASE', 'db_schema')  
     
    if configs['db_type'] == 'sqlite':
        configs['sqlite_db_path'] = configParser.get('DATABASE', 'sqlite_db_path')
        configs['db_is_immutable'] = True
    elif configs['db_type'] in ['mysql', 'postgresql']:
        try:
            configs['db_database'] = configParser.get('DATABASE', 'db_database')
        except NoOptionError:
            configs['db_database'] = configs['db_schema']
        try:
            configs['db_schema_uniprot'] = configParser.get('DATABASE', 'db_schema_uniprot')
        except NoOptionError:
            configs['db_schema_uniprot'] = configs['db_schema']
        configs['db_username'] = configParser.get('DATABASE', 'db_username')
        configs['db_password'] = configParser.get('DATABASE', 'db_password')
        configs['db_url'] = configParser.get('DATABASE', 'db_url')
        configs['db_port'] = configParser.get('DATABASE', 'db_port')
        configs['db_is_immutable'] = False
    else:
        raise Exception("Only `MySQL`, `PostgreSQL`, and `SQLite` databases are supported!")
        
    # From [SETTINGS]
    configs['path_to_archive'] = configParser.get('SETTINGS', 'path_to_archive')
    configs['blast_db_path'] = configParser.get('SETTINGS', 'blast_db_path')
    try:
        configs['remote_blast_db_path'] = configParser.get('SETTINGS', 'remote_blast_db_path')
    except NoOptionError:
        configs['remote_blast_db_path'] = ''
    configs['pdb_path'] = configParser.get('SETTINGS', 'pdb_path')
    configs['bin_path'] = configParser.get('SETTINGS', 'bin_path')

    # From [GET_MODEL]
    configs['modeller_runs'] = configParser.getint('GET_MODEL', 'modeller_runs')

    # From [GET_MUTATION]
    configs['foldx_water'] = configParser.get('GET_MUTATION', 'foldx_water')
    configs['foldx_num_of_runs'] = configParser.getint('GET_MUTATION', 'foldx_num_of_runs')
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
    
    