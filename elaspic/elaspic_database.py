#!/usr/bin/env python
"""
Created on Tue Mar 31 13:00:57 2015
@author: Alexey Strokach
"""
from __future__ import unicode_literals

import os
import six
import random
import argparse

import sqlalchemy as sa

if six.PY3:
    from imp import reload


try:
    base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')
except:
    base_path = os.path.join(os.getcwd(), '../')
code_path = os.path.join(base_path, 'elaspic/')

print('base_path: {}'.format(base_path))
    
    
    
#%%
def create_database(args):
    from elaspic import conf
    from elaspic import sql_db
    conf.read_configuration_file(args.config_file)
    configs = conf.configs.copy()
    my_db = sql_db.MyDatabase(configs)
    my_db.create_database_tables(args.clear_schema, args.keep_uniprot_sequence)
    my_db.logger.info('Done!')
    

def load_data_to_database(args):
    from elaspic import conf
    from elaspic import sql_db
    conf.read_configuration_file(args.config_file)
    configs = conf.configs.copy()
    my_db = sql_db.MyDatabase(configs)


def test_database(args):
    pass


def delete_database(args):
    pass



#%%
def get_parser():
    description = """
    Perform maintenance on the ELASPIC database
    """
    
    parser = argparse.ArgumentParser(description=description)
    subparsers = parser.add_subparsers(title='Database action', help='Database action to perform')
    
    
    ### Create an empty database schema
    parser_create = subparsers.add_parser(
        'create', description='Create an empty database')
    parser_create.add_argument(
        '-c', '--config_file', required=True, 
        help='ELASPIC configuration file')
    parser_create.add_argument(
        '--clear_schema', default=False, 
        help=(
            'Whether or not to clear the schema if it already exists. '
            'WARNING: Choosing `True` will remove all existing data from the schema specified '
            'in your configuration file!!!'))
    parser_create.add_argument(
        '--keep_uniprot_sequence', default=True, 
        help=(
            'Whether or not to leave the `uniprot_sequence` table untouched when clearing the schema.'
            'Only applicable if `--clear_schema` is set to `True`.'))
    parser_create.set_defaults(func=create_database)
    
    
    ### Load data to the database
    parser_load_data = subparsers.add_parser(
        'load_data', description='Load data from text files to the database')
    parser_load_data.add_argument(
        '-c', '--config_file', required=True, 
        help='ELASPIC configuration file')
    parser_load_data.add_argument(
        '--path_to_data', required=True, 
        help='Location of the text files that should be loaded to the database')
    parser_load_data.add_argument(
        '--import_type', choices=['basic', 'complete'], default='complete', 
        help=(
            "Specifies which tables should be loaded. "
            "If `basic`, only the essential tables will be loaded. "
            "(i.e. tables 'domain', 'domain_contact', 'uniprot_domain', 'uniprot_domain_contact', "
            "'uniprot_domain_pair', and 'uniprot_domain_pair_contact'.)"
            "If `complete`, all tables will be loaded."))
    parser_load_data.set_defaults(func=load_data_to_database)
    
    
    ### Test the created database by running several mutations
    parser_test = subparsers.add_parser(
        'test', description='Test the database by running some mutations')
    parser_test.add_argument(
        '-c', '--config_file', required=True, 
        help='ELASPIC configuration file')
    parser_test.add_argument(
        '--recalculate_existing', default=True,
        help='Recalculate existing mutations or calculate new mutations?')
    parser_test.add_argument(
        '--make_provean', default=True,
        help=(
            "Whether or not to run Provean when testing the database. "
            "Only applicable when `--recalculate_existing` is `True`. "))
    parser_test.set_defaults(func=test_database)
    
    
    ### Delete database
    parser_test = subparsers.add_parser(
        'delete', description='Test the database by running some mutations')
    
    return parser


#%%  
def create_clean_schema(configs, clear_schema=True, keep_uniprot_sequence=False, logger=None):
    """Create an empty schema in the database defined by the loaded configuration file
    """
    my_db = MyDatabase(configs, logger)
    clear_schema = True
    keep_uniprot_sequence = False
    
    logger = hf.get_logger()
    result_index = random.randint(0, 100) 
    organism_folder = os.path.join(base_path, 'database/Homo_sapiens_test')


#%% LOAD_DATA

### Load all tables to the database
table_names = [
    'domain', 'domain_contact', 
    'uniprot_sequence', 'provean', 
    'uniprot_domain', 'uniprot_domain_template', 'uniprot_domain_model', 'uniprot_domain_mutation',
    'uniprot_domain_pair', 'uniprot_domain_pair_template', 
    'uniprot_domain_pair_model', 'uniprot_domain_pair_mutation',
]

def load_all_tables_to_db(db_type, configs, table_names):
    for table_name in table_names:
        print(table_name)
        sql_db.load_table_to_db(db_type, configs, table_name)




#%% TEST

def test(DB_TYPE='mysql'):
   
    clear_schema = True
    keep_uniprot_sequence = False
    config_filename = os.path.join(base_path, 'config/', config_filenames[DB_TYPE])
    logger = hf.get_logger()
    result_index = random.randint(0, 100) 
    organism_folder = os.path.join(base_path, 'database/Homo_sapiens_test')
    
    logger.info('*' * 80)
    logger.info('Reading the configuration file...')
    conf.read_configuration_file(config_filename)
    configs = conf.configs.copy()
    configs['organism_folder'] = organism_folder
    configs['result_index'] = result_index
    
    logger.info('*' * 80)
    logger.info('Creating an empty database schema...')
    create_clean_schema(configs, clear_schema, keep_uniprot_sequence, logger)
    
    logger.info('*' * 80)
    logger.info('Loading data from text files to the empty database...')
    load_all_tables_to_db(DB_TYPE, configs)
    
    logger.info('*' * 80)
    logger.info('Running a sample domain mutation...')
    test_uniprot_domain = test_elaspic.TestUniprotDomain()
    test_uniprot_domain.setup_class(configs)
    test_uniprot_domain.test()
    
    logger.info('*' * 80)
    logger.info('Running a sample interface mutation...')
    test_uniprot_domain_pair = test_elaspic.TestUniprotDomainPair()
    test_uniprot_domain_pair.setup_class(configs)
    test_uniprot_domain_pair.test()
    
    logger.info('*' * 80)
    logger.info('Removing the schema that we created for testing')
    delete_schema(configs, logger, clear_schema, keep_uniprot_sequence)



#%% DELETE

def _drop_uniprot_table_and_schema(my_db, configs):
    my_db.engine.execute('drop table {db_schema_uniprot}.uniprot_sequence;'.format(**configs))
    my_db.engine.execute('drop schema {db_schema_uniprot};'.format(**configs))


def delete_schema(configs, logger, clear_schema, keep_uniprot_sequence):
    if not clear_schema:
        return
    
    configs = configs.copy()
    my_db = sql_db.MyDatabase(configs, logger)
    
    for table_name in reversed(table_names):
        if table_name == 'uniprot_sequence':
            continue
        configs['table_name'] = table_name
        my_db.engine.execute('drop table {db_schema}.{table_name};'.format(**configs))
        
    if configs['db_schema'] != configs['db_schema_uniprot']:
        my_db.engine.execute('drop schema {db_schema};'.format(**configs))
        if not keep_uniprot_sequence:
            _drop_uniprot_table_and_schema()
    else:
        if not keep_uniprot_sequence:
            _drop_uniprot_table_and_schema()


def main():
    parser = get_parser()
    args = parser.parse_args()
    print('Running function {}...'.format(args.func.__name__))
    args.func(args)
       




#%%

if __name__ == '__main__':
    main()


