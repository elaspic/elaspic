#!/usr/bin/env python
"""
Created on Tue Mar 31 13:00:57 2015
@author: Alexey Strokach
"""
from __future__ import unicode_literals

import os
import argparse
import subprocess
from contextlib import contextmanager

from elaspic import conf
    
def create_database(args):
    from elaspic import sql_db    
    db = sql_db.MyDatabase()
    db.create_database_tables(args.clear_schema, args.keep_uniprot_sequence)
    db.logger.info('Done!')
    

@contextmanager
def open_gzip(filename):
    try:
        print("Gunzipping file '{}'...".format(filename))
        subprocess.check_call("gunzip '{}'".format(filename), shell=True)
    except Exception as e:
        print('Unzipping the file failed with an error: {}'.format(e))
        raise e
    else:
        yield
    finally:
        print("Gzipping the file back again...")
        subprocess.check_call("gzip '{}'".format(filename.rstrip('.gz')), shell=True)
    
    
def load_data_to_database(args):
    from elaspic import sql_db
    db = sql_db.MyDatabase()
    args.table_folder = args.table_folder.rstrip('/')
    table_names = args.table_names.split(',') if args.table_names else None
    dirpath, dirnames, filenames = next(os.walk(args.table_folder))
    for table in sql_db.Base.metadata.sorted_tables:
        if table_names is not None and table.name not in table_names:
            print("Skipping table '{}' because it was not included in the 'table_names' list..."
                .format(table.name))
            continue
        if '{}.tsv'.format(table.name) in filenames:
            db.copy_table_to_db(table.name, args.table_folder)
            print("Successfully loaded data from file '{}' to table '{}'"
                .format('{}.tsv'.format(table.name), table.name))
        elif '{}.tsv.gz'.format(table.name) in filenames:
            with open_gzip(os.path.join(args.table_folder, '{}.tsv.gz'.format(table.name))):
                db.copy_table_to_db(table.name, args.table_folder.rstrip('/'))
            print("Successfully loaded data from file '{}' to table '{}'"
                .format('{}.tsv.gz'.format(table.name), table.name))   


def test_database(args):
    from elaspic import elaspic_testing
    if args.mutation_type in ['domain', 'both']:
        print('*' * 80)
        print('Running a sample domain mutation...')
        test_uniprot_domain = elaspic_testing.TestUniprotDomain()
        test_uniprot_domain.setup_class(uniprot_domain_id=args.uniprot_domain_id)
        test_uniprot_domain.test()
    if args.mutation_type in ['interface', 'both']:   
        print('*' * 80)
        print('Running a sample interface mutation...')
        test_uniprot_domain_pair = elaspic_testing.TestUniprotDomainPair()
        test_uniprot_domain_pair.setup_class(uniprot_domain_pair_id=args.uniprot_domain_pair_id)
        test_uniprot_domain_pair.test()


def delete_database(args):
    from elaspic import sql_db    
    db = sql_db.MyDatabase()
    db.delete_database_tables(args.drop_schema, args.keep_uniprot_sequence)
    db.logger.info('Done!')
    

def get_parser():
    parser = argparse.ArgumentParser(
        prog='ELASPIC database', 
        description="Perform maintenance tasks on the ELASPIC database.")
    parser.add_argument(
        '-c', '--config_file', required=True, 
        help='ELASPIC configuration file')
    subparsers = parser.add_subparsers(
        title='tasks', 
        help='Maintenance tasks to perform')
    
    
    ### Create an empty database schema
    parser_create = subparsers.add_parser(
        name='create',
        description='Create an empty database')
    parser_create.add_argument(
        '--clear_schema', type=bool, default=False, 
        help=('Whether or not to first drop all existing tables from the database schema. \n'
            'WARNING: Choosing `True` will remove all existing data from the schema specified '
            'in your configuration file!!!'))
    parser_create.add_argument(
        '--keep_uniprot_sequence', type=bool, default=True, 
        help=('Whether or not to leave the `uniprot_sequence` table untouched when clearing the schema. \n'
            'Only applicable if ``--clear_schema`` is set to ``True``.'))
    parser_create.set_defaults(func=create_database)
    
    
    ### Load data to the database
    parser_load_data = subparsers.add_parser(
        name='load_data', 
        description='Load data from text files to the database.')
    parser_load_data.add_argument(
        '--table_folder', default='.', 
        help='Location of text files to be loaded to the database.')
    parser_load_data.add_argument(
        '--table_names', default=None, 
        help=("Names of text files to be loaded to the database. \n"
            "``all`` : load all tables found in the location specified by ``table_folder``."))
    parser_load_data.set_defaults(func=load_data_to_database)
    
    
    ### Test the created database by running several mutations
    parser_test = subparsers.add_parser(
        name='test',
        description='Test the database by running some mutations. ')
    parser_test.add_argument(
        '--mutation_type', choices=['domain', 'interface', 'both'], default='both',
        help=("The type of mutatation that you want to test. \n"
            "``domain`` : Test the impact of a mutation on the folding of a domain. \n"
            "``interface`` : Test the impact of a mutation on the interaction between two domains. \n"
            "``both`` : Test both a domain mutation and an interface mutation (DEFAULT) "))
    parser_test.add_argument(
        '--do_provean', type=bool, default=True,
        help="Whether or not to run Provean when testing the database (DEFAULT True)")
    parser_test.add_argument(
        '--uniprot_domain_id', type=int, default=None,
        help="Unique ID identifying the domain that you want to mutate. ")
    parser_test.add_argument(
        '--uniprot_domain_pair_id', type=int, default=True,
        help="Unique ID identifying the domain pair that you want to mutate. ")
    parser_test.set_defaults(func=test_database)
    
    
    ### Delete database
    parser_delete = subparsers.add_parser(
        name='delete',
        description='Delete the database specified in the configuration file.')
    parser_delete.add_argument(
        '--drop_schema', type=bool, default=False, 
        help=('Whether or not to drop the schema that contains the relevant tables. \n'
            'WARNING: Choosing ``True`` will remove all existing data from the schema specified '
            'in your configuration file!!!'))
    parser_delete.add_argument(
        '--keep_uniprot_sequence', type=bool, default=True, 
        help=('Whether or not to leave the ``uniprot_sequence`` table untouched when clearing the schema. \n'
            'Only applicable if ``--drop_schema`` is set to ``True``.'))
    parser_delete.set_defaults(func=delete_database)
    
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    if 'func' not in args.__dict__:
        args = parser.parse_args(['--help'])
    conf.read_configuration_file(args.config_file)
    print("Running function '{}'...".format(args.func.__name__))
    args.func(args)


if __name__ == '__main__':
    main()


