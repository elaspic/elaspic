#!/usr/bin/env python
"""
Created on Fri Mar  6 18:18:24 2015
@author: Alexey Strokach


In order to run this script from a Spyder console, you first need to import modeller using the
``import modeller`` command. If your paths are set up to use modeller from python 3.4, but
you're actually using python 2.7, you also need to update your ``sys.path`` variable (which is
equivalent to the ``$PYTHONPATH`` variable in bash):

.. code-block: python

    import sys
    sys.path = [
        c if "lib/x86_64-intel8/python3" not in c
        else '/'.join(c.split('/')[:-1]) + "/python2.5"
        for c in sys.path
    ]
    import modeller


"""
#%%
from __future__ import unicode_literals

import os
import argparse
import subprocess
from contextlib import contextmanager

from elaspic import conf


#%% Parse arguments

def get_parser():

    description = """
    Run the ELASPIC pipeline
    """
    
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter)
        
    parser.add_argument(
        '-c', '--config_file', required=True, 
        help='ELASPIC configuration file')
    parser.add_argument(
        '-u', '--uniprot_id', 
        help="The uniprot_id of the protein that you want to mutate (e.g. 'P28223')")
    parser.add_argument(
        '-s', '--pdb_id',
        help="The pdb_id of the structure that you want to mutate (e.g. '4dkl')")
    parser.add_argument(
        '-m', '--mutations', nargs='?', default=['',], 
        help="Mutation(s) that you wish to evaluate (e.g. 'D172E,R173H' or 'A_V10I')")
    parser.add_argument(
        '-p', '--uniprot_domain_pair_ids',  nargs='?', default='',
        help="List of uniprot_domain_pair_ids to analyse "
            "(useful if you want to restrict your analysis to only a handful of domains) ")
    parser.add_argument(
        '-f', '--input_file', 
        help=("A tab separated file of uniprot_ids, mutations, and optionally, uniprot_domain_pair_ids \n"
            "(optional; to be used instead of `--uniprot_id` and `--mutations`)"))
    parser.add_argument(
        '-t', '--run_type', nargs='?', type=int, default=5, choices=[1,2,3,4,5],
        help=('Type of analysis to perform: \n'
            '  1: Calculate Provean only \n'
            '  2: Create homololgy models only \n'
            '  3: Evaluate mutations only \n'
            '  4: Create homology models and evaluate mutations \n'
            '  5: Calculate Provean, create homology models, and evaluate mutations \n'))
    return parser



#%%
def validate_args(args):
    if not os.path.isfile(args.config_file):
        raise Exception('The configuration file {} does not exist!'.format(args.config_file))
        
    if args.input_file and not os.path.isfile(args.input_file):
        raise Exception('The input file {} does not exist!'.format(args.input_file))
        
    if ((args.uniprot_id is None and args.input_file is None) or 
        (args.uniprot_id is not None and args.input_file is not None)):
            raise Exception('One of `--uniprot_id` or `--input_file` must be specified!')


def parse_input_file(input_file):
    """
    Does not work! Do not use!
    """
    uniprot_ids = []
    mutations = []
    uniprot_domain_pair_ids = []
    with open(input_file, 'r') as fh:
        for line in fh:
            # Can skip lines by adding spaces or tabs before them
            if line[0][0] == ' ' or line[0][0] == '\t':
                print('Skipping line: {}'.format(line))
                continue
            row = [ l.strip() for l in line.split('\t') ]
            # Specifying the mutation is optional
            if len(row) == 2:
                uniprot_id, mutation, uniprot_domain_pair_id = row[0], row[1], row[2]
            elif len(row) == 1:
                uniprot_id, mutation, uniprot_domain_pair_id = row[0], row[1], ''
            elif len(row) == 1:
                uniprot_id, mutation, uniprot_domain_pair_id = row[0], '', ''
            uniprot_ids.append(uniprot_id)
            mutations.append(mutation)
            uniprot_domain_pair_ids.append(uniprot_domain_pair_id)
    return uniprot_ids, mutations, uniprot_domain_pair_ids


def elaspic():
    parser = get_parser()
    args = parser.parse_args()
    validate_args(args)
    pipeline = Pipeline(args.config_file)
    if args.input_file:
        uniprot_ids, mutations, uniprot_domain_pair_ids = parse_input_file(args.input_file)
    else:
        uniprot_ids = [args.uniprot_id,]
        mutations = args.mutations if args.mutations else ''
        uniprot_domain_pair_ids = args.uniprot_domain_pair_ids
        
    uniprot_domain_pair_ids_asint = (
        [int(x) for x in uniprot_domain_pair_ids.split(',') if x]
        if uniprot_domain_pair_ids 
        else []
    )
    pipeline(
        uniprot_id=uniprot_ids[0],         
        mutations=mutations, 
        run_type=args.run_type, 
        uniprot_domain_pair_ids=uniprot_domain_pair_ids_asint,
    )

    #~ # Run jobs
    #~ for uniprot_id, mutation in zip(uniprot_ids, mutations):
        #~ print uniprot_id
        #~ print mutation
        #~ print run_type
        #~ uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutation, run_type, n_cores)







#%%
    
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
    args.data_folder = args.data_folder.rstrip('/')
    table_names = args.data_files.split(',') if args.data_files else None
    dirpath, dirnames, filenames = next(os.walk(args.data_folder))
    for table in sql_db.Base.metadata.sorted_tables:
        if table_names is not None and table.name not in table_names:
            print("Skipping table '{}' because it was not included in the 'table_names' list..."
                .format(table.name))
            continue
        if '{}.tsv'.format(table.name) in filenames:
            db.copy_table_to_db(table.name, args.data_folder)
            print("Successfully loaded data from file '{}' to table '{}'"
                .format('{}.tsv'.format(table.name), table.name))
        elif '{}.tsv.gz'.format(table.name) in filenames:
            with open_gzip(os.path.join(args.data_folder, '{}.tsv.gz'.format(table.name))):
                db.copy_table_to_db(table.name, args.data_folder.rstrip('/'))
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
        '--data_folder', default='.', 
        help='Location of text files to be loaded to the database.')
    parser_load_data.add_argument(
        '--data_files', default=None, 
        help=("Names of text files to be loaded to the database. \n"
            "``all`` : load all tables found in the location specified by ``data_folder``."))
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


def elaspic_database():
    parser = get_parser()
    args = parser.parse_args()
    if 'func' not in args.__dict__:
        args = parser.parse_args(['--help'])
    conf.read_configuration_file(args.config_file)
    print("Running function '{}'...".format(args.func.__name__))
    args.func(args)



#%%
if __name__ == '__main__':
    a
    main()

