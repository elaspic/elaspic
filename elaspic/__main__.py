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
# %%
from __future__ import unicode_literals

import os
import argparse
import subprocess
from contextlib import contextmanager
import logging
import logging.config

from elaspic import conf

logger = logging.getLogger(__name__)
configs = conf.Configs()


# %%
def configure_logger():
    level = 'DEBUG' if configs['debug'] else 'INFO'

    LOGGING_CONFIGS = {
        'version': 1,
        'disable_existing_loggers': False,  # this fixes the problem

        'formatters': {
            'standard': {
                'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
            },
            'clean': {
                'format': '%(message)s',
            },
        },
        'handlers': {
            'default': {
                'level': level,
                'class': 'logging.StreamHandler',
                'formatter': 'standard',
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


# %% Parse arguments
def get_elaspic_parser():
    description = """
    Run the ELASPIC pipeline.
    """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        '-c', '--config_file', required=True,
        help='ELASPIC configuration file.')
    parser.add_argument(
        '-u', '--uniprot_id',
        help="The Uniprot ID of the protein that you want to mutate (e.g. 'P28223')."
             "This option relies on a local elaspic database, which has to be specified "
             "in the configuration file.")
    parser.add_argument(
        '-p', '--pdb_file',
        help="Full filename (including path) of the PDB file that you wish to mutate.")
    parser.add_argument(
        '-s', '--sequence_file',
        help="Full filename (including path) of the FASTA file containing the sequence that you "
             "wish to model. If you choose this option, you also have to specify "
             "a template PDB file using the '--pdb-file' option.")
    parser.add_argument(
        '-m', '--mutations', nargs='?', default=[''],
        help="Mutation(s) that you wish to evaluate.\n"
             "If you used '--uniprot_id', mutations must be provided using uniprot coordinates "
             "(e.g. 'D172E,R173H' or 'A_V10I').\n"
             "If you used '--pdb_file', mutations must be provided using the chain "
             "and residue id (e.g. 'A_M1C,B_C20P' to mutate a residue with id '1' on chain A "
             "to Cysteine, and residue with id '20' on chain B to Proline).\n"
             "If you used '--sequence_file', mutations must be provided using the chain "
             "and residue INDEX (e.g. '1_M1C,2_C20P' to mutate the first residue in sequence 1 "
             "to Cysteine, and the 20th residue in sequence 2 to Proline).")
    parser.add_argument(
        '-i', '--uniprot_domain_pair_ids',  nargs='?', default='',
        help="List of uniprot_domain_pair_ids to analyse "
             "(useful if you want to restrict your analysis to only a handful of domains).")
    parser.add_argument(
        '-f', '--input_file',
        help="A tab separated file of uniprot_ids, mutations, "
             "and uniprot_domain_pair_ids (optional).\n"
             "This option can be used instead of `--uniprot_id` and `--mutations` "
             "to input a list of proteins and mutations")
    parser.add_argument(
        '-t', '--run_type', nargs='?', type=int, default=5, choices=[1, 2, 3, 4, 5],
        help=('Type of analysis to perform: \n'
              '  1: Calculate Provean only \n'
              '  2: Create homololgy models only \n'
              '  3: Evaluate mutations only \n'
              '  4: Create homology models and evaluate mutations \n'
              '  5: Calculate Provean, create homology models, and evaluate mutations \n'))
    return parser


def validate_args(args):
    if not os.path.isfile(args.config_file):
        raise Exception('The configuration file {} does not exist!'.format(args.config_file))

    if args.input_file and not os.path.isfile(args.input_file):
        raise Exception('The input file {} does not exist!'.format(args.input_file))

    choose_one = [
        args.uniprot_id is not None, args.input_file is not None, args.pdb_file is not None,
    ]
    if sum(choose_one) != 1:
        raise Exception(
            "One of '--uniprot_id', '--input_file', or '--pdb_file' must be specified!"
        )

    if args.sequence_file and not args.pdb_file:
        raise Exception(
            "A template PDB file must be specified using the '--pdb_file' option, "
            "when you specify a target sequence using the '--sequence_file' option!"
        )


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
            row = [l.strip() for l in line.split('\t')]
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
    parser = get_elaspic_parser()
    args = parser.parse_args()
    validate_args(args)

    if args.input_file:
        conf.read_configuration_file(args.config_file)
        configure_logger()
        # Run database pipeline for each row in file
        from elaspic import database_pipeline
        for uniprot_id, mutations, uniprot_domain_pair_id in \
                zip(*parse_input_file(args.input_file)):
            pipeline = database_pipeline.DatabasePipeline(
                uniprot_id, mutations,
                run_type=args.run_type,
            )
            pipeline.run()
    elif args.uniprot_id:
        conf.read_configuration_file(args.config_file)
        configure_logger()
        # Run database pipeline
        if args.uniprot_domain_pair_ids:
            logger.debug('uniprot_domain_pair_ids: {}'.format(args.uniprot_domain_pair_ids))
            uniprot_domain_pair_ids_asint = (
                [int(x) for x in args.uniprot_domain_pair_ids.split(',') if x]
            )
        else:
            uniprot_domain_pair_ids_asint = []
        # Run database pipeline
        from elaspic import database_pipeline
        pipeline = database_pipeline.DatabasePipeline(
            args.uniprot_id, args.mutations,
            run_type=args.run_type,
            uniprot_domain_pair_ids=uniprot_domain_pair_ids_asint
        )
        pipeline.run()
    elif args.pdb_file:
        conf.read_configuration_file(args.config_file, unique_temp_dir=os.getcwd())
        configure_logger()
        # Run local pipeline
        from elaspic import local_pipeline
        pipeline = local_pipeline.LocalPipeline(
            args.pdb_file, args.sequence_file, args.mutations,
        )
        if args.run_type == 1:
            pipeline.run_all_sequences()
        elif args.run_type == 2:
            pipeline.run_all_models()
        elif args.run_type == 3:
            pipeline.run_all_mutations()
        elif args.run_type == 5:
            pipeline.run()


# %%
def create_database(args):
    from elaspic import database
    db = database.MyDatabase()
    db.create_database_tables(args.clear_schema, args.keep_uniprot_sequence)
    logger.info('Done!')


@contextmanager
def open_gzip(filename):
    """
    Temporarly unzip a file so that it can be processed with tools that do not work with
    compressed archives.
    """
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
    from elaspic import database
    db = database.MyDatabase()
    args.data_folder = args.data_folder.rstrip('/')
    table_names = args.data_files.split(',') if args.data_files else None
    dirpath, dirnames, filenames = next(os.walk(args.data_folder))
    for table in database.Base.metadata.sorted_tables:
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
    raise NotImplementedError


def delete_database(args):
    from elaspic import database
    db = database.MyDatabase()
    db.delete_database_tables(args.drop_schema, args.keep_uniprot_sequence)
    logger.info('Done!')


def get_database_parser():
    parser = argparse.ArgumentParser(
        prog='ELASPIC database',
        description="Perform maintenance tasks on the ELASPIC database.")
    parser.add_argument(
        '-c', '--config_file', required=True,
        help='ELASPIC configuration file')
    subparsers = parser.add_subparsers(
        title='tasks',
        help='Maintenance tasks to perform')

    # Create an empty database schema
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
        help="Whether or not to leave the 'uniprot_sequence' table untouched when clearing "
             "the schema. Only applicable if '--clear_schema' is set to 'True'.")
    parser_create.set_defaults(func=create_database)

    # Load data to the database
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

    # Test the created database by running several mutations
    parser_test = subparsers.add_parser(
        name='test',
        description='Test the database by running some mutations. ')
    parser_test.add_argument(
        '--mutation_type', choices=['domain', 'interface', 'both'], default='both',
        help="The type of mutatation that you want to test. \n"
             "``domain`` : Test the impact of a mutation on the folding of a domain.\n"
             "``interface`` : Test the impact of a mutation on the interaction "
             "between two domains. \n"
             "``both`` : Test both a domain mutation and an interface mutation (DEFAULT) ")
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

    # Delete database
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
        help="Whether or not to leave the ``uniprot_sequence`` table untouched "
             "when clearing the schema.\n"
             "Only applicable if ``--drop_schema`` is set to ``True``.")
    parser_delete.set_defaults(func=delete_database)

    return parser


def elaspic_database():
    parser = get_database_parser()
    args = parser.parse_args()
    if 'func' not in args.__dict__:
        args = parser.parse_args(['--help'])
    conf.read_configuration_file(args.config_file)
    print("Running function '{}'...".format(args.func.__name__))
    args.func(args)


# %%
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] != 'database':
        elaspic_database()
    else:
        elaspic()
