import os
import os.path as op
import json
import argparse
import pickle
import logging
import logging.config

from elaspic import conf, helper

logger = logging.getLogger(__name__)


# #################################################################################################
# ELASPIC RUN
def validate_args(args):
    if args.config_file and not os.path.isfile(args.config_file):
        raise Exception('The configuration file {} does not exist!'.format(args.config_file))

    if ((args.uniprot_id is None and args.structure_file is None) or
            (args.uniprot_id is not None and args.structure_file is not None)):
        raise Exception("""\
One of '-u' ('--uniprot_id') or '-p' ('--structure_file') must be specified!""")

    if (args.uniprot_id and (
            (args.config_file is None) and
            (args.pdb_dir is None or args.blast_db_dir is None or
             args.archive_dir is None))):
        raise Exception("""\
When using the database pipeline, \
you must either provide a configuration file ('-c', '--config_file') or \
'--pdb_dir', '--blast_db_dir', and '--archive_dir'.""")

    if args.sequence_file and not args.structure_file:
        raise Exception("""\
A template PDB file must be specified using the '--structure_file' option, \
when you specify a target sequence using the '--sequence_file' option!""")


def elaspic(args):
    validate_args(args)

    if args.uniprot_id:
        if args.config_file is not None:
            conf.read_configuration_file(args.config_file)
        else:
            conf.read_configuration_file(
                DATABASE={
                    'connection_string': args.connection_string
                },
                EXTERNAL_DIRS={
                    'pdb_dir': args.pdb_dir,
                    'blast_db_dir': args.blast_db_dir,
                    'archive_dir': args.archive_dir,
                })
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
    elif args.structure_file:
        unique_temp_dir = op.abspath(op.join(os.getcwd(), '.elaspic'))
        os.makedirs(unique_temp_dir, exist_ok=True)
        conf.read_configuration_file(
            DEFAULT={
                'unique_temp_dir': unique_temp_dir
            },
            EXTERNAL_DIRS={
                'pdb_dir': args.pdb_dir,
                'blast_db_dir': args.blast_db_dir,
                'archive_dir': args.archive_dir,
            })
        # Run local pipeline
        from elaspic import standalone_pipeline
        pipeline = standalone_pipeline.StandalonePipeline(
            args.structure_file, args.sequence_file, args.mutations,
            mutation_format=args.mutation_format,
        )
        if args.run_type in ['1', 'sequence']:
            pipeline.run_all_sequences()
        elif args.run_type in ['2', 'model']:
            pipeline.run_all_models()
        elif args.run_type in ['3', 'mutations']:
            pipeline.run_all_mutations()
        elif args.run_type in ['5']:
            pipeline.run()
        else:
            raise RuntimeError("Incorrect value for run_type: '{}'".format(args.run_type))


def configure_run_parser(sub_parsers):
    help = "Run ELASPIC"
    description = help + ""
    example = r"""

Examples
--------
$ elaspic run -p 4DKL.pdb -m A_M6A -n 1

$ elaspic run -u P00044 -m M1A -c config_file.ini

$ elaspic run -u P00044 -m M1A \
    --connection_string=mysql://user:pass@localhost/elaspic \
    --pdb_dir=/home/pdb/data/data/structures/divided/pdb \
    --blast_db_dir=/home/ncbi/blast/db \
    --archive_dir=/home/elaspic
"""
    parser = sub_parsers.add_parser(
        'run',
        help=help,
        description=description,
        epilog=example,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        '-c', '--config_file', nargs='?', type=str,
        help='ELASPIC configuration file.')
    parser.add_argument(
        '--connection_string', nargs='?', type=str,
        help=('SQLAlchemy formatted string describing the connection to the database.'))
    parser.add_argument(
        '--pdb_dir', nargs='?', type=str,
        help=("Folder containing PDB files in split format (e.g. 'ab/pdb1ab2.ent.gz')."))
    parser.add_argument(
        '--blast_db_dir', nargs='?', type=str,
        help=("Folder containing NCBI `nr` and `pdbaa` databases."))
    parser.add_argument(
        '--archive_dir', nargs='?', type=str,
        help=('Folder containing precalculated ELASPIC data.'))

    parser.add_argument(
        '-u', '--uniprot_id',
        help="The Uniprot ID of the protein that you want to mutate (e.g. 'P28223')."
             "This option relies on a local elaspic database, which has to be specified "
             "in the configuration file.")
    parser.add_argument(
        '-p', '--structure_file',
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
             "If you used '--structure_file', mutations must be provided using the chain "
             "and residue id (e.g. 'A_M1C,B_C20P' to mutate a residue with id '1' on chain A "
             "to Cysteine, and residue with id '20' on chain B to Proline).\n"
             "If you used '--sequence_file', mutations must be provided using the chain "
             "and residue INDEX (e.g. '1_M1C,2_C20P' to mutate the first residue in sequence 1 "
             "to Cysteine, and the 20th residue in sequence 2 to Proline).")
    parser.add_argument(
        '-n', '--mutation_format', nargs='?', default=None,
        help="Mutation format:\n"
             "  1. {pdb_chain}_{pdb_mutation},...\n"
             "  2. {pdb_chain}_{sequence_mutation},...\n"
             "  3. {sequence_pos}_{sequence_mutation}... (default)\n\n"
             "If `sequence_file` is None, this does not matter "
             "(always {pdb_chain}_{pdb_mutation})."
    )
    parser.add_argument(
        '-i', '--uniprot_domain_pair_ids', nargs='?', default='',
        help="List of uniprot_domain_pair_ids to analyse "
             "(useful if you want to restrict your analysis to only a handful of domains).")
    parser.add_argument(
        '-t', '--run_type', nargs='?', type=str, default='5',
        choices=['1', '2', '3', '4', '5', 'sequence', 'model', 'mutations'],
        help=('Type of analysis to perform: \n'
              '  1 / sequence: Calculate Provean only \n'
              '  2 / model: Create homololgy models only \n'
              '  3 / mutations: Evaluate mutations only \n'
              '  4: Create homology models and evaluate mutations \n'
              '  5: Calculate Provean, create homology models, and evaluate mutations \n'))

    parser.set_defaults(func=elaspic)


# #################################################################################################
# ELASPIC DATABASE

def elaspic_database(args):
    conf.read_configuration_file(args.config_file)
    print("Running function '{}'...".format(args.func.__name__))


def create_database(args):
    conf.read_configuration_file(args.config_file)
    from elaspic import elaspic_database
    db = elaspic_database.MyDatabase()
    db.create_database_tables(args.clear_schema, args.keep_uniprot_sequence)
    logger.info('Done!')


def load_data_to_database(args):
    conf.read_configuration_file(args.config_file)
    from elaspic import elaspic_database
    db = elaspic_database.MyDatabase()
    args.data_folder = args.data_folder.rstrip('/')
    table_names = args.data_files.split(',') if args.data_files else None
    dirpath, dirnames, filenames = next(os.walk(args.data_folder))
    for table in elaspic_database.Base.metadata.sorted_tables:
        if table_names is not None and table.name not in table_names:
            print("Skipping table '{}' because it was not included in the 'table_names' list..."
                  .format(table.name))
            continue
        if '{}.tsv'.format(table.name) in filenames:
            db.copy_table_to_db(table.name, args.data_folder)
            print("Successfully loaded data from file '{}' to table '{}'"
                  .format('{}.tsv'.format(table.name), table.name))
        elif '{}.tsv.gz'.format(table.name) in filenames:
            with helper.decompress(os.path.join(args.data_folder, '{}.tsv.gz'.format(table.name))):
                db.copy_table_to_db(table.name, args.data_folder.rstrip('/'))
            print("Successfully loaded data from file '{}' to table '{}'"
                  .format('{}.tsv.gz'.format(table.name), table.name))


def test_database(args):
    raise NotImplementedError


def delete_database(args):
    conf.read_configuration_file(args.config_file)
    from elaspic import elaspic_database
    db = elaspic_database.MyDatabase()
    db.delete_database_tables(args.drop_schema, args.keep_uniprot_sequence)
    logger.info('Done!')


def configure_database_parser(sub_parsers):
    help = "Perform database maintenance tasks"
    description = help + """
"""
    example = """
Examples:

    elaspic database -c config_file.ini create

"""
    parser = sub_parsers.add_parser(
        'database',
        help=help,
        description=description,
        epilog=example,
    )

    parser.add_argument(
        '-c', '--config_file', required=True,
        help='ELASPIC configuration file')
    subparsers = parser.add_subparsers(
        title='tasks',
        help='Maintenance tasks to perform')
    # parser.set_defaults(func=elaspic_database)

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


def elaspic_train(args):
    import pandas as pd
    from elaspic import DATA_DIR, machine_learning

    # Load data
    with open(op.join(DATA_DIR, 'core_options_p0.json')) as ifh:
        core_options_p0 = json.load(ifh)
    with open(op.join(DATA_DIR, 'core_options_p1.json')) as ifh:
        core_options_p1 = json.load(ifh)
    with open(op.join(DATA_DIR, 'interface_options_p0.json')) as ifh:
        interface_options_p0 = json.load(ifh)
    with open(op.join(DATA_DIR, 'interface_options_p1.json')) as ifh:
        interface_options_p1 = json.load(ifh)

    core_training_set = (
        pd.read_csv(
            op.join(DATA_DIR, 'core_training_set.tsv.gz'), sep='\t')
    )
    interface_training_set = (
        pd.read_csv(
            op.join(DATA_DIR, 'interface_training_set.tsv.gz'), sep='\t')
    )

    # Train predictors
    core_clf_p0 = (
        machine_learning.get_final_predictor(
            core_training_set, core_options_p0['features'], core_options_p0)
    )
    core_clf_p1 = (
        machine_learning.get_final_predictor(
            core_training_set, core_options_p1['features'], core_options_p1)
    )
    interface_clf_p0 = (
        machine_learning.get_final_predictor(
            interface_training_set, interface_options_p0['features'], interface_options_p0)
    )
    interface_clf_p1 = (
        machine_learning.get_final_predictor(
            interface_training_set, interface_options_p1['features'], interface_options_p1)
    )

    # Save predictors and features
    with open(op.join(DATA_DIR, 'ml_features_core_p0.json'), 'w') as ofh:
        json.dump(core_options_p0['features'], ofh)
    with open(op.join(DATA_DIR, 'ml_features_core_p1.json'), 'w') as ofh:
        json.dump(core_options_p1['features'], ofh)
    with open(op.join(DATA_DIR, 'ml_features_interface_p0.json'), 'w') as ofh:
        json.dump(interface_options_p0['features'], ofh)
    with open(op.join(DATA_DIR, 'ml_features_interface_p1.json'), 'w') as ofh:
        json.dump(interface_options_p1['features'], ofh)

    with open(op.join(DATA_DIR, 'ml_clf_core_p0.pickle'), 'wb') as ofh:
        pickle.dump(core_clf_p0, ofh)
    with open(op.join(DATA_DIR, 'ml_clf_core_p1.pickle'), 'wb') as ofh:
        pickle.dump(core_clf_p1, ofh)
    with open(op.join(DATA_DIR, 'ml_clf_interface_p0.pickle'), 'wb') as ofh:
        pickle.dump(interface_clf_p0, ofh)
    with open(op.join(DATA_DIR, 'ml_clf_interface_p1.pickle'), 'wb') as ofh:
        pickle.dump(interface_clf_p1, ofh)


def configure_train_parser(sub_parsers):
    help = "Train the ELASPIC classifiers"
    description = help + """
"""
    example = """
Examples:

    elaspic train

"""
    parser = sub_parsers.add_parser(
        'train',
        help=help,
        description=description,
        epilog=example,
    )
    parser.set_defaults(func=elaspic_train)


def main():
    parser = argparse.ArgumentParser(
        prog='elaspic',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    sub_parsers = parser.add_subparsers(
        title='command',
        help=''
    )
    configure_run_parser(sub_parsers)
    configure_database_parser(sub_parsers)
    configure_train_parser(sub_parsers)
    args = parser.parse_args()
    if 'func' not in args.__dict__:
        args = parser.parse_args(['--help'])
    args.func(args)


if __name__ == '__main__':
    import sys
    sys.exit(main())
