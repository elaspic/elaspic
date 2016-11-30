"""ELASPIC DATABASE
"""
import os
import logging

from kmtools import system_tools

import elaspic

logger = logging.getLogger(__name__)


def elaspic_database(args):
    if args.config_file:
        elaspic.conf.read_configuration_file(args.config_file)
    elif args.connection_string:
        elaspic.conf.read_configuration_file(
            DATABASE={'connection_string': args.connection_string})
    else:
        raise Exception("Either 'config_file' or 'connection_string' must be specified!")
    print("Running function '{}'...".format(args.func.__name__))


def create_database(args):
    global elaspic
    if args.config_file:
        elaspic.conf.read_configuration_file(args.config_file)
    elif args.connection_string:
        elaspic.conf.read_configuration_file(
            DATABASE={'connection_string': args.connection_string})
    else:
        raise Exception("Either 'config_file' or 'connection_string' must be specified!")

    import elaspic.database
    db = elaspic.database.Database()
    db.create_database_tables(args.drop_schema)
    logger.info('Done!')


def load_data_to_database(args):
    global elaspic
    if args.config_file:
        elaspic.conf.read_configuration_file(args.config_file)
    elif args.connection_string:
        elaspic.conf.read_configuration_file(
            DATABASE={'connection_string': args.connection_string})
    else:
        raise Exception("Either 'config_file' or 'connection_string' must be specified!")

    import elaspic.database
    db = elaspic.database.Database()
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
            with system_tools.decompress(
                    os.path.join(args.data_folder, '{}.tsv.gz'.format(table.name))):
                db.copy_table_to_db(table.name, args.data_folder.rstrip('/'))
            print("Successfully loaded data from file '{}' to table '{}'"
                  .format('{}.tsv.gz'.format(table.name), table.name))


def test_database(args):
    raise NotImplementedError


def delete_database(args):
    global elaspic
    if args.config_file:
        elaspic.conf.read_configuration_file(args.config_file)
    elif args.connection_string:
        elaspic.conf.read_configuration_file(
            DATABASE={'connection_string': args.connection_string})
    else:
        raise Exception("Either 'config_file' or 'connection_string' must be specified!")

    import elaspic.database
    db = elaspic.database.Database()
    db.delete_database_tables(args.drop_schema, args.drop_uniprot_sequence)
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
        '-c', '--config_file', nargs='?', type=str,
        help='ELASPIC configuration file.')
    parser.add_argument(
        '--connection_string', nargs='?', type=str,
        help=('SQLAlchemy formatted string describing the connection to the database.'))
    subparsers = parser.add_subparsers(
        title='tasks',
        help='Maintenance tasks to perform')
    # parser.set_defaults(func=elaspic_database)

    # Create an empty database schema
    parser_create = subparsers.add_parser(
        name='create',
        description='Create an empty database')
    parser_create.add_argument(
        '--drop_schema', action='store_true', default=False,
        help=('Whether or not to first drop all existing tables from the database schema. \n'
              'WARNING: Choosing `True` will remove all existing data from the schema specified '
              'in your configuration file!!!'))
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
        '--drop_schema', action='store_true', default=False,
        help=('Whether or not to first drop all existing tables from the database schema. \n'
              'WARNING: Choosing `True` will remove all existing data from the schema specified '
              'in your configuration file!!!'))
    parser_delete.add_argument(
        '--drop_uniprot_sequence', action='store_true', default=False,
        help="Whether or not to leave the 'uniprot_sequence' table untouched when clearing "
             "the schema. Only applicable if '--clear_schema' is set to 'True'.")
    parser_delete.set_defaults(func=delete_database)
