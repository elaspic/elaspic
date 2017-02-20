"""ELASPIC RUN
"""
import os
import os.path as op
import logging
import argparse

import elaspic


logger = logging.getLogger(__name__)


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


def elaspic_run(args):
    validate_args(args)

    # Configure
    if args.config_file is not None:
        _configure_from_file(args)
    elif args.structure_file:
        _configure_standalone(args)
    elif args.uniprot_id:
        _configure_database(args)

    # Run
    if args.structure_file:
        _run_standalone_pipeline(args)
    elif args.uniprot_id:
        _run_database_pipeline(args)


def _configure_from_file(args):
    elaspic.conf.read_configuration_file(args.config_file)


def _configure_standalone(args):
    unique_temp_dir = op.abspath(op.join(os.getcwd(), '.elaspic'))
    os.makedirs(unique_temp_dir, exist_ok=True)
    elaspic.conf.read_configuration_file(
        DEFAULT={
            'unique_temp_dir': unique_temp_dir
        },
        EXTERNAL_DIRS={
            'pdb_dir': args.pdb_dir,
            'blast_db_dir': args.blast_db_dir,
            'archive_dir': args.archive_dir
        })


def _configure_database(args):
    logger.debug("_configure_database(%s)", args)
    elaspic.conf.read_configuration_file(
        DATABASE={
            'db_connection_string': args.connection_string
        },
        EXTERNAL_DIRS={
            'pdb_dir': args.pdb_dir,
            'blast_db_dir': args.blast_db_dir,
            'archive_dir': args.archive_dir,
        })


def _run_standalone_pipeline(args):
    # Run local pipeline
    pipeline = elaspic.StandalonePipeline(
        args.structure_file, args.sequence_file, args.mutations,
        mutation_format=args.mutation_format,
        run_type=args.run_type,
    )
    pipeline.run()


def _run_database_pipeline(args):
    # Run database pipeline
    if args.uniprot_domain_pair_ids:
        logger.debug('uniprot_domain_pair_ids: %s', args.uniprot_domain_pair_ids)
        uniprot_domain_pair_ids_asint = (
            [int(x) for x in args.uniprot_domain_pair_ids.split(',') if x]
        )
    else:
        uniprot_domain_pair_ids_asint = []

    # Run database pipeline
    import elaspic.database
    pipeline = elaspic.database.DatabasePipeline(
        args.uniprot_id, args.mutations,
        run_type=args.run_type,
        uniprot_domain_pair_ids=uniprot_domain_pair_ids_asint
    )
    pipeline.run()


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
        '--config_file', nargs='?', type=str,
        help='ELASPIC configuration file.')
    parser.add_argument(
        '-c', '--connection_string', nargs='?', type=str,
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
        '-v', '--verbose', action='count',
        help=('Specify verbosity level.'))

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
        '-t', '--run_type', nargs='?', type=str, default='all',
        choices=sorted(elaspic.pipeline._Pipeline._valid_run_types),
        help='Type of analysis to perform.')

    parser.set_defaults(func=elaspic_run)
