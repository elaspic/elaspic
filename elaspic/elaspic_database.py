import os
import os.path as op
import subprocess
import datetime
import six
import logging
import shutil
from contextlib import contextmanager
import pandas as pd
import sqlalchemy as sa
from kmtools.db_tools import parse_connection_string, make_connection_string
from . import helper, errors, conf
from .elaspic_database_tables import (
    Base,
    UniprotDomain, UniprotDomainModel, UniprotDomainMutation,
    UniprotDomainPair, UniprotDomainPairModel,
    UniprotDomainPairMutation,
)

logger = logging.getLogger(__name__)


def enable_sqlite_foreign_key_checks(engine):
    from sqlalchemy import event

    # Enable foreign key contraints
    def _fk_pragma_on_connect(dbapi_con, con_record):
        dbapi_con.execute('pragma foreign_keys=ON')
    event.listen(engine, 'connect', _fk_pragma_on_connect)

    # Enable the write-ahead lock so that reads can occur simultaneously with writes
    def _fk_pragma_on_connect(dbapi_con, con_record):
        dbapi_con.execute('PRAGMA journal_mode=WAL')
    event.listen(engine, 'connect', _fk_pragma_on_connect)

    # Set a longer timeout duration
    def _fk_pragma_on_connect(dbapi_con, con_record):
        dbapi_con.execute('pragma busy_timeout=60000')  # 60 sec
    event.listen(engine, 'connect', _fk_pragma_on_connect)


# Get the session that will be used for all future queries.
# `expire_on_commit` so that you keep all the table objects even after the session closes.
Session = sa.orm.sessionmaker(expire_on_commit=False)
# Session = sa.orm.scoped_session(sa.orm.sessionmaker(expire_on_commit=False))


class MyDatabase(object):
    """
    """

    def __init__(self, echo=False):
        self.engine = self.get_engine(echo=echo)
        self.configure_session()

        logger.info(
            "Using precalculated data from the following folder: '{archive_dir}'"
            .format(**conf.CONFIGS)
        )

    def get_engine(self, echo=False):
        """Get an SQLAlchemy engine that can be used to connect to the database."""
        sa_opts = {
            'echo': echo,
            'pool_size': 1,
        }
        if conf.CONFIGS['db_type'] == 'sqlite':
            sa_opts['isolation_level'] = 'READ UNCOMMITTED'
        elif conf.CONFIGS['db_type'] == 'mysql':
            sa_opts['isolation_level'] = 'READ UNCOMMITTED'
            sa_opts['pool_recycle'] = 3600
        elif conf.CONFIGS['db_type'] == 'postgresql':
            sa_opts['pool_recycle'] = 3600
        else:
            raise Exception("Unsupported 'db_type': '{}'!".format(conf.CONFIGS['db_type']))
        engine = sa.create_engine(conf.CONFIGS['connection_string'], **sa_opts)
        if conf.CONFIGS['db_type'] == 'sqlite':
            enable_sqlite_foreign_key_checks(engine)
        logger.info("Opened database connection using engine: '{}'".format(engine))
        return engine

    def configure_session(self):
        """
        Configure the Session class to use the current engine.

        `autocommit` and `autoflush` are enabled for the `sqlite` database in order to improve
        performance.

        """
        global Session
        if conf.CONFIGS['db_type'] == 'sqlite':
            autocommit = False  # True
            autoflush = True  # True
            # retry_on_failure = True
        elif conf.CONFIGS['db_type'] in ['postgresql', 'mysql']:
            autocommit = False
            autoflush = True
            # retry_on_failure = True
        Session.configure(bind=self.engine, autocommit=autocommit, autoflush=autoflush)

    @contextmanager
    def session_scope(self):
        """Provide a transactional scope around a series of operations.

        Enables the following construct: ``with self.session_scope() as session:``.
        """
        session = Session()
        try:
            yield session
            if not conf.CONFIGS['db_is_immutable']:
                session.commit()
        except:
            if not conf.CONFIGS['db_is_immutable']:
                session.rollback()
            raise
        finally:
            session.expunge_all()
            session.close()

    def create_database_schema(self, db_schema):
        """Create ELASPIC database schema."""
        # Create engine without a default schema
        engine = sa.create_engine(
            make_connection_string(**{
                **parse_connection_string(conf.CONFIGS['connection_string']),
                'db_schema': '',
            }))
        sql_command = "CREATE SCHEMA IF NOT EXISTS `{}`;".format(db_schema)
        logger.debug("sql_command: '{}'".format(sql_command))
        engine.execute(sql_command)

    def drop_database_schema(self, db_schema):
        """Drop ELASPIC database schema."""
        # Create engine without a default schema
        engine = sa.create_engine(
            make_connection_string(**{
                **parse_connection_string(conf.CONFIGS['connection_string']),
                'db_schema': '',
            }))
        sql_command = "DROP SCHEMA IF EXISTS {};".format(db_schema)
        logger.debug("sql_command: '{}'".format(sql_command))
        engine.execute(sql_command)

    def create_database_tables(self, drop_schema=False):
        """Create a new database in the schema specified by the ``schema_version`` global variable.

        If ``clear_schema == True``, remove all the tables in the schema first.

        .. warning::

            Using this function with an existing database can lead to loss of data.
            Make sure that you know what you are doing!

        Parameters
        ----------
        clear_schema : bool
            Whether or not to delete all tables in the database schema before creating new tables.
        keep_uniprot_sequence : bool
            Whether or not to keep the `uniprot_sequence` table.
            Only relevant if `clear_schema` is `True`.
        """
        if drop_schema:
            self.drop_database_schema(conf.CONFIGS['db_schema'])

        self.create_database_schema(conf.CONFIGS['db_schema'])

        # Create all tables, creating schema as neccessary
        for table in Base.metadata.sorted_tables:
            table.create(self.engine, checkfirst=True)
        logger.debug('Database tables were created successfully.')

    def delete_database_tables(self, drop_schema=False, drop_uniprot_sequence=False):
        """.

        Parameters
        ----------
        drop_schema : bool
            Whether or not to drop the schema after dropping the tables.
        keep_uniprot_sequence : bool
            Wheter or not to keep the table (and schema) containing uniprot sequences.
        """
        if drop_schema:
            if conf.CONFIGS['db_type'] == 'sqlite':
                os.remove(conf.CONFIGS['db_schema'])
            else:
                self.engine.execute('drop schema {db_schema};'.format(**conf.CONFIGS))
            logger.info(
                "Successfully removed database schema: {db_schema}"
                .format(**conf.CONFIGS))
            return

        # Remove tables one by one
        for table in reversed(Base.metadata.sorted_tables):
            if table.name != 'uniprot_sequence' or drop_uniprot_sequence:
                conf.CONFIGS['table_name'] = table.name
                self.engine.execute(
                    'drop table if exists {db_schema}.{table_name};'.format(**conf.CONFIGS))
                self.engine.execute(
                    'drop table if exists {db_schema}.{table_name};'.format(**conf.CONFIGS))

    # %%
    mysql_load_table_template = (
        r"""mysql --local-infile --host={db_url} --user={db_username} --password={db_password} """
        r"""{table_db_schema} -e "{sql_command}" """
    )

    psql_load_table_template = (
        r"""PGPASSWORD={db_password} psql -h {db_url} -p {db_port} -U {db_username} """
        r"""-d {db_database} -c "{sql_command}" """
    )

    # Need to double up on '\\'
    mysql_command_template = (
        r"""load data local infile '{table_folder}/{table_name}.tsv' """
        r"""into table {table_db_schema}.{table_name} """
        r"""fields terminated by '\t' escaped by '\\\\' lines terminated by '\n'; """
    )

    psql_command_template = (
        r"""\\copy {table_db_schema}.{table_name} """
        r"""from '{table_folder}/{table_name}.tsv' """
        r"""with csv delimiter E'\t' null '\N' escape '\\'; """
    )

    sqlite_table_filename = '{table_folder}/{table_name}.tsv'

    def _load_data_into_sqlite(self, configs):
        table_df = pd.read_csv(
            self.sqlite_table_filename.format(**configs),
            sep='\t', na_values='\\N',
            # escapechar='\\', # escapes the `na_values` character and causes problems
            names=Base.metadata.tables[configs['table_name']].columns.keys())
        table_df.to_sql(configs['table_name'], self.engine, index=False, if_exists='append')

    def _run_create_table_system_command(self, system_command):
        if conf.CONFIGS['debug']:
            logger.debug(system_command)
        result, error_message, return_code = helper.subprocess_check_output(system_command)
        if return_code != 0:
            logger.error(result)
            raise Exception(error_message)

    def copy_table_to_db(self, table_name, table_folder):
        """Copy data from a ``.tsv`` file to a table in the database."""
        cmd_options = conf.CONFIGS.copy()
        cmd_options['table_name'] = table_name
        cmd_options['table_folder'] = table_folder

        def _format_configs():
            cmd_options['table_db_schema'] = cmd_options['db_schema']

        logger.info("Copying '{table_name}' to '{db_type}' database...".format(**cmd_options))
        if cmd_options['db_type'] == 'sqlite':
            self._load_data_into_sqlite(cmd_options)
        elif cmd_options['db_type'] == 'mysql':
            _format_configs()
            cmd_options['sql_command'] = self.mysql_command_template.format(**cmd_options)
            system_command = self.mysql_load_table_template.format(**cmd_options)
            self._run_create_table_system_command(system_command)
        elif cmd_options['db_type'] == 'postgresql':
            _format_configs()
            cmd_options['sql_command'] = self.psql_command_template.format(**cmd_options)
            system_command = self.psql_load_table_template.format(**cmd_options)
            self._run_create_table_system_command(system_command)
        else:
            raise Exception("Unsupported database type: '{}'".format(cmd_options['db_type']))

    @helper.retry_database
    def get_rows_by_ids(self, row_object, row_object_identifiers, row_object_identifier_values):
        """Get the rows from the table `row_object` identified by keys `row_object_identifiers`."""
        with self.session_scope() as session:
            if len(row_object_identifiers) != len(row_object_identifier_values):
                raise Exception(
                    'The number of identifiers and the number of identifier '
                    'values must be the same.')
            if len(row_object_identifiers) > 3:
                raise Exception(
                    'Too many identifiers provied. The function is hard-coded '
                    'to handle at most three identifiers.')
            if len(row_object_identifiers) == 1:
                row_instances = (
                    session.query(row_object)
                    .filter(row_object_identifiers[0] == row_object_identifier_values[0])
                    .all())
            if len(row_object_identifiers) == 2:
                row_instances = (
                    session.query(row_object)
                    .filter(row_object_identifiers[0] == row_object_identifier_values[0])
                    .filter(row_object_identifiers[1] == row_object_identifier_values[1])
                    .all())
            if len(row_object_identifiers) == 3:
                row_instances = (
                    session.query(row_object)
                    .filter(row_object_identifiers[0] == row_object_identifier_values[0])
                    .filter(row_object_identifiers[1] == row_object_identifier_values[1])
                    .filter(row_object_identifiers[2] == row_object_identifier_values[2])
                    .all())
            return row_instances

    @helper.retry_database
    def get_uniprot_domain(self, uniprot_id, copy_data=False):
        """
        """
        with self.session_scope() as session:
            uniprot_domains = (
                session
                .query(UniprotDomain)
                .filter(UniprotDomain.uniprot_id == uniprot_id)
                # .options(sa.orm.joinedload('template').joinedload('model'))
                # .options(sa.orm.joinedload('template', innerjoin=True))
                .limit(100)
                .all()
            )

        archive_dir = conf.CONFIGS['archive_dir']
        archive_type = conf.CONFIGS['archive_type']
        d_idx = 0
        while d_idx < len(uniprot_domains):
            d = uniprot_domains[d_idx]
            if not d.template:
                logger.debug(
                    'Skipping uniprot domain with id {} because it does not '
                    'have a structural template...'
                    .format(d.uniprot_domain_id))
                del uniprot_domains[d_idx]
                continue
            # Copy precalculated Provean data
            if copy_data:
                try:
                    self._copy_provean(
                        d, archive_dir, archive_type)
                except subprocess.CalledProcessError as e:
                    logger.error(e)
                    logger.error('Failed to copy provean supporting set!')
                    d.uniprot_sequence.provean.provean_supset_filename = ''
            # Copy precalculated homology models
            if copy_data:
                try:
                    self._copy_uniprot_domain_data(
                        d, d.path_to_data, archive_dir, archive_type)
                except subprocess.CalledProcessError as e:
                    logger.error(e)
                    logger.error('Failed to copy the domain alignment and / or homology model!')
                    d.template.model.alignment_filename = None
                    d.template.model.model_filename = None
            d_idx += 1

        return uniprot_domains

    @helper.retry_database
    def get_uniprot_domain_pair(self, uniprot_id, copy_data=False, uniprot_domain_pair_ids=[]):
        """
        """
        with self.session_scope() as session:
            uniprot_domain_pairs_query = (
                session.query(UniprotDomainPair)
                .filter(sa.or_(
                    sa.text("uniprot_id_1='{}'".format(uniprot_id)),
                    sa.text("uniprot_id_2='{}'".format(uniprot_id))))
            )
            if uniprot_domain_pair_ids:
                uniprot_domain_pairs_query = (
                    uniprot_domain_pairs_query
                    .filter(UniprotDomainPair.uniprot_domain_pair_id.in_(uniprot_domain_pair_ids))
                )
            uniprot_domain_pairs = (
                uniprot_domain_pairs_query
                # .options(sa.orm.joinedload('template', innerjoin=True).joinedload('model'))
                # .options(sa.orm.joinedload('template', innerjoin=True))
                .limit(100)
                .all()
            )

        # The above SQL query may result in duplicates if we have homodimers.
        # So we need to remove possible dimers.
        _seen = set()
        uniprot_domain_pairs = [
            d for d in uniprot_domain_pairs
            if d.uniprot_domain_pair_id not in _seen and
            not _seen.add(d.uniprot_domain_pair_id)
        ]
        #
        archive_dir = conf.CONFIGS['archive_dir']
        archive_type = conf.CONFIGS['archive_type']
        d_idx = 0
        while d_idx < len(uniprot_domain_pairs):
            d = uniprot_domain_pairs[d_idx]
            if not d.template:
                logger.debug(
                    'Skipping uniprot domain pair with id {} because it does not '
                    'have a structural template...'
                    .format(d.uniprot_domain_pair_id))
                del uniprot_domain_pairs[d_idx]
                continue
            # Copy precalculated Provean data
            if copy_data:
                if d.uniprot_id_1 == uniprot_id:
                    ud = d.uniprot_domain_1
                elif d.uniprot_id_2 == uniprot_id:
                    ud = d.uniprot_domain_2
                try:
                    self._copy_provean(ud, archive_dir, archive_type)
                except subprocess.CalledProcessError as e:
                    logger.error(e)
                    logger.error('Failed to copy provean supporting set!')
                    d.uniprot_sequence.provean.provean_supset_filename = ''
            # Copy precalculated homology models
            if copy_data:
                try:
                    self._copy_uniprot_domain_pair_data(
                        d, d.path_to_data, archive_dir, archive_type)
                except subprocess.CalledProcessError as e:
                    logger.error(e)
                    logger.error('Failed to copy domain pair alignments and / or homology model!')
                    d.template.model.alignment_filename_1 = None
                    d.template.model.alignment_filename_2 = None
                    d.template.model.model_filename = None
            d_idx += 1

        return uniprot_domain_pairs

    def _copy_uniprot_domain_data(self, d, path_to_data, archive_dir, archive_type):
        if path_to_data is None:
            logger.error('Cannot copy uniprot domain data because `path_to_data` is None')
            return
        if (d.template and
                d.template.model and
                d.template.model.alignment_filename and
                d.template.model.model_filename):
            try:
                tmp_save_path = op.join(conf.CONFIGS['archive_temp_dir'], path_to_data)
                archive_save_path = op.join(archive_dir, path_to_data)
                args = [tmp_save_path] + d.template.model.alignment_filename.split('/')[:-1]
                path_to_alignment = op.join(*args)
                subprocess.check_call(
                    "umask ugo=rwx; mkdir -m 777 -p '{}'"
                    .format(path_to_alignment), shell=True)
                shutil.copyfile(
                    op.join(archive_save_path, d.template.model.alignment_filename),
                    op.join(tmp_save_path, d.template.model.alignment_filename),
                )
                shutil.copyfile(
                    op.join(archive_save_path, d.template.model.model_filename),
                    op.join(tmp_save_path, d.template.model.model_filename),
                )
            except FileNotFoundError:
                if archive_type != '7zip':
                    raise
                # Extract files from a 7zip archive
                filenames = [
                    op.join(path_to_data, d.template.model.alignment_filename),
                    op.join(path_to_data, d.template.model.model_filename),
                ]
                path_to_7z = os.path.join(
                    archive_dir,
                    'uniprot_domain',
                    'uniprot_domain.7z'
                )
                self._extract_files_from_7zip(path_to_7z, filenames)

    def _copy_uniprot_domain_pair_data(self, d, path_to_data, archive_dir, archive_type):
        if path_to_data is None:
            logger.error('Cannot copy uniprot domain pair data because `path_to_data` is None')
            return
        if (d.template and
                d.template.model and
                d.template.model.alignment_filename_1 and
                d.template.model.alignment_filename_2 and
                d.template.model.model_filename):
            try:
                tmp_save_path = op.join(conf.CONFIGS['archive_temp_dir'], path_to_data)
                archive_save_path = op.join(archive_dir, path_to_data)
                path_to_alignment_1 = op.join(
                    tmp_save_path,
                    *d.template.model.alignment_filename_1.split('/')[:-1])
                path_to_alignment_2 = op.join(
                    tmp_save_path,
                    *d.template.model.alignment_filename_2.split('/')[:-1])
                os.makedirs(path_to_alignment_1, exist_ok=True)
                os.makedirs(path_to_alignment_2, exist_ok=True)
                shutil.copyfile(
                    op.join(archive_save_path, d.template.model.alignment_filename_1),
                    op.join(tmp_save_path, d.template.model.alignment_filename_1))
                shutil.copyfile(
                    op.join(archive_save_path, d.template.model.alignment_filename_2),
                    op.join(tmp_save_path, d.template.model.alignment_filename_2))
                shutil.copyfile(
                    op.join(archive_save_path, d.template.model.model_filename),
                    op.join(tmp_save_path, d.template.model.model_filename))
            except FileNotFoundError:
                if archive_type != '7zip':
                    raise
                # Extract files from a 7zip archive
                filenames = [
                    op.join(path_to_data, d.template.model.alignment_filename_1),
                    op.join(path_to_data, d.template.model.alignment_filename_2),
                    op.join(path_to_data, d.template.model.model_filename),
                ]
                path_to_7zip = os.path.join(
                    archive_dir,
                    'uniprot_domain_pair',
                    'uniprot_domain_pair.7z'
                )
                self._extract_files_from_7zip(path_to_7zip, filenames)

    def _copy_provean(self, ud, archive_dir, archive_type):
        if not (ud.uniprot_sequence and
                ud.uniprot_sequence.provean and
                ud.uniprot_sequence.provean.provean_supset_filename):
            logger.warning('Provean supset is missing for domain: {}'.format(ud.uniprot_domain_id))
            return

        try:
            # archive_type != '7zip' or extraction failed
            os.makedirs(
                op.dirname(op.join(
                    conf.CONFIGS['archive_temp_dir'],
                    get_uniprot_base_path(ud),
                    ud.uniprot_sequence.provean.provean_supset_filename)),
                exist_ok=True
            )
            shutil.copyfile(
                op.join(
                    archive_dir,
                    get_uniprot_base_path(ud),
                    ud.uniprot_sequence.provean.provean_supset_filename),
                op.join(
                    conf.CONFIGS['archive_temp_dir'],
                    get_uniprot_base_path(ud),
                    ud.uniprot_sequence.provean.provean_supset_filename))
            shutil.copyfile(
                op.join(
                    archive_dir,
                    get_uniprot_base_path(ud),
                    ud.uniprot_sequence.provean.provean_supset_filename + '.fasta'),
                op.join(
                    conf.CONFIGS['archive_temp_dir'],
                    get_uniprot_base_path(ud),
                    ud.uniprot_sequence.provean.provean_supset_filename + '.fasta'))
        except FileNotFoundError:
            if archive_type != '7zip':
                raise
            filenames = [
                op.join(
                    get_uniprot_base_path(ud),
                    ud.uniprot_sequence.provean.provean_supset_filename),
                op.join(
                    get_uniprot_base_path(ud),
                    ud.uniprot_sequence.provean.provean_supset_filename + '.fasta'),
            ]
            path_to_7zip = os.path.join(
                archive_dir,
                'provean',
                'provean.7z'
            )
            logger.debug('path_to_7zip: {}'.format(path_to_7zip))
            try:
                self._extract_files_from_7zip(path_to_7zip, filenames)
                return
            except errors.Archive7zipFileNotFoundError as e:
                logger.debug(e)
                logger.debug(
                    "Could not extract these files: {} from this archive: {}\n"
                    "Checking if they are in the archive directory..."
                    .format(filenames, path_to_7zip)
                )

    @helper.retry_archive
    def _extract_files_from_7zip(self, path_to_7zip, filenames_in):
        """Extract files to `config['archive_temp_dir']`."""
        logger.debug("Extracting the following files: {}".format(filenames_in))
        filenames = [
            f for f in filenames_in
            if not op.isfile(op.join(conf.CONFIGS['archive_temp_dir'], f))
        ]
        if not filenames:
            logger.debug("All files already been extracted. Done!")
            return
        system_command = "7za x '{path_to_7zip}' '{files}' -y".format(
            path_to_7zip=path_to_7zip,
            files="' '".join(filenames)
        )
        logger.debug('System command:\n{}'.format(system_command))
        p = helper.run(system_command, cwd=conf.CONFIGS['archive_temp_dir'])
        if 'No files to process' in p.stdout:
            raise errors.Archive7zipFileNotFoundError(p.stdout, p.stderr, p.returncode)
        if p.returncode != 0:
            raise errors.Archive7zipError(p.stdout, p.stderr, p.returncode)

    @helper.retry_database
    def get_uniprot_mutation(self, d, mutation, uniprot_id=None, copy_data=False):
        """
        """
        if isinstance(d, UniprotDomain):
            with self.session_scope() as session:
                uniprot_mutation = (
                    session.query(UniprotDomainMutation)
                    .filter(
                        (UniprotDomainMutation.uniprot_domain_id == d.uniprot_domain_id) &
                        (UniprotDomainMutation.mutation == mutation))
                    .scalar())
        elif isinstance(d, UniprotDomainPair) and isinstance(uniprot_id, six.string_types):
            with self.session_scope() as session:
                uniprot_mutation = (
                    session.query(UniprotDomainPairMutation)
                    .filter(
                        (UniprotDomainPairMutation.uniprot_id == uniprot_id) &
                        (UniprotDomainPairMutation.uniprot_domain_pair_id ==
                         d.uniprot_domain_pair_id) &
                        (UniprotDomainPairMutation.mutation == mutation))
                    .scalar())
        else:
            logger.debug('d: {}\td type: {}'.format(d, type(d)))
            raise Exception('Not enough arguments, or the argument types are incorrect!')

        if copy_data:
            try:
                self._copy_mutation_data(
                    uniprot_mutation, d.path_to_data, conf.CONFIGS['archive_dir'])
            except (subprocess.CalledProcessError,
                    errors.Archive7zipError,
                    errors.Archive7zipFileNotFoundError) as e:
                logger.error(e)
                uniprot_mutation.model_filename_wt = None
        return uniprot_mutation

    def _copy_mutation_data(self, mutation, path_to_data, archive_dir):
        if mutation and mutation.model_filename_wt:
            if archive_dir.endswith('.7z'):
                # Extract files from a 7zip archive
                filenames = [
                    op.join(path_to_data, mutation.model_filename_wt),
                    op.join(path_to_data, mutation.model_filename_mut),
                ]
                self._extract_files_from_7zip(archive_dir, filenames)
            else:
                tmp_save_path = op.join(conf.CONFIGS['archive_temp_dir'], path_to_data)
                archive_save_path = op.join(archive_dir, path_to_data)
                path_to_mutation = op.dirname(op.join(tmp_save_path, mutation.model_filename_wt))
                subprocess.check_call(
                    "umask ugo=rwx; mkdir -m 777 -p '{}'".format(path_to_mutation), shell=True)
                shutil.copyfile(
                    op.join(archive_save_path, mutation.model_filename_wt),
                    op.join(tmp_save_path, mutation.model_filename_wt))
                shutil.copyfile(
                    op.join(archive_save_path, mutation.model_filename_mut),
                    op.join(tmp_save_path, mutation.model_filename_mut))

    @helper.retry_database
    def remove_model(self, d):
        """Remove a model from the database.

        Do this if you realized that the model you built is incorrect
        or that some of the data is missing.

        Raises
        ------
        errors.ModelHasMutationsError
            The model you are trying to delete has precalculated mutations,
            so it can't be that bad. Delete those mutations and try again.
        """
        with self.session_scope() as session:
            if isinstance(d, UniprotDomain):
                session.execute(
                    'delete from {0}.uniprot_domain_model where uniprot_domain_id = {1}'
                    .format(conf.CONFIGS['db_schema'], d.uniprot_domain_id))
            elif isinstance(d, UniprotDomainPair):
                session.execute(
                    'delete from {0}.uniprot_domain_pair_model where uniprot_domain_pair_id = {1}'
                    .format(conf.CONFIGS['db_schema'], d.uniprot_domain_pair_id))
            else:
                raise Exception("'d' is of incorrect type!")

    # %% Add objects to the database
    @helper.retry_database
    def merge_row(self, row_instance):
        """Add a list of rows (`row_instances`) to the database."""
        if not conf.CONFIGS['db_is_immutable']:
            with self.session_scope() as session:
                if not isinstance(row_instance, (tuple, list)):
                    session.merge(row_instance)
                else:
                    for instance in row_instance:
                        session.merge(instance)

    def merge_provean(self, provean, provean_supset_file, path_to_data):
        """Add provean score to the database."""
        assert op.isfile(provean_supset_file)
        assert op.isfile(provean_supset_file + '.fasta')

        archive_dir = conf.CONFIGS['archive_dir']
        logger.debug(
            'Moving provean supset to the output folder: {}'
            .format(op.join(archive_dir, path_to_data)))
        helper.makedirs(op.join(archive_dir, path_to_data), mode=0o777)
        helper.copyfile(
            provean_supset_file,
            op.join(archive_dir, path_to_data, provean.provean_supset_filename),
            mode=0o666)
        helper.copyfile(
            provean_supset_file + '.fasta',
            op.join(archive_dir, path_to_data, provean.provean_supset_filename + '.fasta'),
            mode=0o666)
        self.merge_row(provean)

    def merge_model(self, d, files_dict={}):
        """Add MODELLER models to the database."""
        # Save a copy of the alignment to the export folder
        if files_dict:
            archive_dir = conf.CONFIGS['archive_dir']
            archive_save_path = op.join(archive_dir, d.path_to_data)
            helper.makedirs(archive_save_path, mode=0o777)
            if d.template.model.model_filename is not None:
                # Save alignments
                if isinstance(d.template.model, UniprotDomainModel):
                    helper.copyfile(
                        files_dict['alignment_files'][0],
                        op.join(archive_save_path, d.template.model.alignment_filename),
                        mode=0o666
                    )
                elif isinstance(d.template.model, UniprotDomainPairModel):
                    helper.copyfile(
                        files_dict['alignment_files'][0],
                        op.join(archive_save_path, d.template.model.alignment_filename_1),
                        mode=0o666
                    )
                    helper.copyfile(
                        files_dict['alignment_files'][1],
                        op.join(archive_save_path, d.template.model.alignment_filename_2),
                        mode=0o666
                    )
                # Save the modelled structure
                helper.copyfile(
                    files_dict['model_file'],
                    op.join(archive_save_path, d.template.model.model_filename),
                    mode=0o666
                )
        self.merge_row(d.template.model)

    def merge_mutation(self, mut, path_to_data=False):
        """
        """
        mut.mut_date_modified = datetime.datetime.utcnow()
        if path_to_data and (mut.model_filename_wt is not None):
            archive_temp_save_dir = op.join(conf.CONFIGS['archive_temp_dir'], path_to_data)
            archive_save_dir = op.join(conf.CONFIGS['archive_dir'], path_to_data)
            # Save Foldx structures
            if mut.model_filename_wt and mut.model_filename_mut:
                helper.makedirs(
                    op.dirname(op.join(archive_save_dir, mut.model_filename_wt)),
                    mode=0o777)
                helper.copyfile(
                    op.join(archive_temp_save_dir, mut.model_filename_wt),
                    op.join(archive_save_dir, mut.model_filename_wt),
                    mode=0o666
                )
                helper.copyfile(
                    op.join(archive_temp_save_dir, mut.model_filename_mut),
                    op.join(archive_save_dir, mut.model_filename_mut),
                    mode=0o666
                )
        self.merge_row(mut)


def get_uniprot_base_path(d=None, uniprot_name=None, uniprot_id=None):
    """Return the name of the subfolder for storing protein information.

    Parameters
    ----------
    d : UniprotDomain | UniprotDomainPair | None
    uniprot_name : str
    uniprot_id : str
    """
    if isinstance(d, UniprotDomain):
        uniprot_name = d.uniprot_sequence.uniprot_name
        uniprot_id = d.uniprot_id
    elif isinstance(d, UniprotDomainPair):
        uniprot_name = d.uniprot_domain_1.uniprot_sequence.uniprot_name
        uniprot_id = d.uniprot_domain_1.uniprot_id
    else:
        assert uniprot_name is not None and uniprot_id is not None

    # TODO: Screw the splitting of uniprot ids!
    uniprot_base_path = (
        '{organism_name}/{uniprot_id_part_1}/{uniprot_id_part_2}/{uniprot_id_part_3}/'
        .format(
            organism_name=uniprot_name.split('_')[-1].lower(),
            uniprot_id_part_1=uniprot_id[:3],
            uniprot_id_part_2=uniprot_id[3:5],
            uniprot_id_part_3=uniprot_id,))
    return uniprot_base_path


def get_uniprot_domain_path(d=None, **vargs):
    """Return the name of the subfolder for storing protein *domain* information.

    Parameters
    ----------
    ### Domain
    d : UniprotDomain | None
    pfam_clan : str
    alignment_def : str

    #### Domain pair
    d : UniprotDomainPair | None
    pfam_clan_1 : str
    alignment_def_1 : str
    pfam_clan_2 : str
    alignment_def_2 : str
    uniprot_id_2 : str
    """
    if d is not None:
        if isinstance(d, UniprotDomain):
            vargs['pfam_clan'] = d.pfam_clan
            vargs['alignment_def'] = d.alignment_def
        elif isinstance(d, UniprotDomainPair):
            vargs['pfam_clan_1'] = d.uniprot_domain_1.pfam_clan
            vargs['alignment_def_1'] = d.uniprot_domain_1.alignment_def
            vargs['pfam_clan_2'] = d.uniprot_domain_2.pfam_clan
            vargs['alignment_def_2'] = d.uniprot_domain_2.alignment_def
            vargs['uniprot_id_2'] = d.uniprot_domain_2.uniprot_id
        else:
            raise Exception

    if 'pfam_clan' in vargs:
        vargs['alignment_def'] = vargs['alignment_def'].replace(':', '-')
        uniprot_domain_path = (
            '{pfam_clan:.36}.{alignment_def}/'
            .format(**vargs)
        )
    elif 'pfam_clan_1' in vargs:
        vargs['alignment_def_1'] = vargs['alignment_def_1'].replace(':', '-')
        vargs['alignment_def_2'] = vargs['alignment_def_2'].replace(':', '-')
        uniprot_domain_path = (
            '{pfam_clan_1:.36}.{alignment_def_1}/'
            '{pfam_clan_2:.36}.{alignment_def_2}/'
            '{uniprot_id_2}/'
            .format(**vargs)
        )
    else:
        raise Exception

    return uniprot_domain_path
