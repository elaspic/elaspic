# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import next
from builtins import object

import os
import re
import stat
import urllib.request
import urllib.error
import subprocess
import json
import datetime
import time
import pickle
import six
import functools
import inspect
from contextlib import contextmanager
from collections import deque

import pandas as pd
import sqlalchemy as sa
import sqlalchemy.exc as sa_exc
from retrying import retry

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#import parse_pfamscan
from . import conf
from . import helper_functions as hf
from . import errors as error
from .database_tables import (
    Base, 
    Domain, DomainContact, UniprotSequence, Provean,
    UniprotDomain, UniprotDomainTemplate, UniprotDomainModel, UniprotDomainMutation,
    UniprotDomainPair, UniprotDomainPairTemplate, UniprotDomainPairModel, UniprotDomainPairMutation,
)


#%%
def decorate_all_methods(decorator):
    """Decorate all methods of a class with `decorator`.
    """
    def apply_decorator(cls):
        for k, f in cls.__dict__.items():
            if inspect.isfunction(f):
                setattr(cls, k, decorator(f))
        return cls
    return apply_decorator


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
        dbapi_con.execute('pragma busy_timeout=60000') # 60 sec
    event.listen(engine, 'connect', _fk_pragma_on_connect)
    
    

#%%
#: Get the session that will be used for all future queries.
#: `expire_on_commit` so that you keep all the table objects even after the session closes.
Session = sa.orm.sessionmaker(expire_on_commit=False)
#Session = scoped_session(sa.orm.sessionmaker(expire_on_commit=False))
        
        
#%%
class MyDatabase(object):
    """
    """
    
    def __init__(self, configs=conf.configs, logger=hf.get_logger()):
        """
        Parameters
        ----------
        configs : dict
            ELASPIC configuration options specified in :py:data:`elaspic.conf.configs`
        """
        self.configs = configs.copy()
        self.logger = logger        
        self.engine = self.get_engine()
        self.Session = self.configure_session()
        
        self.logger.info(
            "Using precalculated data from the following folder: '{path_to_archive}'"
            .format(**self.configs)
        )


    #%%
    def get_engine(self):
        """
        Get an SQLAlchemy engine that can be used to connect to the database.
        """
        if self.configs['db_type'] == 'sqlite':
            info_message = (
                "Connected to a {db_type} database in the following location: {sqlite_db_path}"
                .format(**self.configs)
            )
            # set `isolation_level` to 'READ UNCOMMITTED' so that reads are non-blocking
            # (required for SCINET)
            engine = sa.create_engine(
                '{db_type}:///{sqlite_db_path}'.format(**self.configs),
                 isolation_level='READ UNCOMMITTED'
            )
            enable_sqlite_foreign_key_checks(engine)
        elif self.configs['db_type'] in ['postgresql', 'mysql']:
            info_message = (
                "Connected to a {db_type} database: {db_url}:{db_port}/{db_database}{db_socket}"
                .format(**self.configs)
            )
            engine = sa.create_engine(
                '{db_type}://{db_username}:{db_password}@{db_url}:{db_port}/{db_database}{db_socket}'
                .format(**self.configs)
            ) # echo=True
        else:
            raise Exception("Unsupported `db_type`: '{}'".format(self.configs['db_type']))
        self.logger.info(info_message)
        return engine
    
    
    def configure_session(self):
        """
        Configure the Session class to use the current engine.
        
        `autocommit` and `autoflush` are enabled for the `sqlite` database in order to improve
        performance.
        """
        if self.configs['db_type'] == 'sqlite':
            autocommit = False #True
            autoflush = False #True
            retry_on_failure = True
        elif self.configs['db_type'] in ['postgresql', 'mysql']:
            autocommit = False
            autoflush = False
            retry_on_failure = True
        Session.configure(bind=self.engine, autocommit=autocommit, autoflush=autoflush)
        if retry_on_failure:
            decorator = (
                decorate_all_methods(
                    retry(retry_on_exception=lambda exc: isinstance(exc, sa_exc.OperationalError),
                          wait_exponential_multiplier=1000, # start with one second delay
                          wait_exponential_max=600000) # go up to 10 minutes
                )
            )
            return decorator(Session)
        else:
            return Session


    @contextmanager
    def session_scope(self):
        """ 
        Provide a transactional scope around a series of operations.
        Enables the following construct: ``with self.session_scope() as session:``.
        """
        session = self.Session()
        try:
            yield session
            if not self.configs['db_is_immutable']:
                session.commit()
        except:
            if not self.configs['db_is_immutable']:
                session.rollback()
            raise
        finally:
            session.close()
            
            
    #%%    
    def create_database_tables(self, clear_schema=False, keep_uniprot_sequence=True):
        """
        Create a new database in the schema specified by the ``schema_version`` global variable.
        If ``clear_schema == True``, remove all the tables in the schema first.
        
        DANGER!!! 
        Using this function with an existing database can lead to loss of data.
        Make sure that you know what you are doing.
        
        Parameters
        ----------
        clear_schema : bool
            Whether or not to delete all tables in the database schema before creating new tables.
        keep_uniprot_sequence : bool
            Whether or not to keep the `uniprot_sequence` table. 
            Only relevant if `clear_schema` is `True`.
        """
        # 
        if clear_schema:
            self.delete_database_tables(drop_schema=False, keep_uniprot_sequence=keep_uniprot_sequence)
            self.logger.debug('Database schema was cleared successfully.')
        
        # Create all tables, creating schema as neccessary
        for table in Base.metadata.sorted_tables:
            try:
                table.create(self.engine, checkfirst=True)
            except (sa.exc.OperationalError, sa.exc.ProgrammingError) as e:
                self.logger.error(str(e))
                if re.search('schema .* does not exist', str(e)):
                    missing_schema = str(e)[str(e).find('schema "')+8:str(e).find('" does not exist')]
                elif 'Unknown database ' in str(e):
                    missing_schema = str(e)[str(e).find("Unknown database '")+18:str(e).find("'\")")]
                else:
                    raise e
                sql_command = 'create schema {};'.format(missing_schema)
                self.logger.warning("Creating missing schema with system command: '{}'".format(sql_command))
                self.engine.execute(sql_command)
                table.create(self.engine, checkfirst=True)
        self.logger.debug('Database tables were created successfully.')
    
    
    def delete_database_tables(self, drop_schema=False, keep_uniprot_sequence=True):
        """
        Parameters
        ----------
        drop_schema : bool
            Whether or not to drop the schema after dropping the tables.
        keep_uniprot_sequence : bool
            Wheter or not to keep the table (and schema) containing uniprot sequences.
        """       
        if self.configs['db_type'] == 'sqlite':
            os.remove(self.configs['sqlite_db_path'])
            self.logger.info("Successfully removed the sqlite database file: {sqlite_db_path}".format(**self.configs))
            return
        
        # Remove tables one by one
        for table in reversed(Base.metadata.sorted_tables):
            if table.name != 'uniprot_sequence':
                self.configs['table_name'] = table.name
                self.engine.execute('drop table {db_schema}.{table_name};'.format(**self.configs))
            elif not keep_uniprot_sequence:
                self.configs['table_name'] = table.name
                self.engine.execute('drop table {db_schema_uniprot}.{table_name};'.format(**self.configs))
                
        # Remove the database schema
        uniprot_on_diff_schema = self.configs['db_schema'] != self.configs['db_schema_uniprot']
        if drop_schema and uniprot_on_diff_schema:
            self.engine.execute('drop schema {db_schema};'.format(**self.configs))
        if drop_schema and uniprot_on_diff_schema and not keep_uniprot_sequence:
            self.engine.execute('drop schema {db_schema_uniprot};'.format(**self.configs))
        if drop_schema and not uniprot_on_diff_schema and not keep_uniprot_sequence:
            self.engine.execute('drop schema {db_schema};'.format(**self.configs))
        
        self.logger.info("Successfully removed the {db_type} database schema: {db_schema}".format(**self.configs))    



    #%%
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
        if self.configs['debug']:
            self.logger.debug(system_command)
        result, error_message, return_code = hf.popen(system_command)
        if return_code != 0:
            self.logger.error(result)
            raise Exception(error_message)
        
        
    def copy_table_to_db(self, table_name, table_folder):
        """
        Copy data from a ``.tsv`` file to a table in the database.
        """
        configs = self.configs.copy()
        configs['table_name'] = table_name
        configs['table_folder'] = table_folder        
        
        def _format_configs():
            if table_name == 'uniprot_sequence':
                configs['table_db_schema'] = configs['db_schema_uniprot']
            else:
                configs['table_db_schema'] = configs['db_schema']
        
        self.logger.info("Copying '{table_name}' to '{db_type}' database...".format(**configs))
        if configs['db_type'] == 'sqlite':
            self._load_data_into_sqlite(configs)
        elif configs['db_type'] == 'mysql':
            _format_configs()
            configs['sql_command'] = self.mysql_command_template.format(**configs)
            system_command = self.mysql_load_table_template.format(**configs)
            self._run_create_table_system_command(system_command)
        elif configs['db_type'] == 'postgresql':
            _format_configs()
            configs['sql_command'] = self.psql_command_template.format(**configs)
            system_command = self.psql_load_table_template.format(**configs)
            self._run_create_table_system_command(system_command)
        else:
            raise Exception("Unsupported database type: '{}'".format(configs['db_type']))



    #%% Get objects from the database
    def get_rows_by_ids(self, row_object, row_object_identifiers, row_object_identifier_values):
        """ Get the rows from the table *row_object* identified by keys
        *row_object_identifiers* with values *row_object_identifier_values*
        """
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


    def get_domain(self, pfam_names, subdomains=False):
        """ 
        Returns pdbfam-based definitions of all pfam domains in the pdb.
        """
        with self.session_scope() as session:
            domain_set = set()
            for pfam_name in pfam_names:
                if not subdomains:
                    domain = (
                        session.query(Domain)
                        .filter(Domain.pfam_name==pfam_name)
                        .distinct().all() )
                else:
                    domain = (
                        session.query(Domain).filter(
                            (Domain.pfam_name.like(pfam_name)) |
                            (Domain.pfam_name.like(pfam_name+'+%')) |
                            (Domain.pfam_name.like(pfam_name+'\_%')) | # need an escape character because _ matches any single character
                            (Domain.pfam_name.like('%+'+pfam_name)) |
                            (Domain.pfam_name.like('%+'+pfam_name+'+%')) |
                            (Domain.pfam_name.like('%+'+pfam_name+'\_%')) ) # need an escape character because _ matches any single character
                        .distinct().all() )
                domain_set.update(domain)
        if not domain_set:
            self.logger.debug('No domain definitions found for pfam: %s' % str(pfam_names))
        return list(domain_set)


    def get_domain_contact(self, pfam_names_1, pfam_names_2, subdomains=False):
        """ 
        Returns domain-domain interaction information from pdbfam.
        Note that the produced dataframe may not have the same order as the keys.
        """
        with self.session_scope() as session:
            domain_contact_1 = self._get_domain_contact(pfam_names_1, pfam_names_2, session, subdomains)
            domain_contact_2 = self._get_domain_contact(pfam_names_2, pfam_names_1, session, subdomains)

        if not len(domain_contact_1) and not len(domain_contact_2):
            self.logger.debug('No domain contact template found for domains %s, %s' % (str(pfam_names_1), str(pfam_names_2),))

        return [domain_contact_1, domain_contact_2]


    def _get_domain_contact(self, pfam_names_1, pfam_names_2, session, subdomains):
        """
        """
        domain_1 = sa.orm.aliased(Domain)
        domain_2 = sa.orm.aliased(Domain)
        domain_contact_set = set()
        for pfam_name_1 in pfam_names_1:
            for pfam_name_2 in pfam_names_2:
                if not subdomains:
                    domain_contact = (
                        session.query(DomainContact)
                        # .join(domain_1, DomainContact.cath_id_1==domain_1.cath_id)
                        .filter(domain_1.pfam_name==pfam_name_1)
                        # .join(domain_2, DomainContact.cath_id_2==domain_2.cath_id)
                        .filter(domain_2.pfam_name==pfam_name_2)
                        .distinct().all() )
                else:
                    domain_contact = (
                        session.query(DomainContact)
                        # .join(domain_1, DomainContact.cath_id_1==domain_1.cath_id)
                        .filter(
                            (domain_1.pfam_name.like(pfam_name_1)) |
                            (domain_1.pfam_name.like(pfam_name_1+'+%')) |
                            (domain_1.pfam_name.like(pfam_name_1+'\_%')) | # need an escape character because _ matches any single character
                            (domain_1.pfam_name.like('%+'+pfam_name_1)) |
                            (domain_1.pfam_name.like('%+'+pfam_name_1+'+%')) |
                            (domain_1.pfam_name.like('%+'+pfam_name_1+'\_%')) ) # need an escape character because _ matches any single character
                        # .join(domain_2, DomainContact.cath_id_2==domain_2.cath_id)
                        .filter(
                            (domain_2.pfam_name.like(pfam_name_2)) |
                            (domain_2.pfam_name.like(pfam_name_2+'+%')) |
                            (domain_2.pfam_name.like(pfam_name_2+'\_%')) | # need an escape character because _ matches any single character
                            (domain_2.pfam_name.like('%+'+pfam_name_2)) |
                            (domain_2.pfam_name.like('%+'+pfam_name_2+'+%')) |
                            (domain_2.pfam_name.like('%+'+pfam_name_2+'\_%')) ) # need an escape character because _ matches any single character
                        .distinct().all() )
                domain_contact_set.update(domain_contact)
        return list(domain_contact_set)


    def get_uniprot_domain(self, uniprot_id, copy_data=False):
        """
        """
        with self.session_scope() as session:
            uniprot_domains = (
                session
                .query(UniprotDomain)
                .filter(UniprotDomain.uniprot_id == uniprot_id)
                # .options(sa.orm.joinedload('template').sa.orm.joinedload('model'))
                .options(sa.orm.joinedload('template', innerjoin=True))
                .all()
            )

        path_to_archive = self.configs['path_to_archive']
        d_idx = 0
        while d_idx < len(uniprot_domains):
            d = uniprot_domains[d_idx]
            if not d.template:
                self.logger.debug(
                    'Skipping uniprot domain with id {} because it does not '
                    'have a structural template...'
                    .format(d.uniprot_domain_id))
                del uniprot_domains[d_idx]
                continue
            # Copy precalculated Provean data
            if copy_data:
                try:
                    self._copy_provean(d, path_to_archive)
                except subprocess.CalledProcessError as e:
                    self.logger.error(e)
                    self.logger.error('Failed to copy provean supporting set!')
                    d.uniprot_sequence.provean.provean_supset_filename = ''
            # Copy precalculated homology models
            if copy_data:
                try:
                    self._copy_uniprot_domain_data(d, d.path_to_data, path_to_archive)
                except subprocess.CalledProcessError as e:
                    self.logger.error(e)
                    self.logger.error('Failed to copy the domain alignment and / or homology model!')
                    d.template.model.alignment_filename = None
                    d.template.model.model_filename = None
            d_idx += 1

        return uniprot_domains


    def get_uniprot_domain_pair(self, uniprot_id, copy_data=False, uniprot_domain_pair_ids=[]):
        """
        """
        with self.session_scope() as session:
            uniprot_domain_pairs_query = (
                session.query(UniprotDomainPair)
                .filter(sa.or_(
                    sa.text("uniprot_domain_1.uniprot_id='{}'".format(uniprot_id)),
                    sa.text("uniprot_domain_2.uniprot_id='{}'".format(uniprot_id))))
            )
            if uniprot_domain_pair_ids:
                uniprot_domain_pairs_query = (
                    uniprot_domain_pairs_query
                    .filter(UniprotDomainPair.uniprot_domain_pair_id.in_(uniprot_domain_pair_ids))
                )
            uniprot_domain_pairs = (
                uniprot_domain_pairs_query
                # .options(sa.orm.joinedload('template').sa.orm.joinedload('model'))
                .options(sa.orm.joinedload('template', innerjoin=True))
                .all()
            )

        path_to_archive = self.configs['path_to_archive']
        d_idx = 0
        while d_idx < len(uniprot_domain_pairs):
            d = uniprot_domain_pairs[d_idx]
            if not d.template:
                self.logger.debug(
                    'Skipping uniprot domain pair with id {} because it does not '
                    'have a structural template...'
                    .format(d.uniprot_domain_pair_id))
                del uniprot_domain_pairs[d_idx]
                continue
            # Copy precalculated Provean data
            if copy_data:
                if d.uniprot_domain_1.uniprot_id == uniprot_id:
                    ud = d.uniprot_domain_1
                elif d.uniprot_domain_2.uniprot_id == uniprot_id:
                    ud = d.uniprot_domain_2
                try:
                    self._copy_provean(ud, path_to_archive)
                except subprocess.CalledProcessError as e:
                    self.logger.error(e)
                    self.logger.error('Failed to copy provean supporting set!')
                    d.uniprot_sequence.provean.provean_supset_filename = ''
            # Copy precalculated homology models
            if copy_data:
                try:
                    self._copy_uniprot_domain_pair_data(d, d.path_to_data, path_to_archive)
                except subprocess.CalledProcessError as e:
                    self.logger.error(e)
                    self.logger.error('Failed to copy domain pair alignments and / or homology model!')
                    d.template.model.alignment_filename_1 = None
                    d.template.model.alignment_filename_2 = None
                    d.template.model.model_filename = None
            d_idx += 1

        return uniprot_domain_pairs


    def _copy_uniprot_domain_data(self, d, path_to_data, path_to_archive):
        if path_to_data is None:
            self.logger.error('Cannot copy uniprot domain data because `path_to_data` is None')
            return
        if (d.template != None and
            d.template.model != None and
            d.template.model.alignment_filename != None and
            d.template.model.model_filename != None):
                if path_to_archive.endswith('.7z'):
                    # Extract files from a 7zip archive
                    filenames = [
                        path_to_data + d.template.model.alignment_filename,
                        path_to_data + d.template.model.model_filename,
                    ]
                    self._extract_files_from_7zip(path_to_archive, filenames)
                else:
                    tmp_save_path = self.configs['temp_path'] + path_to_data
                    archive_save_path = path_to_archive + path_to_data
                    path_to_alignment = tmp_save_path + '/'.join(d.template.model.alignment_filename.split('/')[:-1]) + '/'
                    subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(path_to_alignment), shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        archive_save_path + d.template.model.alignment_filename,
                        tmp_save_path + d.template.model.alignment_filename), shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        archive_save_path + d.template.model.model_filename,
                        tmp_save_path + d.template.model.model_filename), shell=True)


    def _copy_uniprot_domain_pair_data(self, d, path_to_data, path_to_archive):
        if path_to_data is None:
            self.logger.error('Cannot copy uniprot domain pair data because `path_to_data` is None')
            return
        if (d.template != None and
            d.template.model != None and
            d.template.model.alignment_filename_1 != None and
            d.template.model.alignment_filename_2 != None and
            d.template.model.model_filename != None):
                if path_to_archive.endswith('.7z'):
                    # Extract files from a 7zip archive
                    filenames = [
                        path_to_data + d.template.model.alignment_filename_1,
                        path_to_data + d.template.model.alignment_filename_2,
                        path_to_data + d.template.model.model_filename,
                    ]
                    self._extract_files_from_7zip(path_to_archive, filenames)
                else:
                    tmp_save_path = self.configs['temp_path'] + path_to_data
                    archive_save_path = path_to_archive + path_to_data
                    path_to_alignment_1 = tmp_save_path + '/'.join(d.template.model.alignment_filename_1.split('/')[:-1]) + '/'
                    path_to_alignment_2 = tmp_save_path + '/'.join(d.template.model.alignment_filename_2.split('/')[:-1]) + '/'
                    subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(path_to_alignment_1), shell=True)
                    subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(path_to_alignment_2), shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        archive_save_path + d.template.model.alignment_filename_1,
                        tmp_save_path + d.template.model.alignment_filename_1), shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        archive_save_path + d.template.model.alignment_filename_2,
                        tmp_save_path + d.template.model.alignment_filename_2), shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        archive_save_path + d.template.model.model_filename,
                        tmp_save_path + d.template.model.model_filename), shell=True)


    def _copy_provean(self, ud, path_to_archive):
        if (ud.uniprot_sequence and
            ud.uniprot_sequence.provean and
            ud.uniprot_sequence.provean.provean_supset_filename):
                if path_to_archive.endswith('.7z'):
                    # Extract files from a 7zip archive
                    filenames = [
                        get_uniprot_base_path(ud) + ud.uniprot_sequence.provean.provean_supset_filename,
                        get_uniprot_base_path(ud) + ud.uniprot_sequence.provean.provean_supset_filename + '.fasta',
                    ]
                    self._extract_files_from_7zip(path_to_archive, filenames)
                else:
                    # Compy files from the archive folders
                    subprocess.check_call(
                        "umask ugo=rwx; mkdir -m 777 -p '{}'".format(
                            os.path.dirname(
                                self.configs['temp_path'] + get_uniprot_base_path(ud) +
                                ud.uniprot_sequence.provean.provean_supset_filename)),
                        shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        path_to_archive + get_uniprot_base_path(ud) +
                            ud.uniprot_sequence.provean.provean_supset_filename,
                        self.configs['temp_path'] + get_uniprot_base_path(ud) +
                            ud.uniprot_sequence.provean.provean_supset_filename), shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        path_to_archive + get_uniprot_base_path(ud) +
                            ud.uniprot_sequence.provean.provean_supset_filename + '.fasta',
                        self.configs['temp_path'] + get_uniprot_base_path(ud) +
                            ud.uniprot_sequence.provean.provean_supset_filename + '.fasta'), shell=True)


    @retry(retry_on_exception=lambda exc: type(exc) == error.Archive7zipError, # isinstance(exc, error.Archive7zipError)
           wait_exponential_multiplier=1000, 
           wait_exponential_max=60000)
    def _extract_files_from_7zip(self, path_to_archive, filenames):
        """Extract files to `config['temp_path']`
        """
        system_command = "7za x '{path_to_archive}' '{files}' -y".format(
            path_to_archive=path_to_archive, 
            files="' '".join(filenames)
        )
        self.logger.debug(
            'Extracting files from 7zip archive using the following system command:\n{}'
            .format(system_command))
        result, error_message, return_code = (
            hf.run_subprocess_locally_full(self.configs['temp_path'], system_command)
        )
        
        def log_error():
            self.logger.error(
                '\n result:{}\n error_message:{}\n return_code:{}'
                .format(result, error_message, return_code)
            )
            
        if 'No files to process' in result:
            log_error()
            raise error.Archive7zipFileNotFoundError(result, error_message, return_code)
        
        if return_code:
            log_error()
            raise error.Archive7zipError(result, error_message, return_code)


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
                        .scalar() )
        elif isinstance(d, UniprotDomainPair) and isinstance(uniprot_id, six.string_types):
            with self.session_scope() as session:
                uniprot_mutation = (
                    session.query(UniprotDomainPairMutation)
                        .filter(
                            (UniprotDomainPairMutation.uniprot_id == uniprot_id) &
                            (UniprotDomainPairMutation.uniprot_domain_pair_id == d.uniprot_domain_pair_id) &
                            (UniprotDomainPairMutation.mutation == mutation))
                        .scalar() )
        else:
            self.logger.debug('d: {}\td type: {}'.format(d, type(d)))
            raise Exception('Not enough arguments, or the argument types are incorrect!')

        if copy_data:
            try:
                self._copy_mutation_data(uniprot_mutation, d.path_to_data, self.configs['path_to_archive'])
            except (subprocess.CalledProcessError,
                    error.Archive7zipError,
                    error.Archive7zipFileNotFoundError) as e:
                self.logger.error(e)
                uniprot_mutation.model_filename_wt = None
        return uniprot_mutation


    def _copy_mutation_data(self, mutation, path_to_data, path_to_archive):
        if mutation and mutation.model_filename_wt:
            if path_to_archive.endswith('.7z'):
                # Extract files from a 7zip archive
                filenames = [
                    path_to_data + mutation.model_filename_wt,
                    path_to_data + mutation.model_filename_mut,
                ]
                self._extract_files_from_7zip(path_to_archive, filenames)
            else:
                tmp_save_path = self.configs['temp_path'] + path_to_data
                archive_save_path = path_to_archive + path_to_data
                path_to_mutation = os.path.dirname(tmp_save_path + mutation.model_filename_wt)
                subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(path_to_mutation), shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    archive_save_path + mutation.model_filename_wt,
                    tmp_save_path + mutation.model_filename_wt), shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    archive_save_path + mutation.model_filename_mut,
                    tmp_save_path + mutation.model_filename_mut), shell=True)


    def remove_model(self, d):
        """
        """
        if isinstance(d, UniprotDomain):
            with self.session_scope() as session:
                session.execute(
                    'delete from {0}.uniprot_domain_model where uniprot_domain_id = {1}'
                    .format(self.configs['db_schema'], d.uniprot_domain_id))
        elif isinstance(d, UniprotDomainPair):
            with self.session_scope() as session:
                session.execute(
                    'delete from {0}.uniprot_domain_pair_model where uniprot_domain_pair_id = {1}'
                    .format(self.configs['db_schema'], d.uniprot_domain_pair_id))
        else:
            raise Exception('Not enough arguments, or the argument types are incorrect!')


    #%% Add objects to the database
    def merge_row(self, row_instance):
        """Adds a list of rows (`row_instances`) to the database.
        """
        if not self.configs['db_is_immutable']:
            with self.session_scope() as session:
                if not isinstance(row_instance, list):
                    session.merge(row_instance)
                else:
                    deque( (session.merge(row) for row in row_instance), maxlen=0 )


    def merge_provean(self, provean, uniprot_base_path=False):
        """Adds provean score to the database.
        """
        if (uniprot_base_path and
            provean.provean_supset_filename and
            os.path.isfile(self.configs['temp_path'] + uniprot_base_path +
                provean.provean_supset_filename) and
            os.path.isfile(self.configs['temp_path'] + uniprot_base_path +
                provean.provean_supset_filename + '.fasta')):
            # ...
            path_to_archive = self.configs['path_to_archive']
            self.logger.debug(
                'Moving provean supset to the output folder: {}'
                .format(path_to_archive + uniprot_base_path))
            subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(
                path_to_archive + uniprot_base_path), shell=True)
            subprocess.check_call("cp -f '{}' '{}'".format(
                self.configs['temp_path'] + uniprot_base_path + provean.provean_supset_filename,
                path_to_archive + uniprot_base_path + provean.provean_supset_filename), shell=True)
            subprocess.check_call("cp -f '{}' '{}'".format(
                self.configs['temp_path'] + uniprot_base_path + provean.provean_supset_filename + '.fasta',
                path_to_archive + uniprot_base_path + provean.provean_supset_filename + '.fasta'), shell=True)
        self.merge_row(provean)


    def merge_model(self, d, path_to_data=False):
        """Adds MODELLER models to the database.
        """
        # Save a copy of the alignment to the export folder
        if path_to_data:
            path_to_archive = self.configs['path_to_archive']
            archive_save_path = path_to_archive + path_to_data
            tmp_save_path = self.configs['temp_path'] + path_to_data
            # Save the row corresponding to the model as a serialized sqlalchemy object
            subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(archive_save_path), shell=True)
            # Don't need to dump template. Templates are precalculated
            # pickle_dump(sa_ext_serializer.dumps(d.template), archive_save_path + 'template.pickle')
            # pickle_dump(sa_ext_serializer.dumps(d.template.model), archive_save_path + 'model.pickle')
            # Save the modelled structure
            if d.template.model.model_filename is not None:
                # Save alignments
                if isinstance(d.template.model, UniprotDomainModel):
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        tmp_save_path + d.template.model.alignment_filename,
                        archive_save_path + d.template.model.alignment_filename), shell=True)
                elif isinstance(d.template.model, UniprotDomainPairModel):
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        tmp_save_path + d.template.model.alignment_filename_1,
                        archive_save_path + d.template.model.alignment_filename_1), shell=True)
                    subprocess.check_call("cp -f '{}' '{}'".format(
                        tmp_save_path + d.template.model.alignment_filename_2,
                        archive_save_path + d.template.model.alignment_filename_2), shell=True)
                # Save the model
                subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(archive_save_path), shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    tmp_save_path + d.template.model.model_filename,
                    archive_save_path + d.template.model.model_filename), shell=True)
        self.merge_row([d.template, d.template.model])


    def merge_mutation(self, mut, path_to_data=False):
        """
        """
        mut.mut_date_modified = datetime.datetime.utcnow()
        if path_to_data and (mut.model_filename_wt is not None):
            path_to_archive = self.configs['path_to_archive']
            archive_save_path = path_to_archive + path_to_data
            tmp_save_path = self.configs['temp_path'] + path_to_data
            archive_save_subpath = mut.model_filename_wt.split('/')[0] + '/'
            # Save the row corresponding to the mutation as a serialized sqlalchemy object
            subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(
                archive_save_path + archive_save_subpath), shell=True)
            # pickle_dump(sa_ext_serializer.dumps(mut), archive_save_path + archive_save_subpath + 'mutation.pickle')
            if mut.model_filename_wt and mut.model_filename_mut:
                # Save Foldx structures
                subprocess.check_call("cp -f '{}' '{}'".format(
                    tmp_save_path + mut.model_filename_wt,
                    archive_save_path + mut.model_filename_wt), shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    tmp_save_path + mut.model_filename_mut,
                    archive_save_path + mut.model_filename_mut), shell=True)
        self.merge_row(mut)


    #%%
    def get_uniprot_sequence(self, uniprot_id, check_external=True):
        """        
        Parameters
        ----------
        uniprot_id : str
            Uniprot ID of the protein
        check_external : bool
            Whether or not to look online if the protein sequence is not found in the local database
        
        Returns
        -------
        SeqRecord :
            Contains the sequence of the specified uniprot
        """
        with self.session_scope() as session:
            uniprot_sequence = session\
                .query(UniprotSequence)\
                .filter(UniprotSequence.uniprot_id==uniprot_id)\
                .all()

        if len(uniprot_sequence) == 1:
            uniprot_sequence = uniprot_sequence[0]

        elif len(uniprot_sequence) > 1:
            self.logger.error('Type(uniprot_sequence): {}'.format(type(uniprot_sequence)))
            self.logger.error('uniprot_sequence: {}'.format(type(uniprot_sequence)))
            raise Exception('Several uniprot sequences returned!? This should never happen!')

        elif len(uniprot_sequence) == 0:
            username = hf.get_username()
            if (username.strip() == 'joan' # on Scinet
                    or not check_external): # don't bother with external sequences
                print (
                    "Couldn't find a sequence for uniprot {}, and not bothering to look for it online"
                    .format(uniprot_id))
                return None
            else:
                self.logger.debug('Fetching sequence for uniprot {} from an online server'.format(uniprot_id))
                print('Fetching sequence for uniprot {} from an online server'.format(uniprot_id))
                address = 'http://www.uniprot.org/uniprot/{}.fasta'.format(uniprot_id)
                try:
                    handle = urllib.request.urlopen(address)
                    sequence = next(SeqIO.parse(handle, "fasta"))
                except (StopIteration, urllib.error.HTTPError) as e:
                    self.logger.debug('{}: {}'.format(type(e), str(e)))
                    print('{}: {}'.format(type(e), str(e)))
                    return None
                uniprot_sequence = UniprotSequence()
                uniprot_sequence.uniprot_id = uniprot_id
                sp_or_trembl, uniprot_id_2, uniprot_name = sequence.name.split('|')
                if uniprot_id != uniprot_id_2:
                    print (
                        'Uniprot id of the fasta file ({}) does not match the '
                        'uniprot id of the query ({}). Skipping...'
                        .format(uniprot_id_2, uniprot_id))
                    return None
                uniprot_sequence.db = sp_or_trembl
                uniprot_sequence.uniprot_name = uniprot_name
                uniprot_sequence.uniprot_description = sequence.description
                uniprot_sequence.uniprot_sequence = str(sequence.seq)
                self.add_uniprot_sequence(uniprot_sequence)

        uniprot_seqrecord = SeqRecord(
            seq=Seq(str(uniprot_sequence.uniprot_sequence)),
            id=uniprot_sequence.uniprot_id,
            name=uniprot_sequence.uniprot_name)

        return uniprot_seqrecord


    def add_uniprot_sequence(self, uniprot_sequence):
        """ Add new sequences to the database.
        :param uniprot_sequence: UniprotSequence object
        :rtype: None
        """
        with self.session_scope() as session:
            session.add(uniprot_sequence)


    #%% Domain and domain contact
    def add_domain(self, d):
        with self.session_scope() as session:
            if isinstance(d, Domain):
                (
                    session
                    .query(Domain)
                    .filter(Domain.cath_id == d.cath_id)
                    .update({Domain.domain_errors: d.domain_errors})
                )
            elif isinstance(d, DomainContact):
                (
                    session
                    .query(DomainContact)
                    .filter(DomainContact.domain_contact_id == d.domain_contact_id)
                    .update({DomainContact.domain_contact_errors: d.domain_contact_errors})
                )


    def add_domain_errors(self, t, error_string):
        with self.session_scope() as session:
            if isinstance(t, UniprotDomain):
                domain = (
                    session
                    .query(Domain)
                    .filter(Domain.cath_id==t.cath_id)
                    .as_scalar()
                )
                domain.domain_errors = error_string
                session.merge(domain)
            elif isinstance(t, UniprotDomainPair):
                domain_contact = (
                    session
                    .query(DomainContact)
                    .filter(DomainContact.cath_id_1==t.cath_id_1)
                    .filter(DomainContact.cath_id_2==t.cath_id_2)
                    .all()[0]
                )
                domain_contact.domain_contact_errors = error_string
                session.merge(domain_contact)
            else:
                raise Exception('Wrong type for template!!!')



    #%%
    def _split_domain(self, domain):
        """
        Converts a string of domain boundaries to integers.
        The separator is '-', and it can happen that both or one boundary is negative, 
        making this tricky::

            -150-200,   meaning from -150 to 200
            -150--100,  meaning from -150 to -100

        Parameters
        ----------
        domain : str
            A string of domain boundaries
        
        Returns
        -------
        list
            Domain boundaries converted to integers
        """
        # split the domain boundaries, keep eventual minus signs
        if domain[0] == '-' and len(domain[1:].split('-')) == 2:
            domain = ['-' + domain[1:].split('-')[0], domain[1:].split('-')[1]]
        elif domain[0] == '-' and len(domain[1:].split('-')) > 2:
            domain = ['-' + domain[1:].split('-')[0], '-' + domain[1:].split('-')[-1]]
        else:
            domain = [domain.split('-')[0], domain.split('-')[1]]
        # strip the letters
        if domain[0][-1] in hf.uppercase:
            domain[0] = domain[0][:-1]
        if domain[1][-1] in hf.uppercase:
            domain[1] = domain[1][:-1]
        domain = [int(domain[0]), int(domain[1])]
        return domain


    def _split_domain_semicolon(self, domains):
        """ 
        Unlike `split_domain()`, this function returns a tuple of tuples of strings,
        preserving letter numbering (e.g. 10B).
        """
        x = domains
        return tuple([ tuple([r.strip() for r in ro.split(':')]) for ro in x.split(',') ])


    def _split_interface_aa(self, interface_aa):
        """
        """
        if interface_aa and (interface_aa != '') and (interface_aa != 'NULL'):
            if interface_aa[-1] == ',':
                interface_aa = interface_aa[:-1]

            x  = interface_aa
            return_tuple = tuple([int(r.strip()) for r in x.split(',')])

        else:
            return_tuple = []

        return return_tuple


    #%%
    def get_alignment(self, model, path_to_data):
        """
        """

        tmp_save_path = self.configs['temp_path'] + path_to_data
        archive_save_path = self.configs['path_to_archive'] + path_to_data

        if isinstance(model, UniprotDomainModel):

            # Load previously-calculated alignments
            if os.path.isfile(tmp_save_path + model.alignment_filename):
                alignment = AlignIO.read(tmp_save_path + model.alignment_filename, 'clustal')
            elif os.path.isfile(archive_save_path + model.alignment_filename):
                alignment = AlignIO.read(archive_save_path + model.alignment_filename, 'clustal')
            else:
                raise error.NoPrecalculatedAlignmentFound(archive_save_path, model.alignment_filename)

            return [alignment, None]

        elif isinstance(model, UniprotDomainPairModel):

            # Read alignment from the temporary folder
            if (os.path.isfile(tmp_save_path + model.alignment_filename_1)
            and os.path.isfile(tmp_save_path + model.alignment_filename_2)):
                alignment_1 = AlignIO.read(tmp_save_path + model.alignment_filename_1, 'clustal')
                alignment_2 = AlignIO.read(tmp_save_path + model.alignment_filename_2, 'clustal')
            # Read alignment from the export database
            elif (os.path.isfile(archive_save_path + model.alignment_filename_1)
            and os.path.isfile(archive_save_path + model.alignment_filename_2)):
                alignment_1 = AlignIO.read(archive_save_path + model.alignment_filename_1, 'clustal')
                alignment_2 = AlignIO.read(archive_save_path + model.alignment_filename_2, 'clustal')
            else:
                raise error.NoPrecalculatedAlignmentFound(archive_save_path, model.alignment_filename_1)

            return [alignment_1, alignment_2]


    def load_db_from_archive(self):
        """
        TODO: In the future I should move back to using json...
        """
        data = [
            ['human/*/*/*/*/template.json', UniprotDomainTemplate],
            ['human/*/*/*/*/model.json', UniprotDomainModel],
            ['human/*/*/*/*/*/mutation.json', UniprotDomainMutation],
            ['human/*/*/*/*/*/*/template.json', UniprotDomainPairTemplate],
            ['human/*/*/*/*/*/*/model.json', UniprotDomainPairModel],
            ['human/*/*/*/*/*/*/*/mutation.json', UniprotDomainPairMutation],
        ]

        for d in data:
            result, __, __ = hf.popen('ls ' + self.configs['path_to_archive'] + d[0])
            filenames = [fname for fname in result.split('\n') if fname != '']
            for filename in filenames:
                with open(filename, 'r') as fh:
                    row = json.load(fh)
                try:
                    self.session.merge(d[1](**row))
                except TypeError as e:
                    self.logger.debug(
                        'Error merging {}.\n'
                        'Probably from an older version of the database.\n'
                        'Skipping...'.format(filename))
                    self.logger.debug('\t', e)
                self.logger.debug('Merged %s' % filename)
            self.session.commit()
            self.logger.debug('Committed changes\n\n\n')



#%% Elaspic-specific helper functions
def get_uniprot_base_path(d):
    """ The uniprot id is cut into several chunks to create folders that will
    hold a manageable number of pdbs.
    """
    if isinstance(d, UniprotDomain):
        uniprot_id = d.uniprot_id
        uniprot_name = d.uniprot_sequence.uniprot_name
    elif isinstance(d, UniprotDomainPair):
        uniprot_id = d.uniprot_domain_1.uniprot_id
        uniprot_name = d.uniprot_domain_1.uniprot_sequence.uniprot_name
    elif isinstance(d, dict):
        uniprot_id = d['uniprot_id']
        uniprot_name = d['uniprot_name']
    else:
        raise Exception('Input parameter type is not supported!')

    uniprot_base_path = (
        '{organism_name}/{uniprot_id_part_1}/{uniprot_id_part_2}/{uniprot_id_part_3}/'
        .format(
            organism_name=uniprot_name.split('_')[-1].lower(),
            uniprot_id_part_1=uniprot_id[:3],
            uniprot_id_part_2=uniprot_id[3:5],
            uniprot_id_part_3=uniprot_id,))
    return uniprot_base_path


def get_uniprot_domain_path(d):
    """ Return the path to individual domains or domain pairs.
    """
    if isinstance(d, UniprotDomain):
        uniprot_domain_path = (
            '{pfam_clan:.36}.{alignment_def}/'
            .format(
                pfam_clan=d.pfam_clan,
                alignment_def=d.alignment_def.replace(':','-'),))
    elif isinstance(d, UniprotDomainPair):
        uniprot_domain_path = (
            '{pfam_clan_1:.36}.{alignment_def_1}/{pfam_clan_2:.36}.{alignment_def_2}/{uniprot_id_2}/'
            .format(
                pfam_clan_1 = d.uniprot_domain_1.pfam_clan,
                alignment_def_1 = d.uniprot_domain_1.alignment_def.replace(':','-'),
                pfam_clan_2 = d.uniprot_domain_2.pfam_clan,
                alignment_def_2 = d.uniprot_domain_2.alignment_def.replace(':','-'),
                uniprot_id_2 = d.uniprot_domain_2.uniprot_id,))
    return uniprot_domain_path


def scinet_cleanup(folder, destination, name=None):
    """
    zip and copy the results from the ramdisk to /scratch
    """
    print('saving the result in', folder)
    os.chdir(folder)
    if name == None:
        output_name = 'result_' + time.strftime("%Y_%m_%d_at_%Hh_%Mm") + '.tar.bz2'
    else:
        output_name = name + '_' + time.strftime("%Y_%m_%d_at_%Hh_%Mm") + '.tar.bz2'
    copy = 'cp ' + output_name + ' ' + destination + output_name
#    copy = 'cp ' + output_name + ' $SCRATCH/niklas-pipeline/' + output_name
#    copy = 'cp ' + output_name + ' /home/niklas/tmp/' + output_name
    system_command = 'tar -acf ' + output_name + ' * && ' + copy

    result, error_message, return_code = hf.popen(system_command)
    if return_code != 0:
        print('error_message:', error_message)
    return


def pickle_dump(obj, filename):
    mode = stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH | stat.S_IWOTH
    umask_original = os.umask(0)
    try:
        handle = os.fdopen(os.open(filename, os.O_WRONLY | os.O_CREAT, mode), 'wb')
    finally:
        os.umask(umask_original)
    pickle.dump(obj, handle, pickle.HIGHEST_PROTOCOL)
    handle.close()
