# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 18:58:50 2012

@author: kimlab
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import next
from builtins import object

import os
import stat
import urllib.request
import urllib.error
import subprocess
import json
import datetime
import time
import pickle
from contextlib import contextmanager
from collections import deque

import six
from sqlalchemy import or_
from sqlalchemy import create_engine
from sqlalchemy import Column, Index
from sqlalchemy import Integer, Float, String, Text, DateTime
from sqlalchemy import ForeignKey
from sqlalchemy.orm import sessionmaker, relationship, backref, aliased, joinedload
from sqlalchemy.ext.declarative import declarative_base, DeferredReflection
from sqlalchemy.ext.serializer import dumps

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#import parse_pfamscan
from . import conf
from . import helper_functions as hf
from . import errors as error


#%% Constants

# Default sizes for creating varchar fields
SHORT = 15
MEDIUM = 255
LONG = 16384

naming_convention = {
  "ix": 'ix_%(column_0_label)s',
  "uq": "uq_%(table_name)s_%(column_0_name)s",
  "ck": "ck_%(table_name)s_%(constraint_name)s",
  "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
  "pk": "pk_%(table_name)s"
}


#%% Some database-specific parameters that SQLAlchemy can't figure out

DB_TYPE = conf.configs.get('db_type')
DB_DATABASE = conf.configs.get('db_database')
DB_SCHEMA = conf.configs.get('db_schema')
DB_SCHEMA_UNIPROT = conf.configs.get('db_schema_uniprot')

db_specific_properties = {
    'mysql': {
        'BINARY_COLLATION': 'utf8_bin',
        'STRING_COLLATION': 'utf8_unicode_ci',
        'schema_version_tuple': {'schema': DB_SCHEMA},
        'uniprot_kb_schema_tuple': {'schema': DB_SCHEMA_UNIPROT},
    },
    'postgresql': {
        'BINARY_COLLATION': 'en_US.utf8',
        'STRING_COLLATION': 'en_US.utf8',
        'schema_version_tuple': {'schema': DB_SCHEMA},
        'uniprot_kb_schema_tuple': {'schema': DB_SCHEMA_UNIPROT},
    },
    'sqlite': {
        'BINARY_COLLATION': 'RTRIM',
        'STRING_COLLATION': 'NOCASE',
        'schema_version_tuple': {'sqlite_autoincrement': True},
        'uniprot_kb_schema_tuple': {'sqlite_autoincrement': True},
    }, 
}


def get_db_specific_param(key):
    if DB_TYPE is None:
        error_message = (
            'The `DB_TYPE` has not been set. Do not know what database is being used!'
        )
        raise Exception(error_message)
    if (DB_TYPE in ['mysql', 'postgresql'] and 
        (DB_DATABASE is None or DB_SCHEMA is None or DB_SCHEMA_UNIPROT is None)):
            error_message = (
                'Both the `DB_SCHEMA` and `DB_SCHEMA_UNIPROT` have to be specified when using '
                'a MySQL or PostgreSQL database!'
            )
            raise Exception(error_message)
    return db_specific_properties[DB_TYPE][key]
    

def get_table_args(table_name, index_columns=[], db_specific_params=[]):
    """
    Returns a tuple of additional table arguments.
    """
    table_args = []
    # Create indexes over several columns
    for columns in index_columns:
        if type(columns) == tuple:
            column_names, unique = columns
        elif type(columns) == list:
            column_names = columns
            unique = False
        index_name = (
            'ix_{table_name}_{column_0_name}'
            .format(table_name=table_name, column_0_name=column_names[0])[:255]
        )
        table_args.append(Index(index_name, *column_names, unique=unique))
    # Other table parameters, such as schemas, etc.
    for db_specific_param in db_specific_params:
        table_args.append(get_db_specific_param(db_specific_param))
    return tuple(table_args)
   


#%%
Base = declarative_base()
Base.metadata.naming_conventions = naming_convention

class Domain(Base):
    """ Table containing pdbfam domain definitions for all pdbs()
    """
    __tablename__ = 'domain'
    _indexes = [
        (['pdb_id', 'pdb_chain', 'pdb_pdbfam_name', 'pdb_pdbfam_idx'], True),
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
    
    cath_id = Column(
        String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')), 
        primary_key=True)
    pdb_id = Column(String(SHORT), nullable=False)
    pdb_chain = Column(String(SHORT), nullable=False)
    pdb_domain_def = Column(String(MEDIUM), nullable=False)
    pdb_pdbfam_name = Column(String(LONG), nullable=False, index=True)
    pdb_pdbfam_idx = Column(Integer)
    domain_errors = Column(Text)


class DomainContact(Base):
    """ Table containing interactions between all pdbfam domains in the pdb
    """
    __tablename__ = 'domain_contact'
    _indexes = [
        (['cath_id_1', 'cath_id_2'], True),
        (['cath_id_2', 'cath_id_1'], True),
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
    
    domain_contact_id = Column(Integer, primary_key=True)
    cath_id_1 = Column(
        None, ForeignKey(Domain.cath_id, onupdate='cascade', ondelete='cascade'), nullable=False)
    cath_id_2 = Column(
        None, ForeignKey(Domain.cath_id, onupdate='cascade', ondelete='cascade'), nullable=False)
#    cath_id_2 = Column(
#        String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')), 
#        nullable=False)
    min_interchain_distance = Column(Float)
    contact_volume = Column(Float)
    contact_surface_area = Column(Float)
    atom_count_1 = Column(Integer)
    atom_count_2 = Column(Integer)
    number_of_contact_residues_1 = Column(Integer)
    number_of_contact_residues_2 = Column(Integer)
    contact_residues_1 = Column(Text)
    contact_residues_2 = Column(Text)
    crystal_packing = Column(Float)
    domain_contact_errors = Column(Text)

    # Relationships
    domain_1 = relationship(
        Domain, primaryjoin=cath_id_1==Domain.cath_id, cascade='expunge', lazy='joined')
#    # the second domain may be a ligand or a peptide, and so the foreign key constraint does not work
    domain_2 = relationship(
        Domain, primaryjoin=cath_id_2==Domain.cath_id, cascade='expunge', lazy='joined')


class UniprotSequence(Base):
    """ Table containing the entire Swissprot + Trembl database as well as any
    additional sequences that were added to the database.
    """
    __tablename__ = 'uniprot_sequence'
    __table_args__ = get_table_args(__tablename__, [], ['uniprot_kb_schema_tuple'])

    db = Column(String(SHORT), nullable=False)
    uniprot_id = Column(String(SHORT), primary_key=True)
    uniprot_name = Column(String(SHORT), nullable=False)
    protein_name = Column(String(MEDIUM))
    organism_name = Column(String(MEDIUM), index=True)
    gene_name = Column(String(MEDIUM), index=True)
    protein_existence = Column(Integer)
    sequence_version = Column(Integer)
    uniprot_sequence = Column(Text, nullable=False)


class Provean(Base):
    __tablename__ = 'provean'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_id = Column(
        None, ForeignKey(
            UniprotSequence.uniprot_id, 
            onupdate='cascade', ondelete='cascade'), 
        primary_key=True)
    provean_supset_filename = Column(String(MEDIUM))
    provean_supset_length = Column(Integer)
    provean_errors = Column(Text)
    provean_date_modified = Column(
        DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)

    # Relationships
    uniprot_sequence = relationship(
        UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('provean', uselist=False, cascade='expunge', lazy='joined'))


class UniprotDomain(Base):
    __tablename__ = 'uniprot_domain'
    if 'training' in DB_SCHEMA:
        # The database used for storing training data has an extra column `max_seq_identity`,
        # because we want to make homology models at different sequence identities.
        max_seq_identity = Column(Integer)
        _indexes = [
            (['uniprot_id', 'alignment_def', 'max_seq_identity'], True),
        ]
        _create_uniprot_id_index = False
    else:
        _indexes = []
        _create_uniprot_id_index = True
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
    
    uniprot_domain_id = Column(Integer, nullable=False, primary_key=True, autoincrement=True)
    uniprot_id = Column(
        None, ForeignKey(
            UniprotSequence.uniprot_id, 
            onupdate='cascade', ondelete='cascade'),
        index=_create_uniprot_id_index, nullable=False)
    pdbfam_name = Column(String(LONG), index=True, nullable=False)
    pdbfam_idx = Column(Integer, nullable=False)
    pfam_clan = Column(Text)
    alignment_def = Column(String(MEDIUM))
    pfam_names = Column(String(LONG))
    alignment_subdefs = Column(Text)
    path_to_data = Column(Text)
    
    # Relationships
    uniprot_sequence = relationship(
        UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('uniprot_domain', cascade='expunge')) # many to one



class UniprotDomainPair(Base):
    __tablename__ = 'uniprot_domain_pair'
    _indexes = [
            (['uniprot_domain_id_1', 'uniprot_domain_id_2'], True),
            (['uniprot_domain_id_2', 'uniprot_domain_id_1'], True),
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_domain_pair_id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_domain_id_1 = Column(
        None, ForeignKey(
            UniprotDomain.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    uniprot_domain_id_2 = Column(
        None, ForeignKey(
            UniprotDomain.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    rigids = Column(Text) # Interaction references from iRefIndex
    domain_contact_ids = Column(Text) # interaction references from PDBfam
    path_to_data = Column(Text)

    # Relationships
    uniprot_domain_1 = relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_1==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one
    uniprot_domain_2 = relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_2==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one


class UniprotDomainTemplate(Base):
    __tablename__ = 'uniprot_domain_template'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_domain_id = Column(
        None, ForeignKey(
            UniprotDomain.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    template_errors = Column(Text)
    cath_id = Column(
        None, ForeignKey(
            Domain.cath_id, 
            onupdate='cascade', ondelete='cascade'), 
        index=True, nullable=False)
    domain_start = Column(Integer, index=True)
    domain_end = Column(Integer, index=True)
    domain_def = Column(String(MEDIUM))
    alignment_identity = Column(Float)
    alignment_coverage = Column(Float)
    alignment_score = Column(Float)
    t_date_modified = Column(
        DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)
        
    # Relationships
    uniprot_domain = relationship(
        UniprotDomain, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('template', uselist=False, cascade='expunge', lazy='joined')) # one to one
    domain = relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('uniprot_domain', cascade='expunge')) # many to one



class UniprotDomainModel(Base):
    __tablename__ = 'uniprot_domain_model'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_domain_id = Column(
        None, ForeignKey(
            UniprotDomainTemplate.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    model_errors = Column(Text)
    alignment_filename = Column(String(MEDIUM))
    model_filename = Column(String(MEDIUM))
    chain = Column(String(SHORT))
    norm_dope = Column(Float)
    sasa_score = Column(Text)
    model_domain_def = Column(String(MEDIUM))
    m_date_modified = Column(DateTime, default=datetime.datetime.utcnow,
                             onupdate=datetime.datetime.utcnow, nullable=False)
                             
    # Relationships
    template = relationship(
        UniprotDomainTemplate, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('model', uselist=False, cascade='expunge', lazy='joined')) # one to one



class UniprotDomainMutation(Base):
    __tablename__ = 'uniprot_domain_mutation'
    _indexes = [
        ['uniprot_id', 'mutation'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_id = Column(
        None, ForeignKey(
            UniprotSequence.uniprot_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    uniprot_domain_id = Column(
        None, ForeignKey(
            UniprotDomainModel.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True, index=True)
    mutation = Column(String(SHORT), index=True, nullable=False, primary_key=True)
    mutation_errors = Column(Text)
    model_filename_wt = Column(String(MEDIUM))
    model_filename_mut = Column(String(MEDIUM))
    chain_modeller = Column(String(SHORT))
    mutation_modeller = Column(String(SHORT))
    stability_energy_wt = Column(Text)
    stability_energy_mut = Column(Text)
    physchem_wt = Column(Text)
    physchem_wt_ownchain = Column(Text)
    physchem_mut = Column(Text)
    physchem_mut_ownchain = Column(Text)
    matrix_score = Column(Float)
    secondary_structure_wt = Column(Text)
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(Text)
    solvent_accessibility_mut = Column(Float)
    provean_score = Column(Float)
    ddg = Column(Float, index=True)
    mut_date_modified = Column(
        DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)
        
    # Relationships
    model = relationship(
        UniprotDomainModel, cascade='expunge', uselist=False, lazy='joined',
        backref=backref('mutations', cascade='expunge')) # many to one



class UniprotDomainPairTemplate(Base):
    __tablename__ = 'uniprot_domain_pair_template'
    _indexes = [
        ['cath_id_1', 'cath_id_2'],
        ['cath_id_2', 'cath_id_1'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_domain_pair_id = Column(
        None, ForeignKey(
            UniprotDomainPair.uniprot_domain_pair_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    domain_contact_id = Column(
        None, ForeignKey(
            DomainContact.domain_contact_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False)
    cath_id_1 = Column(
        None, ForeignKey(
            Domain.cath_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    cath_id_2 = Column(
        None, ForeignKey(
            Domain.cath_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)

    identical_1 = Column(Float)
    conserved_1 = Column(Float)
    coverage_1 = Column(Float)
    score_1 = Column(Float)

    identical_if_1 = Column(Float)
    conserved_if_1 = Column(Float)
    coverage_if_1 = Column(Float)
    score_if_1 = Column(Float)

    identical_2 = Column(Float)
    conserved_2 = Column(Float)
    coverage_2 = Column(Float)
    score_2 = Column(Float)

    identical_if_2 = Column(Float)
    conserved_if_2 = Column(Float)
    coverage_if_2 = Column(Float)
    score_if_2 = Column(Float)

    score_total = Column(Float)
    score_if_total = Column(Float)
    score_overall = Column(Float)

    t_date_modified = Column(
        DateTime, default=datetime.datetime.utcnow,
        onupdate=datetime.datetime.utcnow, nullable=False)
    template_errors = Column(Text)

    # Relationships
    domain_pair = relationship(
        UniprotDomainPair, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('template', uselist=False, cascade='expunge', lazy='joined')) # one to one
    domain_contact = relationship(
        DomainContact, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('uniprot', cascade='expunge')) # one to one
    domain_1 = relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        primaryjoin=(cath_id_1==Domain.cath_id)) # many to one
    domain_2 = relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        primaryjoin=(cath_id_2==Domain.cath_id)) # many to one



class UniprotDomainPairModel(Base):
    __tablename__ = 'uniprot_domain_pair_model'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_domain_pair_id = Column(
        None, ForeignKey(
            UniprotDomainPairTemplate.uniprot_domain_pair_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    model_errors = Column(Text)
    alignment_filename_1 = Column(String(MEDIUM))
    alignment_filename_2 = Column(String(MEDIUM))
    model_filename = Column(String(MEDIUM))
    chain_1 = Column(String(SHORT))
    chain_2 = Column(String(SHORT))
    norm_dope = Column(Float)
    interface_area_hydrophobic = Column(Float)
    interface_area_hydrophilic = Column(Float)
    interface_area_total = Column(Float)
    interface_dg = Column(Float)
    interacting_aa_1 = Column(Text)
    interacting_aa_2 = Column(Text)
    model_domain_def_1 = Column(String(MEDIUM))
    model_domain_def_2 = Column(String(MEDIUM))
    m_date_modified = Column(
        DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)
        
    # Relationships
    template = relationship(
        UniprotDomainPairTemplate, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('model', uselist=False, cascade='expunge', lazy='joined')) # one to one
        


class UniprotDomainPairMutation(Base):
    __tablename__ = 'uniprot_domain_pair_mutation'
    _indexes = [
        ['uniprot_id', 'mutation'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_id = Column(None, ForeignKey(
        UniprotSequence.uniprot_id, onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    uniprot_domain_pair_id = Column(None, ForeignKey(
        UniprotDomainPairModel.uniprot_domain_pair_id, onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    mutation = Column(String(SHORT),
        nullable=False, primary_key=True)
    mutation_errors = Column(Text)
    model_filename_wt = Column(String(MEDIUM))
    model_filename_mut = Column(String(MEDIUM))
    chain_modeller = Column(String(SHORT))
    mutation_modeller = Column(String(SHORT))
    analyse_complex_energy_wt = Column(Text)
    stability_energy_wt = Column(Text)
    analyse_complex_energy_mut = Column(Text)
    stability_energy_mut = Column(Text)
    physchem_wt = Column(Text)
    physchem_wt_ownchain = Column(Text)
    physchem_mut = Column(Text)
    physchem_mut_ownchain = Column(Text)
    matrix_score = Column(Float)
    secondary_structure_wt = Column(Text)
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(Text)
    solvent_accessibility_mut = Column(Float)
    contact_distance_wt = Column(Float)
    contact_distance_mut = Column(Float)
    provean_score = Column(Float)
    ddg = Column(Float, index=False)
    mut_date_modified = Column(DateTime, default=datetime.datetime.utcnow,
                               onupdate=datetime.datetime.utcnow, nullable=False)
    # Relationships
    model = relationship(
        UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('mutations', cascade='expunge')) # many to one



#%%
# Get the session that will be used for all future queries
# Expire on commit so that you keep all the table objects even after the session closes.
Session = sessionmaker(expire_on_commit=False)
#Session = scoped_session(sessionmaker(expire_on_commit=False))

class MyDatabase(object):
    """
    """
    
    def __init__(self, configs=conf.configs, logger=None):
        """
        """
        if logger is None:
            self.logger = hf.get_logger()
        else:
            self.logger = logger
            
        self.schema_version = configs['db_schema']
        
        # Choose which database to use
        if configs['db_type'] == 'sqlite':
            self.logger.info(
                "Connecting to an {dt_type} database in the following location: {sqlite_db_path}..."
                .format(**configs)
            )
            autocommit = True
            autoflush = True
            engine = create_engine(
                '{db_type}:///{sqlite_db_path}'.format(**configs),
                isolation_level='READ UNCOMMITTED'
            )
        elif configs['db_type'] in ['postgresql', 'mysql']:
            self.logger.info(
                "Connecting to a {db_type} database ..."
                .format(**configs)
            )
            autocommit = False
            autoflush = False
            engine = create_engine(
                '{db_type}://{db_username}:{db_password}@{db_url}:{db_port}/{db_schema}'
                .format(**configs)
            ) # echo=True
        else:
            raise Exception("Only mysql, postgresql, and sqlite are currently supported!")
                      
        Session.configure(bind=engine, autocommit=autocommit, autoflush=autoflush)
        self.Session = Session
        self.engine = engine
        self.autocommit = autocommit
        self.db_is_immutable = configs['db_is_immutable']
        self.temp_path = configs['temp_path']
        self.path_to_archive = configs['path_to_archive']
        
        self.logger.info(
            "Using precalculated data from the following folder: '{path_to_archive}'"
            .format(**configs)
        )


    @contextmanager
    def session_scope(self):
        """ Provide a transactional scope around a series of operations.
        So you can use: `with self.session_scope() as session:`
        """
        session = self.Session()
        try:
            yield session
            if not self.db_is_immutable:
                session.commit()
        except:
            if not self.db_is_immutable:
                session.rollback()
            raise
        finally:
            session.close()

    
    def create_database_tables(self, clear_schema=False, keep_uniprot_sequence=True):
        """
        Create a new database in the schema specified by the ``schema_version`` global variable.
        If ``clear_schema`` == True, remove all the tables in the schema first.
        
        DANGER!!! 
        Using this function with an existing database can lead to loss of data.
        Make sure that you know what you are doing.
        
        Parameters
        ----------
        engine : sa.Engine
            SQLAlchemy engine connected to the database of interest.
        clear_schema : bool
            Whether or not to drop all the tables in the database before creating new tables.
        keep_uniprot_sequence : bool
            Whether or not to keep the `uniprot_sequence` table. 
            Only relevant if `clear_schema` is `True`.
        """
        metadata_tables = Base.metadata.tables.copy()
        if keep_uniprot_sequence:
            uniprot_sequence_table = [c for c in metadata_tables.keys() if 'uniprot_sequence' in c]
            assert len(uniprot_sequence_table) == 1
            del metadata_tables[uniprot_sequence_table[0]]
        if clear_schema:
            Base.metadata.drop_all(self.engine, metadata_tables.values())
            self.logger.debug('Database schema was cleared successfully.')
        Base.metadata.create_all(self.engine, metadata_tables.values())
        self.logger.debug('Database schema was created successfully.')



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
        domain_1 = aliased(Domain)
        domain_2 = aliased(Domain)
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
                # .options(joinedload('template').joinedload('model'))
                .options(joinedload('template', innerjoin=True))
                .all()
            )

        d_idx = 0
        while d_idx < len(uniprot_domains):
            d = uniprot_domains[d_idx]
            if not d.template:
                self.logger.debug(
                    'Skipping uniprot domain with id {} because it does not '
                    'have a structural template...'.format(d.uniprot_domain_id))
                del uniprot_domains[d_idx]
                continue
            if copy_data:
                try:
                    self._copy_uniprot_domain_data(d, d.path_to_data)
                except subprocess.CalledProcessError as e:
                    self.logger.error(e)
                    d.template.model.alignment_filename = None
                    d.template.model.model_filename = None
            d_idx += 1

        return uniprot_domains


    def get_uniprot_domain_pair(self, uniprot_id, copy_data=False):
        """
        """
        with self.session_scope() as session:
            uniprot_domain_pairs = (
                session.query(UniprotDomainPair)
                .filter(or_(
                    "uniprot_domain_1.uniprot_id='{}'".format(uniprot_id),
                    "uniprot_domain_2.uniprot_id='{}'".format(uniprot_id)))
                # .options(joinedload('template').joinedload('model'))
                .options(joinedload('template', innerjoin=True))
                .all()
            )

        d_idx = 0
        while d_idx < len(uniprot_domain_pairs):
            d = uniprot_domain_pairs[d_idx]
            if not d.template:
                self.logger.debug(
                    'Skipping uniprot domain pair with id {} because it does not '
                    'have a structural template...'.format(d.uniprot_domain_pair_id))
                del uniprot_domain_pairs[d_idx]
                continue
            if copy_data:
                try:
                    self._copy_uniprot_domain_pair_data(d, d.path_to_data, uniprot_id)
                except subprocess.CalledProcessError as e:
                    self.logger.error(e)
                    d.template.model.alignment_filename_1 = None
                    d.template.model.alignment_filename_2 = None
                    d.template.model.model_filename = None
            d_idx += 1

        return uniprot_domain_pairs


    def _copy_uniprot_domain_data(self, d, path_to_data):
        if path_to_data is None:
            self.logger.error('Cannot copy uniprot domain data because `path_to_data` is None')
            return
        if (d.template != None and
            d.template.model != None and
            d.template.model.alignment_filename != None and
            d.template.model.model_filename != None):
                tmp_save_path = self.temp_path + path_to_data
                archive_save_path = self.path_to_archive + path_to_data
                path_to_alignment = tmp_save_path + '/'.join(d.template.model.alignment_filename.split('/')[:-1]) + '/'
                subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(path_to_alignment), shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    archive_save_path + d.template.model.alignment_filename,
                    tmp_save_path + d.template.model.alignment_filename), shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    archive_save_path + d.template.model.model_filename,
                    tmp_save_path + d.template.model.model_filename), shell=True)
        # Copy Provean supporting set
        try:
            self._copy_provean(d)
        except subprocess.CalledProcessError as e:
            self.logger.error('Failed to copy provean supporting set!')
            self.logger.error(e)
            d.uniprot_sequence.provean.provean_supset_filename = ''


    def _copy_uniprot_domain_pair_data(self, d, path_to_data, uniprot_id):
        if path_to_data is None:
            self.logger.error('Cannot copy uniprot domain pair data because `path_to_data` is None')
            return
        if (d.template != None and
            d.template.model != None and
            d.template.model.alignment_filename_1 != None and
            d.template.model.alignment_filename_2 != None and
            d.template.model.model_filename != None):
                tmp_save_path = self.temp_path + path_to_data
                archive_save_path = self.path_to_archive + path_to_data
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
        # Copy Provean supporting set
        if d.uniprot_domain_1.uniprot_id == uniprot_id:
            self._copy_provean(d.uniprot_domain_1)
        elif d.uniprot_domain_2.uniprot_id == uniprot_id:
            self._copy_provean(d.uniprot_domain_2)


    def _copy_provean(self, ud):
        if (ud.uniprot_sequence and
            ud.uniprot_sequence.provean and
            ud.uniprot_sequence.provean.provean_supset_filename):
                subprocess.check_call(
                    "umask ugo=rwx; mkdir -m 777 -p '{}'".format(
                        os.path.dirname(
                            self.temp_path + get_uniprot_base_path(ud) +
                            ud.uniprot_sequence.provean.provean_supset_filename)),
                    shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    self.path_to_archive + get_uniprot_base_path(ud) +
                        ud.uniprot_sequence.provean.provean_supset_filename,
                    self.temp_path + get_uniprot_base_path(ud) +
                        ud.uniprot_sequence.provean.provean_supset_filename), shell=True)
                subprocess.check_call("cp -f '{}' '{}'".format(
                    self.path_to_archive + get_uniprot_base_path(ud) +
                        ud.uniprot_sequence.provean.provean_supset_filename + '.fasta',
                    self.temp_path + get_uniprot_base_path(ud) +
                        ud.uniprot_sequence.provean.provean_supset_filename + '.fasta'), shell=True)


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
        elif isinstance(d, UniprotDomainPair) and isinstance(uniprot_id, str):
            with self.session_scope() as session:
                uniprot_mutation = (
                    session.query(UniprotDomainPairMutation)
                        .filter(
                            (UniprotDomainPairMutation.uniprot_id == uniprot_id) &
                            (UniprotDomainPairMutation.uniprot_domain_pair_id == d.uniprot_domain_pair_id) &
                            (UniprotDomainPairMutation.mutation == mutation))
                        .scalar() )
        else:
            raise Exception('Not enough arguments, or the argument types are incorrect!')

        if copy_data:
            try:
                self._copy_mutation_data(uniprot_mutation, d.path_to_data)
            except subprocess.CalledProcessError as e:
                self.logger.error(e)
                uniprot_mutation.model_filename_wt = None
        return uniprot_mutation


    def _copy_mutation_data(self, mutation, path_to_data):
        if mutation and mutation.model_filename_wt:
            tmp_save_path = self.temp_path + path_to_data
            archive_save_path = self.path_to_archive + path_to_data
            path_to_mutation = tmp_save_path + '/'.join(mutation.model_filename_wt.split('/')[:-1]) + '/'
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
                    .format(self.schema_version, d.uniprot_domain_id))
        elif isinstance(d, UniprotDomainPair):
            with self.session_scope() as session:
                session.execute(
                    'delete from {0}.uniprot_domain_pair_model where uniprot_domain_pair_id = {1}'
                    .format(self.schema_version, d.uniprot_domain_pair_id))
        else:
            raise Exception('Not enough arguments, or the argument types are incorrect!')


    #%% Add objects to the database
    def merge_row(self, row_instance):
        """Adds a list of rows (`row_instances`) to the database.
        """
        if not self.db_is_immutable:
            with self.session_scope() as session:
                if not isinstance(row_instance, list):
                    session.merge(row_instance)
                else:
                    deque( (session.merge(row) for row in row_instance), maxlen=0 )


    def merge_provean(self, provean, uniprot_base_path):
        """Adds provean score to the database.
        """
        if (provean.provean_supset_filename and
                os.path.isfile(self.temp_path + uniprot_base_path +
                    provean.provean_supset_filename) and
                os.path.isfile(self.temp_path + uniprot_base_path +
                    provean.provean_supset_filename + '.fasta') ):
            self.logger.debug('Moving provean supset to the output folder: {}'.format(self.path_to_archive + uniprot_base_path))
            subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(self.path_to_archive + uniprot_base_path), shell=True)
            subprocess.check_call("cp -f '{}' '{}'".format(
                self.temp_path + uniprot_base_path + provean.provean_supset_filename,
                self.path_to_archive + uniprot_base_path + provean.provean_supset_filename), shell=True)
            subprocess.check_call("cp -f '{}' '{}'".format(
                self.temp_path + uniprot_base_path + provean.provean_supset_filename + '.fasta',
                self.path_to_archive + uniprot_base_path + provean.provean_supset_filename + '.fasta'), shell=True)
        self.merge_row(provean)


    def merge_model(self, d, path_to_data=False):
        """Adds MODELLER models to the database.
        """
        # Save a copy of the alignment to the export folder
        if path_to_data:
            tmp_save_path = self.temp_path + path_to_data
            archive_save_path = self.path_to_archive + path_to_data
            # Save the row corresponding to the model as a serialized sqlalchemy object
            subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(archive_save_path), shell=True)
            # Don't need to dump template. Templates are precalculated
            # pickle_dump(dumps(d.template), archive_save_path + 'template.pickle')
            pickle_dump(dumps(d.template.model), archive_save_path + 'model.pickle')
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
            tmp_save_path = self.temp_path + path_to_data
            archive_save_path = self.path_to_archive + path_to_data
            archive_save_subpath = mut.model_filename_wt.split('/')[0] + '/'
            # Save the row corresponding to the mutation as a serialized sqlalchemy object
            subprocess.check_call("umask ugo=rwx; mkdir -m 777 -p '{}'".format(
                archive_save_path + archive_save_subpath), shell=True)
            pickle_dump(dumps(mut), archive_save_path + archive_save_subpath + 'mutation.pickle')
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

        tmp_save_path = self.temp_path + path_to_data
        archive_save_path = self.path_to_archive + path_to_data

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
            childProcess = hf.popen_py2i3('ls ' + self.path_to_archive + d[0])
            result, __ = childProcess.communicate()
            if six.PY3:
                result = str(result, encoding='utf-8')
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

    childProcess = hf.popen_py2i3(system_command)
    result, error_message = childProcess.communicate()
    if six.PY3:
        result = str(result, encoding='utf-8')
        error_message = str(error_message, encoding='utf-8')
    if childProcess.returncode != 0:
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
