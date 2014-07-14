# -*- coding: utf-8 -*-
"""
Created on Sun Feb  3 15:07:51 2013

@author: alexey
"""
import os
import pandas as pd
import urllib2
import subprocess
import json
from string import uppercase
import datetime
import logging
from contextlib import contextmanager
import cPickle as pickle

from sqlalchemy import create_engine
from sqlalchemy import Column, Index, UniqueConstraint
from sqlalchemy import Integer, Float, String, Boolean, Text, DateTime, Sequence
from sqlalchemy import ForeignKey
from sqlalchemy.orm import sessionmaker, relationship, backref, aliased, scoped_session, joinedload
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.serializer import loads, dumps

from Bio import Seq
from Bio import SeqIO
from Bio import AlignIO

#import parse_pfamscan
import helper_functions as hf
import errors as error


###############################################################################
# Not the best place to define this, bot collation changes depending on the
# database type...
#sql_flavor = 'sqlite_file'
sql_flavor = 'postgresql'

if sql_flavor.split('_')[0] == 'sqlite':
    binary_collation = 'RTRIM' # same as binary, except that trailing space characters are ignored.
    string_collation = 'NOCASE'
elif sql_flavor.split('_')[0] == 'mysql':
    binary_collation = 'utf8_bin'
    string_collation = 'utf8_unicode_ci'
elif sql_flavor.split('_')[0] == 'postgresql':
    binary_collation = 'en_US.utf8'
    string_collation = 'en_US.utf8'
else:
    raise Exception('Unknown database type!')

Base = declarative_base()

# Get the session that will be used for all future queries
# Expire on commit so that you keep all the table objects even after the
# session closes.
Session = sessionmaker(expire_on_commit=False)
#Session = scoped_session(sessionmaker(expire_on_commit=False))


###############################################################################
class Domain(Base):
    """ Table containing pdbfam domain definitions for all pdbs
    """
    __tablename__ = 'domain'
    __table_args__ = ({'schema': 'elaspic_v2'},)

    cath_id = Column(String(15, collation=binary_collation), primary_key=True)
    pdb_id = Column(String, nullable=False)
    pdb_type = Column(String, nullable=True)
    pdb_resolution = Column(Float, nullable=True)
    pdb_chain = Column(String, nullable=False)
    pdb_domain_def = Column(String, nullable=False)
    pdb_pdbfam_name = Column(String, nullable=False)
    pdb_pdbfam_idx = Column(Integer)
    domain_errors = Column(String)


class DomainContact(Base):
    """ Table containing interactions between all pdbfam domains in the pdb
    """
    __tablename__ = 'domain_contact'
    __table_args__ = (
        Index('cath_id_1_2', 'cath_id_1', 'cath_id_1', unique=True),
        {'sqlite_autoincrement': True, 'schema': 'elaspic_v2'},
    )

    domain_contact_id = Column(Integer, Sequence('domain_contact_domain_contact_id_seq'), primary_key=True)
    cath_id_1 = Column(None, ForeignKey(Domain.cath_id), index=True, nullable=False)
    cath_id_2 = Column(None, ForeignKey(Domain.cath_id), index=True, nullable=False)
    min_interchain_distance = Column(Float)
    contact_volume = Column(Float)
    contact_surface_area = Column(Float)
    atom_count_1 = Column(Integer)
    atom_count_2 = Column(Integer)
    contact_residues_1 = Column(Text)
    contact_residues_2 = Column(Text)
    crystal_packing = Column(Float)
    domain_contact_errors = Column(String)

    # Relationships
    domain_1 = relationship(Domain, primaryjoin=cath_id_1==Domain.cath_id, cascade='expunge', lazy='joined')
    domain_2 = relationship(Domain, primaryjoin=cath_id_2==Domain.cath_id, cascade='expunge', lazy='joined')


class UniprotSequence(Base):
    """ Table containing the entire Swissprot + Trembl database as well as any
    additional sequences that were added to the database.
    """
    __tablename__ = 'uniprot_sequence'
    __table_args__ = ({'schema': 'elaspic_v2'},)

    origin = Column(String, nullable=False)
    uniprot_id = Column(String, primary_key=True, nullable=False)
    uniprot_name = Column(String, nullable=False)
    protein_name = Column(String)
    organism_name = Column(String)
    gene_name = Column(String)
    protein_existence = Column(Integer)
    sequence_version = Column(Integer)
    seq = Column(Text, nullable=False)


class UniprotDomain(Base):
    __tablename__ = 'uniprot_domain'
    __table_args__ = (
        UniqueConstraint('uniprot_id', 'pdbfam_name', 'alignment_def', name='unique_uniprot_domain'),
        {'schema': 'elaspic_v2'},
    )

    uniprot_domain_id = Column(Integer, Sequence('uniprot_domain_uniprot_domain_id_seq'), primary_key=True)
    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False)
    uniprot_name = Column(Text)
    pdbfam_name = Column(String, index=True, nullable=False)
    pdbfam_idx = Column(Integer, nullable=False)
    pfam_clan = Column(String)
    alignment_def = Column(String)
    pfam_names = Column(String)
    alignment_subdefs = Column(String)
    cath_id = Column(String)
    pdb_id = Column(String)
    pdb_chain = Column(String)
    pdb_pdbfam_name = Column(String)
    pdb_pdbfam_idx = Column(Integer)
    domain_start = Column(Integer)
    domain_end = Column(Integer)
    domain_def = Column(String)
    alignment_identity = Column(Float)
    alignment_coverage = Column(Float)
    alignment_score = Column(Float)
    path_to_data = Column(Text)
    uniprot_domain_errors = Column(Text)

    # Relationships
    uniprot_sequence = relationship(
        UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('domain', cascade='expunge')) # many to one


class UniprotProvean(Base):
    __tablename__ = 'uniprot_provean'
    __table_args__ = ({'schema': 'elaspic_v2'})

    uniprot_id = Column(None, ForeignKey(UniprotDomain.uniprot_id), primary_key=True)
    provean_supset_filename = Column(String)
    provean_supset_length = Column(Integer)
    provean_errors = Column(Text)
    provean_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    uniprot_domains = relationship(
        UniprotDomain, cascade='expunge', lazy='joined',
        backref=backref('uniprot_provean', uselist=False, cascade='expunge', lazy='joined'))


class UniprotDomainPair(Base):
    __tablename__ = 'uniprot_domain_pair'
    __table_args__ = (
        Index('uniprot_domain_id_1_2', 'uniprot_domain_id_1', 'uniprot_domain_id_2', unique=True),
        {'sqlite_autoincrement': True, 'schema': 'elaspic_v2'},
    )

    uniprot_domain_pair_id = Column(Integer, Sequence('uniprot_domain_pair_id_sequence'), primary_key=True)
    uniprot_domain_id_1 = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), index=True, nullable=False)
    uniprot_id_1 = Column(String)
    uniprot_domain_id_2 = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), index=True, nullable=False)
    uniprot_id_2 = Column(String)
    path_to_data = Column(Text)
    rigids = Column(Text)
    domain_contact_ids = Column(Text)
    pdb_id = Column(String)
    cath_id_1 = Column(String)
    pdb_chain_1 = Column(String)
    pdb_pdbfam_name_1 = Column(String)
    pdb_pdbfam_idx_1 = Column(Integer)
    domain_def_1 = Column(String)
    alignment_identity_1 = Column(Float)
    alignment_coverage_1 = Column(Float)
    alignment_score_1 = Column(Float)
    alignment_if_identity_1 = Column(Float)
    alignment_if_coverage_1 = Column(Float)
    cath_id_2 = Column(String)
    pdb_chain_2 = Column(String)
    pdb_pdbfam_name_2 = Column(String)
    pdb_pdbfam_idx_2 = Column(Integer)
    domain_def_2 = Column(String)
    alignment_identity_2 = Column(Float)
    alignment_coverage_2 = Column(Float)
    alignment_score_2 = Column(Float)
    alignment_if_identity_2 = Column(Float)
    alignment_if_coverage_2 = Column(Float)

    # Relationships
    uniprot_domain_1 = relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_1==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one
    uniprot_domain_2 = relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_2==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one


class UniprotDomainModel(Base):
    __tablename__ = 'uniprot_domain_model'
    __table_args__ = ({'schema': 'elaspic_v2'},)

    uniprot_domain_id = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), primary_key=True)
    model_errors = Column(Text)
    model_filename = Column(String)
    chain = Column(String)
    norm_dope = Column(Float)
    sasa_score = Column(Text)
    m_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    uniprot_domain = relationship(
        UniprotDomain, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('model', uselist=False, cascade='expunge')) # one to one


class UniprotDomainMutation(Base):
    __tablename__ = 'uniprot_domain_mutation'
    __table_args__ = ({'schema': 'elaspic_v2'},)

    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False, primary_key=True)
    uniprot_domain_id = Column(None, ForeignKey(UniprotDomainModel.uniprot_domain_id), primary_key=True)
    mutation = Column(String, nullable=False, primary_key=True)
    mutation_errors = Column(Text)
    model_filename_wt = Column(String)
    model_filename_mut = Column(String)
    chain_modeller = Column(String)
    mutation_modeller = Column(String)
    stability_energy_wt = Column(Text)
    stability_energy_mut = Column(Text)
    physchem_wt = Column(Text)
    physchem_wt_ownchain = Column(Text)
    physchem_mut = Column(Text)
    physchem_mut_ownchain = Column(Text)
    matrix_score = Column(Float)
    secondary_structure_wt = Column(String)
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(String)
    solvent_accessibility_mut = Column(Float)
    provean_score = Column(Float)
    ddg = Column(Float)
    mut_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    model = relationship(
        UniprotDomainModel, cascade='expunge', uselist=False, lazy='joined',
        backref=backref('mutations', cascade='expunge')) # many to one


class UniprotDomainPairModel(Base):
    __tablename__ = 'uniprot_domain_pair_model'
    __table_args__ = ({'schema': 'elaspic_v2'},)

    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPair.uniprot_domain_pair_id), primary_key=True)
    model_errors = Column(Text)
    model_filename = Column(String)
    chain_1 = Column(String)
    chain_2 = Column(String)
    norm_dope = Column(Float)
    interface_area_hydrophobic = Column(Float)
    interface_area_hydrophilic = Column(Float)
    interface_area_total = Column(Float)
    interface_dg = Column(Float)
    interacting_aa_1 = Column(Text)
    interacting_aa_2 = Column(Text)
    m_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    domain_pair = relationship(
        UniprotDomainPair, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('model', uselist=False, cascade='expunge')) # one to one


class UniprotDomainPairMutation(Base):
    __tablename__ = 'uniprot_domain_pair_mutation'
    __table_args__ = ({'schema': 'elaspic_v2'},)

    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False, primary_key=True)
    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPairModel.uniprot_domain_pair_id), primary_key=True)
    mutation = Column(String, nullable=False, primary_key=True)
    mutation_errors = Column(Text)
    model_filename_wt = Column(String)
    model_filename_mut = Column(String)
    chain_modeller = Column(String)
    mutation_modeller = Column(String)
    analyse_complex_energy_wt = Column(Text)
    stability_energy_wt = Column(Text)
    analyse_complex_energy_mut = Column(Text)
    stability_energy_mut = Column(Text)
    physchem_wt = Column(String)
    physchem_wt_ownchain = Column(String)
    physchem_mut = Column(String)
    physchem_mut_ownchain = Column(String)
    matrix_score = Column(Float)
    secondary_structure_wt = Column(String)
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(String)
    solvent_accessibility_mut = Column(Float)
    contact_distance_wt = Column(Float)
    contact_distance_mut = Column(Float)
    provean_score = Column(Float)
    ddg = Column(Float)
    mut_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    model = relationship(
        UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('mutations', cascade='expunge')) # many to one


###############################################################################
class MyDatabase(object):
    """
    """
    def __init__(
            self, sql_flavor=sql_flavor, is_immutable=False,
            path_to_temp='/tmp/', path_to_archive='/home/kimlab1/database_data/elaspic/',
            path_to_sqlite_db='', clear_database=False, logger=None):

        # Choose which database to use
        if sql_flavor == 'sqlite':
            autocommit=True
            autoflush=True
            engine = create_engine('sqlite://')
        elif sql_flavor == 'sqlite_file':
            autocommit=True
            autoflush=True
            engine = create_engine('sqlite:///' + path_to_sqlite_db, isolation_level='READ UNCOMMITTED')
        elif sql_flavor == 'postgresql':
            autocommit=False
            autoflush=False
            engine = create_engine('postgresql://elaspic:elaspic@192.168.6.19:5432/kimlab') # , echo=True

        # Commented out because this is dangerous. Clear all tables in the database.
        # if clear_database:
        #     Base.metadata.drop_all(engine)
        # Base.metadata.create_all(engine)

        if logger is None:
            logger = logging.getLogger(__name__)
            logger.handlers = []
            logger.setLevel(logging.DEBUG)
            handler = logging.StreamHandler()
            handler.setLevel(logging.DEBUG)
            logger.addHandler(handler)
        self.logger = logger

        Session.configure(bind=engine, autocommit=autocommit, autoflush=autoflush)
        self.Session = Session
        self.autocommit = autocommit
        self.is_immutable = is_immutable
        self.path_to_temp = path_to_temp
        self.path_to_archive = path_to_archive


    @contextmanager
    def session_scope(self):
        """ Provide a transactional scope around a series of operations.
        So you can use: `with self.session_scope() as session:`
        """
        session = self.Session()
        try:
            yield session
            session.commit()
        except:
            session.rollback()
            raise
        finally:
            session.close()


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


    def merge_row(self, row_instance):
        """ Add a list of rows *row_instances* to the database
        """
        with self.session_scope() as session:
            session.merge(row_instance)


    def get_domain(self, pfam_names, subdomains=False):
        """ Contains pdbfam-based definitions of all pfam domains in the pdb
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
        """ Keeps the domain-domain interaction information from pdbfam
        Note that the produced dataframe may not have the same order as the keys
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
            uniprot_domains = session\
                .query(UniprotDomain)\
                .filter(UniprotDomain.uniprot_id == uniprot_id)\
                .options(joinedload(UniprotDomain.model).joinedload('mutations'))\
                .all()

        if copy_data:
            for uniprot_domain in uniprot_domains:
                self._copy_uniprot_domain_data(uniprot_domain)

        return uniprot_domains


    def get_uniprot_domain_pair(self, uniprot_id, copy_data=False):
        """
        """
        with self.session_scope() as session:
            uniprot_domain_pairs = (
                session.query(UniprotDomainPair)
                .filter(
                    (UniprotDomainPair.uniprot_id_1==uniprot_id) |
                    (UniprotDomainPair.uniprot_id_2==uniprot_id) )
                .options(joinedload(UniprotDomainPair.model).joinedload('mutations'))
                .all() )

        if copy_data:
            for uniprot_domain_pair in uniprot_domain_pairs:
                self._copy_uniprot_domain_pair_data(uniprot_domain_pair)

        return uniprot_domain_pairs


    def _copy_uniprot_domain_data(self, uniprot_domain):
        if (uniprot_domain.path_to_data
                and uniprot_domain.model
                and uniprot_domain.model.alignment_filename
                and uniprot_domain.model.model_filename):
            tmp_save_path = self.path_to_temp + uniprot_domain.path_to_data
            archive_save_path = self.path_to_archive + uniprot_domain.path_to_data
            path_to_alignment = tmp_save_path + '/'.join(uniprot_domain.model.model_filename.split('/')[:-1]) + '/'
            subprocess.check_call('mkdir -p {}'.format(path_to_alignment), shell=True)
            subprocess.check_call('cp -f {} {}'.format(
                archive_save_path + uniprot_domain.model.model_filename,
                tmp_save_path + uniprot_domain.model.model_filename), shell=True)
            subprocess.check_call('cp -f {} {}'.format(
                archive_save_path + uniprot_domain.model.model_filename,
                tmp_save_path + uniprot_domain.model.model_filename), shell=True)
            if uniprot_domain.model.mutations:
                for mutation in uniprot_domain.model.mutations:
                    self._copy_mutation_data(mutation, uniprot_domain.path_to_data)


    def _copy_uniprot_domain_pair_data(self, uniprot_domain_pair):
        if (uniprot_domain_pair.path_to_data
                and uniprot_domain_pair.model
                and uniprot_domain_pair.model.alignment_filename_1
                and uniprot_domain_pair.model.alignment_filename_2
                and uniprot_domain_pair.model.model_filename):
            tmp_save_path = self.path_to_temp + uniprot_domain_pair.path_to_data
            archive_save_path = self.path_to_archive + uniprot_domain_pair.path_to_data
            path_to_alignment_1 = tmp_save_path + '/'.join(uniprot_domain_pair.model.alignment_filename_1.split('/')[:-1]) + '/'
            path_to_alignment_2 = tmp_save_path + '/'.join(uniprot_domain_pair.model.alignment_filename_2.split('/')[:-1]) + '/'
            subprocess.check_call('mkdir -p {}'.format(path_to_alignment_1), shell=True)
            subprocess.check_call('mkdir -p {}'.format(path_to_alignment_2), shell=True)
            subprocess.check_call('cp -f {} {}'.format(
                archive_save_path + uniprot_domain_pair.model.alignment_filename_1,
                tmp_save_path + uniprot_domain_pair.model.alignment_filename_1), shell=True)
            subprocess.check_call('cp -f {} {}'.format(
                archive_save_path + uniprot_domain_pair.model.alignment_filename_2,
                tmp_save_path + uniprot_domain_pair.model.alignment_filename_2), shell=True)
            subprocess.check_call('cp -f {} {}'.format(
                archive_save_path + uniprot_domain_pair.model.model_filename,
                tmp_save_path + uniprot_domain_pair.model.model_filename), shell=True)
            if uniprot_domain_pair.model.mutations:
                for mutation in uniprot_domain_pair.model.mutations:
                    self._copy_mutation_data(mutation, uniprot_domain_pair.path_to_data)


    def _copy_mutation_data(self, mutation, path_to_data):
        tmp_save_path = self.path_to_temp + path_to_data
        archive_save_path = self.path_to_archive + path_to_data
        path_to_mutation = tmp_save_path + '/'.join(mutation.model_filename_wt.split('/')[:-1]) + '/'
        subprocess.check_call('mkdir -p {}'.format(path_to_mutation), shell=True)
        subprocess.check_call('cp -f {} {}'.format(
            archive_save_path + mutation.model_filename_wt,
            tmp_save_path + mutation.model_filename_wt), shell=True)
        subprocess.check_call('cp -f {} {}'.format(
            archive_save_path + mutation.model_filename_mut,
            tmp_save_path + mutation.model_filename_mut), shell=True)





    def get_uniprot_sequence(self, uniprot_id, check_external=False):
        """ Return a Biopython SeqRecord object containg the sequence for the
        specified uniprot.
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
                print 'Fetching sequence for uniprot {} from an online server'.format(uniprot_id)
                address = 'http://www.uniprot.org/uniprot/{}.fasta'.format(uniprot_id)
                try:
                    handle = urllib2.urlopen(address)
                    sequence = next(SeqIO.parse(handle, "fasta"))
                except (StopIteration, urllib2.HTTPError) as e:
                    self.logger.debug('{}: {}'.format(type(e), str(e)))
                    print '{}: {}'.format(type(e), str(e))
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

        uniprot_seqrecord = SeqIO.SeqRecord(
            seq=Seq.Seq(str(uniprot_sequence.uniprot_sequence)),
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






    def get_uniprot_mutation(self, model, uniprot_id, mutation, path_to_data=False):
        """
        """
        with self.session_scope() as session:
            if type(model) == UniprotDomainModel:
                uniprot_mutation = session\
                    .query(UniprotDomainMutation)\
                    .filter(UniprotDomainMutation.uniprot_domain_id==model.uniprot_domain_id)\
                    .filter(UniprotDomainMutation.uniprot_id==uniprot_id)\
                    .filter(UniprotDomainMutation.mutation==mutation)\
                    .all()
            elif type(model) == UniprotDomainPairModel:
                uniprot_mutation = session\
                    .query(UniprotDomainPairMutation)\
                    .filter(UniprotDomainPairMutation.uniprot_domain_pair_id==model.uniprot_domain_pair_id)\
                    .filter(UniprotDomainPairMutation.uniprot_id==uniprot_id)\
                    .filter(UniprotDomainPairMutation.mutation==mutation)\
                    .all()
            else:
                raise Exception('Wrong format for model!')

#        if path_to_data:
#            tmp_save_path = self.path_to_temp + path_to_data
#            archive_save_path = self.path_to_archive + path_to_data
        if path_to_data:
            for mut in uniprot_mutation:
                if mut.model_filename_wt:
                    tmp_save_path = self.path_to_temp + path_to_data
                    archive_save_path = self.path_to_archive + path_to_data
                    mutation_save_subpath = mut.model_filename_wt.split('/')[0] + '/'
                    if (self.path_to_temp != self.path_to_archive):
                        # Not running on SciNet and have structures to save
                        subprocess.check_call(
                            'mkdir -p ' +
                            tmp_save_path + mutation_save_subpath, shell=True)
                        subprocess.check_call(
                            'cp -f ' +
                            archive_save_path + mut.model_filename_wt + ' ' +
                            tmp_save_path + mut.model_filename_wt, shell=True)
                        subprocess.check_call(
                            'cp -f ' +
                            archive_save_path + mut.model_filename_mut + ' ' +
                            tmp_save_path + mut.model_filename_mut, shell=True)

        return uniprot_mutation


    ###########################################################################
    def add_domain(self, d):
        with self.session_scope() as session:
            if isinstance(d, Domain):
                session\
                    .query(Domain)\
                    .filter(Domain.cath_id == d.cath_id)\
                    .update({Domain.domain_errors: d.domain_errors})
            elif isinstance(d, DomainContact):
                session\
                    .query(DomainContact)\
                    .filter(DomainContact.domain_contact_id == d.domain_contact_id)\
                    .update({DomainContact.domain_contact_errors: d.domain_contact_errors})


    def add_domain_errors(self, t, error_string):
        with self.session_scope() as session:
            if isinstance(t, UniprotDomainTemplate):
                domain = session\
                        .query(Domain)\
                        .filter(Domain.cath_id==t.cath_id)\
                        .as_scalar()
                domain.domain_errors = error_string
                session.merge(domain)
            elif isinstance(t, UniprotDomainPairTemplate):
                domain_contact = session\
                    .query(DomainContact)\
                    .filter(DomainContact.cath_id_1==t.cath_id_1)\
                    .filter(DomainContact.cath_id_2==t.cath_id_2)\
                    .all()[0]
                domain_contact.domain_contact_errors = error_string
                session.merge(domain_contact)
            else:
                raise Exception('Wrong type for template!!!')


    def merge_domain(self, d, path_to_data=False):
        """
        """
        # Save a copy of the alignment to the export folder
        if path_to_data:
            tmp_save_path = self.path_to_temp + path_to_data
            archive_save_path = self.path_to_archive + path_to_data
            subprocess.check_call('mkdir -p ' + archive_save_path, shell=True)
            # Save the row d as an sqlalchemy serial object
            d_serial = loads(d)
            pickle.dump(d_serial, open(archive_save_path + 'alignments.json', 'w'))
            # Save the alignments
            if self.path_to_temp != self.path_to_archive: # Not running on SciNet
                if isinstance(d, UniprotDomain):
                    subprocess.check_call('cp -f ' + tmp_save_path + d.alignment_filename +
                                            ' ' + archive_save_path + d.alignment_filename, shell=True)
#                    if d.provean_supset_filename:
#                        subprocess.check_call('cp -f ' + tmp_save_path + t.provean_supset_filename +
#                                            ' ' + archive_save_path + t.provean_supset_filename, shell=True)
                if isinstance(d, UniprotDomainPair):
                    subprocess.check_call('cp -f ' + tmp_save_path + d.alignment_filename_1 +
                                            ' ' + archive_save_path + d.alignment_filename_1, shell=True)
                    subprocess.check_call('cp -f ' + tmp_save_path + d.alignment_filename_2 +
                                            ' ' + archive_save_path + d.alignment_filename_2, shell=True)
        if not self.is_immutable:
            with self.session_scope() as session:
                session.merge(d)


    def add_model(self, uniprot_model, path_to_data=False):
        """
        """
        uniprot_model.m_date_modified = datetime.datetime.utcnow()






        # Save a copy of the alignment to the export folder
        if path_to_data:
            tmp_save_path = self.path_to_temp + path_to_data
            archive_save_path = self.path_to_archive + path_to_data

            with open(archive_save_path + 'model.json', 'w') as fh:
                json.dump(row2dict(uniprot_model), fh, indent=4, separators=(',', ': '))

            if (self.path_to_temp != self.path_to_archive) and (uniprot_model.model_filename is not None):
                # Not running on SciNet and have a structure to save
                subprocess.check_call('mkdir -p ' + archive_save_path, shell=True)
                subprocess.check_call('cp -f ' + tmp_save_path + uniprot_model.model_filename +
                                        ' ' + archive_save_path + uniprot_model.model_filename, shell=True)

        if not self.is_immutable:
            with self.session_scope() as session:
                session.merge(uniprot_model)






    def add_uniprot_mutation(self, uniprot_mutation, path_to_data=False):
        """
        """
        uniprot_mutation.mut_date_modified = datetime.datetime.utcnow()

        if path_to_data and (uniprot_mutation.model_filename_wt is not None):
            tmp_save_path = self.path_to_temp + path_to_data
            archive_save_path = self.path_to_archive + path_to_data
            archive_save_subpath = uniprot_mutation.model_filename_wt.split('/')[0] + '/'

            if not os.path.isdir(archive_save_path + archive_save_subpath):
                os.mkdir(archive_save_path + archive_save_subpath)
            with open(archive_save_path + archive_save_subpath + 'mutation.json', 'w') as fh:
                json.dump(row2dict(uniprot_mutation), fh, indent=4, separators=(',', ': '))

            if (self.path_to_temp != self.path_to_archive):
                # Not running on SciNet and have structures to save
                subprocess.check_call('mkdir -p ' + archive_save_path + archive_save_subpath, shell=True)
                subprocess.check_call('cp -f ' + tmp_save_path + uniprot_mutation.model_filename_wt +
                                        ' ' + archive_save_path + uniprot_mutation.model_filename_wt, shell=True)
                subprocess.check_call('cp -f ' + tmp_save_path + uniprot_mutation.model_filename_mut +
                                        ' ' + archive_save_path + uniprot_mutation.model_filename_mut, shell=True)

        if not self.is_immutable:
            with self.session_scope() as session:
                session.merge(uniprot_mutation)


    ###########################################################################
    def _split_domain(self, domain):
        """
        Takes a string of two domain boundaries and returns a list with int
        The separator is '-' and it can happen that both or one boundary is
        negative, i.e.

            -150-200,   meaning from -150 to 200
            -150--100,  meaning from -150 to -100, etc.

        NOTE! Currently the icode (see Biopython) is disregarded. That means
        that if the numbering is 3B, all '3's are taken. That is the letters
        are stripped! One might want to improve that behaviour.
        """
        # split the domain boundaries, keep eventual minus signs
        if domain[0] == '-' and len(domain[1:].split('-')) == 2:
            domain = ['-' + domain[1:].split('-')[0], domain[1:].split('-')[1]]
        elif domain[0] == '-' and len(domain[1:].split('-')) > 2:
            domain = ['-' + domain[1:].split('-')[0], '-' + domain[1:].split('-')[-1]]
        else:
            domain = [domain.split('-')[0], domain.split('-')[1]]
        # strip the letters
        if domain[0][-1] in uppercase:
            domain[0] = domain[0][:-1]
        if domain[1][-1] in uppercase:
            domain[1] = domain[1][:-1]
        domain = [int(domain[0]), int(domain[1])]
        return domain


    def _split_domain_semicolon(self, domains):
        """ Unlike split_domain(), this function returns a tuple of tuples of strings,
        preserving letter numbering (e.g. 10B)
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


#    def close(self):
#        if not self.autocommit:
#            self.session.commit()
#        self.session.close()


    ###########################################################################
    def get_alignment(self, uniprot_template, path_to_data):
        """
        """

        tmp_save_path = self.path_to_temp + path_to_data
        archive_save_path = self.path_to_archive + path_to_data

        if isinstance(uniprot_template, UniprotDomainTemplate):

            # Load previously-calculated alignments
            if os.path.isfile(tmp_save_path + uniprot_template.alignment_filename):
                alignment = AlignIO.read(tmp_save_path + uniprot_template.alignment_filename, 'clustal')
            elif os.path.isfile(archive_save_path + uniprot_template.alignment_filename):
                alignment = AlignIO.read(archive_save_path + uniprot_template.alignment_filename, 'clustal')
            else:
                raise error.NoPrecalculatedAlignmentFound(archive_save_path, uniprot_template.alignment_filename)

            return [alignment, None]

        elif isinstance(uniprot_template, UniprotDomainPairTemplate):

            # Read alignment from the temporary folder
            if (os.path.isfile(tmp_save_path + uniprot_template.alignment_filename_1)
            and os.path.isfile(tmp_save_path + uniprot_template.alignment_filename_2)):
                alignment_1 = AlignIO.read(tmp_save_path + uniprot_template.alignment_filename_1, 'clustal')
                alignment_2 = AlignIO.read(tmp_save_path + uniprot_template.alignment_filename_2, 'clustal')
            # Read alignment from the export database
            elif (os.path.isfile(archive_save_path + uniprot_template.alignment_filename_1)
            and os.path.isfile(archive_save_path + uniprot_template.alignment_filename_2)):
                alignment_1 = AlignIO.read(archive_save_path + uniprot_template.alignment_filename_1, 'clustal')
                alignment_2 = AlignIO.read(archive_save_path + uniprot_template.alignment_filename_2, 'clustal')
            else:
                raise error.NoPrecalculatedAlignmentFound(archive_save_path, uniprot_template.alignment_filename_1)

            return [alignment_1, alignment_2]


    ###########################################################################
    def load_db_from_csv(self):
        """
        """
        def pd_strip(text):
            # Strip tailing whitespace
            try:
                return text.strip()
            except AttributeError:
                return text

        ## Pupulate sql database from text files
        # Files from while the data will be loaded
        path_to_db = '/home/kimlab1/strokach/working/pipeline/db/'

        domain_infile = path_to_db + 'domain.txt'
        domain_contact_infile = path_to_db + 'domain_contact.txt'
        uniprot_sequence_infile = path_to_db + 'uniprot_sprot_human.txt'
        uniprot_domain_infile = path_to_db + 'uniprot_domain.txt'
        uniprot_domain_pair_infile = path_to_db + 'uniprot_domain_pair.txt'


        # Table `domain`
        names = ['cath_id', 'pdb_id', 'pdb_type', 'pdb_resolution', 'pdb_chain', 'pdb_domain_def', 'pfam_autopfam', 'pfam_name']
        domain_df = pd.read_csv(domain_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None, )
        domain_df['pdb_resolution'] = domain_df['pdb_resolution'].apply(lambda x: float(x))
        domain_df.drop_duplicates(['cath_id',])
        for idx, row in domain_df.iterrows():
            self.session.add(Domain(**row.to_dict()))
            if idx % 10000 == 0:
#                self.session.flush()
                self.logger.debug(idx)
        self.session.commit()
        self.logger.debug('Finished populating table domain')


        # Table `domain_contact`
        names = ['domain_contact_id', 'cath_id_1', 'contact_residues_1', 'cath_id_2', 'contact_residues_2']
        domain_contact_df = pd.read_csv(domain_contact_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
        domain_contact_df = domain_contact_df.dropna() # only a couple of rows are droppeds
        for idx, row in domain_contact_df.iterrows():
            self.session.add(DomainContact(**row.to_dict()))
            if idx % 10000 == 0:
#                self.session.flush()
                self.logger.debug(idx)
        self.session.commit()
        self.logger.debug('Finished populating table domain_contact')


        # Table `uniprot_sequence`
        names = ['uniprot_id', 'uniprot_name', 'uniprot_description', 'uniprot_sequence']
        uniprot_sequence_df = pd.read_csv(uniprot_sequence_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
        for idx, row in uniprot_sequence_df.iterrows():
            self.session.add(UniprotSequence(**row.to_dict()))
            if idx % 10000 == 0:
#                self.session.flush()
                self.logger.debug(idx)
        self.session.commit()
        self.logger.debug('Finished populating table uniprot_sequence')


        # Table `uniprot_domain`
        if os.path.isfile(uniprot_domain_infile):
#            names = ['uniprot_domain_id', 'uniprot_id', 'pfam_name', 'alignment_def']
            uniprot_domain_df_with_id = pd.read_csv(uniprot_domain_infile, sep='\t', na_values='\N', index_col=False)
            uniprot_domain_df_with_id['alignment_defs'] = uniprot_domain_df_with_id['alignment_def']
            uniprot_domain_df_with_id['alignment_def'] = uniprot_domain_df_with_id['alignment_defs'].apply(lambda x: encode_domain(decode_domain(x)))

            tmp = uniprot_domain_df_with_id.merge(uniprot_sequence_df, how='left', left_on='uniprot_id', right_on='uniprot_id', suffixes=('_domain', ''))
            uniprot_domain_df_with_id['organism_name'] = tmp['uniprot_name'].apply(lambda x: x.split('_')[-1])

            uniprot_domain_df_with_id['path_to_data'] = (
                uniprot_domain_df_with_id['organism_name'].apply(lambda x: x.lower()) + '/' +
                uniprot_domain_df_with_id['uniprot_id'].apply(lambda x: x[0:3]) + '/' +
                uniprot_domain_df_with_id['uniprot_id'].apply(lambda x: x[3:5]) + '/' +
                uniprot_domain_df_with_id['uniprot_id'] + '/' +
                uniprot_domain_df_with_id['pfam_name'] + '*' +
                uniprot_domain_df_with_id['alignment_def'].apply(lambda x: x.replace(':','-')) + '/')

            for idx, row in uniprot_domain_df_with_id.iterrows():
                self.session.add(UniprotDomain(**row.to_dict()))
                if idx % 10000 == 0:
#                    self.session.flush()
                    self.logger.debug(idx)
            self.session.commit()
            self.logger.debug('Finished populating table uniprot_domain')
        else:
            pass
#            pfam_parser = parse_pfamscan.make_uniprot_domain_database()
#            pfam_parser.run()
#            uniprot_domain_df = pfam_parser.get_dataframe()
#            uniprot_domain_df_with_id.to_csv(uniprot_domain_infile, sep='\t', na_rep='\N', index=False)


        # Table `uniprot_domain_pair`
        if os.path.isfile(uniprot_domain_pair_infile):
            uniprot_domain_pair_df_with_id = pd.read_csv(uniprot_domain_pair_infile, sep='\t', na_values='\N', index_col=False)

            temp = uniprot_domain_pair_df_with_id\
                .merge(uniprot_domain_df_with_id, how='left', left_on='uniprot_domain_id_1', right_on='uniprot_domain_id')\
                .merge(uniprot_domain_df_with_id, how='left', left_on='uniprot_domain_id_2', right_on='uniprot_domain_id', suffixes=('_1', '_2'))
            temp['path_to_data'] = (temp['path_to_data_1'] + temp['pfam_name_2'] + '*' + temp['alignment_def_2'].apply(lambda x: x.replace(':','-')) + '/' + temp['uniprot_id_2'] + '/')
            temp = temp[['uniprot_domain_pair_id', 'path_to_data']]
            uniprot_domain_pair_df_with_id = uniprot_domain_pair_df_with_id.merge(temp, how='left')

#            uniprot_domain_pair_df_with_id.to_sql('uniprot_domain_pair', conn, flavor=sql_flavor, if_exists='append')
            for idx, row in uniprot_domain_pair_df_with_id.iterrows():
                self.session.add(UniprotDomainPair(**row.to_dict()))
                if idx % 10000 == 0:
#                    self.session.flush()
                    self.logger.debug(idx)
            self.session.commit()
            self.logger.debug('Finished populating table domain')
        else:
            pass
#            pfam_parser = parse_pfamscan.make_uniprot_domain_pair_database(domain_df, domain_contact_df, uniprot_domain_df_with_id,
#                                infile='/home/kimlab1/strokach/working/databases/biogrid/pairs_of_interacting_uniprots_human.tsv')
#            uniprot_domain_pair_df = pfam_parser.get_dataframe()
##            uniprot_domain_pair_df.to_sql('uniprot_domain_pair', conn, flavor=sql_flavor, if_exists='append')
##            uniprot_domain_pair_df_with_id = pd.read_sql('SELECT * from uniprot_domain_pair', conn)
#            uniprot_domain_pair_df_with_id.to_csv(uniprot_domain_pair_infile, sep='\t', na_values='\N', index=False)


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
            childProcess = subprocess.Popen('ls ' + self.path_to_archive + d[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, __ = childProcess.communicate()
            filenames = [fname for fname in result.split('\n') if fname != '']
            for filename in filenames:
                with open(filename, 'r') as fh:
                    row = json.load(fh)
                try:
                    self.session.merge(d[1](**row))
                except TypeError as e:
                    self.logger.debug('Error merging %s.\nProbably from an older version of the database. Skipping...' % filename)
                    self.logger.debug('\t', e)
                self.logger.debug('Merged %s' % filename)
            self.session.commit()
            self.logger.debug('Committed changes\n\n\n')


###############################################################################
if __name__ == '__main__':
#    return
    # run to generate an initial state database (with no precalculatios)
    raise Exception
    print sql_flavor
    db = MyDatabase('/home/kimlab1/strokach/working/pipeline/db/pipeline.db',
                    path_to_archive='/home/kimlab1/database_data/elaspic/',
                    sql_flavor=sql_flavor,
                    clear_database=False)
#    db.load_db_from_csv()
    db.load_db_from_archive()
    db.session.close()

