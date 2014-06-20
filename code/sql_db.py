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
import tarfile
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy import Column, Index, UniqueConstraint
from sqlalchemy import Integer, Float, String, Boolean, Text, DateTime
from sqlalchemy import ForeignKey
from sqlalchemy.orm import sessionmaker, relationship, backref, aliased, scoped_session, joinedload
from sqlalchemy.ext.declarative import declarative_base

import Bio

#import parse_pfamscan
import errors as error
import helper_functions as hf

"""
Naming convention for classes includes an underscore since peewee converts
database table names to lowercase.
"""
###############################################################################
# Helper functions for dealing with sql objects

def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:bz2") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def decode_domain(domains, merge=True, return_string=False):
    """ Unlike split_domain(), this function returns a tuple of tuples of strings,
    preserving letter numbering (e.g. 10B)
    """
    if not domains:
        return None

    if domains[-1] == ',':
        domains = domains[:-1]
    x = domains
    if return_string:
        domain_fragments = [ [r.strip() for r in ro.split(':')] for ro in x.split(',') ]
    else:
        domain_fragments = [ [int(r.strip()) for r in ro.split(':')] for ro in x.split(',') ]
    domain_merged = [ domain_fragments[0][0], domain_fragments[-1][-1] ]
    if merge:
        return domain_merged
    else:
        return domain_fragments


def encode_domain(domains, merged=True):
    x = domains
    if merged:
        return ':'.join([str(r) for r in x])
    else:
        return ','.join([':'.join([str(r) for r in ro]) for ro in x])


def decode_aa_list(interface_aa):
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


def row2dict(row):
    d = {}
    for column in row.__table__.columns:
        d[column.name] = getattr(row, column.name)
        if type(d[column.name]) == datetime.datetime:
            d[column.name] = d[column.name].strftime('%Y-%m-%d %H-%M-%S-%f')

    return d


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
Session = sessionmaker(expire_on_commit=False)
#Session = scoped_session(sessionmaker(expire_on_commit=False))


#############################################################################

class Domain(Base):
    __tablename__ = 'domain'
    __table_args__ = ({'schema': 'elaspic'},)
    cath_id = Column(String(15, collation=binary_collation), primary_key=True)
    pdb_id = Column(String(4, collation=string_collation), nullable=False)
    pdb_type = Column(String(31, collation=string_collation), nullable=True)
    pdb_resolution = Column(Float, nullable=True)
    pdb_chain = Column(String(1, collation=string_collation), nullable=False)
    pdb_domain_def = Column(String(255, collation=string_collation), nullable=False)
    cdhit_cluster = Column(Integer, nullable=True)
    cdhit_cluster_idx = Column(Float, nullable=True)
    cdhit_cluster_length = Column(Integer, nullable=True)
    cdhit_cluster_identity = Column(Float, nullable=True)
    pfam_autopfam = Column(Integer, nullable=False)
    pfam_name = Column(String(255, collation=string_collation), nullable=False)
    domain_errors = Column(String(255, collation=string_collation))


class DomainContact(Base):
    __tablename__ = 'domain_contact'
    __table_args__ = (Index('cath_id_1_2', 'cath_id_1', 'cath_id_1', unique=True),
                      {'sqlite_autoincrement': True,
                      'schema': 'elaspic'},)
    domain_contact_id = Column(Integer, primary_key=True)
    cath_id_1 = Column(None, ForeignKey('elaspic.domain.cath_id'), index=True, nullable=False)
    cath_id_2 = Column(None, ForeignKey('elaspic.domain.cath_id'), index=True, nullable=False)
    min_interchain_distance = Column(Float)
    contact_volume = Column(Float)
    contact_surface_area = Column(Float)
    atom_count_1 = Column(Integer)
    atom_count_2 = Column(Integer)
    contact_residues_1 = Column(Text)
    contact_residues_2 = Column(Text)
    # Relationships
    domain_1 = relationship(Domain, primaryjoin=cath_id_1==Domain.cath_id, cascade='expunge', lazy='joined')
    domain_2 = relationship(Domain, primaryjoin=cath_id_2==Domain.cath_id, cascade='expunge', lazy='joined')
    domain_contact_errors = Column(String(255, collation=string_collation))

###############################################################################

class UniprotSequence(Base):
    __tablename__ = 'uniprot_sequence'
    __table_args__ = ({'schema': 'elaspic'},)

    db = Column(String(12, collation=string_collation), nullable=False)
    uniprot_id = Column(String(50, collation=string_collation), primary_key=True, nullable=False)
    uniprot_name = Column(String(255, collation=string_collation), nullable=False)
    protein_name = Column(String(255, collation=string_collation))
    organism_name = Column(String(255, collation=string_collation))
    gene_name = Column(String(255, collation=string_collation))
    protein_existence = Column(Integer)
    sequence_version = Column(Integer)
    uniprot_sequence = Column(Text, nullable=False)


class UniprotDomain(Base):
    __tablename__ = 'uniprot_domain'
    __table_args__ = (UniqueConstraint('uniprot_id', 'pfam_name', 'envelope_def', name='uniprot_domain_tuple'),
                      {'schema': 'elaspic'},)

    uniprot_domain_id = Column(Integer, primary_key=True)
    db = Column(String)
    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False)
    uniprot_name = Column(String)
    organism_name = Column(String(255, collation=string_collation), nullable=False)
    pfam_id = Column(String)
    pfam_ids = Column(Text)
    pfam_name = Column(String, index=True, nullable=False) #seq_id
    pfam_names = Column(Text)
    clan_id = Column(String)
    clan_ids = Column(Text)
    clan_name = Column(String)
    clan_names = Column(Text)
    envelope_def = Column(String(255, collation=string_collation), nullable=False)
    envelope_defs = Column(Text)

    path_to_data = Column(Text)

    # Relationships
    uniprot_sequence = relationship(UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=backref('domain', cascade='expunge')) # many to one


class UniprotDomainPair(Base):
    __tablename__ = 'uniprot_domain_pair'
    __table_args__ = (Index('uniprot_domain_id_1_2', "uniprot_domain_id_1", "uniprot_domain_id_2", unique=True),
                      {'sqlite_autoincrement': True,
                      'schema': 'elaspic'},)

    uniprot_domain_pair_id = Column(Integer, primary_key=True)
    uniprot_domain_id_1 = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), index=True, nullable=False)
    uniprot_domain_id_2 = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), index=True, nullable=False)
    path_to_data = Column(Text)
    rigids = Column(Text)
    domain_contact_ids = Column(Text)

    # Relationships
    uniprot_domain_1 = relationship(UniprotDomain,
        primaryjoin=uniprot_domain_id_1==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one
    uniprot_domain_2 = relationship(UniprotDomain,
        primaryjoin=uniprot_domain_id_2==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one


###############################################################################

class UniprotDomainTemplate(Base):
    __tablename__ = 'uniprot_domain_template'
    __table_args__ = ({'schema': 'elaspic'},)

    uniprot_domain_id = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), primary_key=True)
    template_errors = Column(Text)
    alignment_filename = Column(String(255, collation=binary_collation))
    provean_supset_filename = Column(String(255, collation=string_collation))

    # failed_cath_id = Column(Text)
    cath_id = Column(None, ForeignKey(Domain.cath_id), index=True)
    domain_def = Column(String(255, collation=string_collation))
    alignment_id = Column(String(255, collation=string_collation))
    alignment_score = Column(Integer)
    alignment_identity = Column(Float)

    #
    t_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    uniprot_domain = relationship(UniprotDomain, uselist=False, cascade='expunge', lazy='joined',
                                  backref=backref('template', uselist=False, cascade='expunge')) # one to one
    domain = relationship(Domain, uselist=False, cascade='expunge', lazy='joined')


class UniprotDomainModel(Base):
    __tablename__ = 'uniprot_domain_model'
    __table_args__ = ({'schema': 'elaspic'},)

    uniprot_domain_id = Column(None, ForeignKey(UniprotDomainTemplate.uniprot_domain_id), primary_key=True)
    model_errors = Column(Text)
    model_filename = Column(String(255, collation=binary_collation))

    chain = Column(String(1, collation=string_collation))
    norm_dope = Column(Float)
    het_flag = Column(Boolean)
    switch_chain = Column(Boolean)

    sasa_score = Column(Text)

    #
    m_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    template = relationship(UniprotDomainTemplate, uselist=False, cascade='expunge', lazy='joined',
                            backref=backref('model', uselist=False, cascade='expunge')) # one to one


class UniprotDomainMutation(Base):
    __tablename__ = 'uniprot_domain_mutation'
    __table_args__ = ({'schema': 'elaspic'},)

    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False, primary_key=True)
    uniprot_domain_id = Column(None, ForeignKey(UniprotDomainModel.uniprot_domain_id), primary_key=True)
    mutation = Column(String(8, collation=string_collation), nullable=False, primary_key=True)
    mutation_errors = Column(Text)

    model_filename_wt = Column(String(255, collation=binary_collation))
    model_filename_mut = Column(String(255, collation=binary_collation))
    #
    chain_modeller = Column(String(1, collation=string_collation))
    mutation_modeller = Column(String(8, collation=string_collation))
    #
    Stability_energy_wt = Column(Text)
    Stability_energy_mut = Column(Text)
    physChem_wt = Column(String(255, collation=binary_collation))
    physChem_wt_ownChain = Column(String(255, collation=binary_collation))
    physChem_mut = Column(String(255, collation=binary_collation))
    physChem_mut_ownChain = Column(String(255, collation=binary_collation))
    matrix_score = Column(Float)
    secondary_structure_wt = Column(String(1, collation=string_collation))
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(String(1, collation=string_collation))
    solvent_accessibility_mut = Column(Float)
    provean_score = Column(Float)
    ddG = Column(Float)
    mut_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    model = relationship(UniprotDomainModel, cascade='expunge', uselist=False, lazy='joined',
                         backref=backref('mutations', cascade='expunge')) # many to one


###############################################################################

class UniprotDomainPairTemplate(Base):
    __tablename__ = 'uniprot_domain_pair_template'
    __table_args__ = ({'schema': 'elaspic'},)

    # Columns
    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPair.uniprot_domain_pair_id), primary_key=True)
    template_errors = Column(Text)
    alignment_filename_1 = Column(String(255, collation=binary_collation))
    alignment_filename_2 = Column(String(255, collation=binary_collation))

#    domain_contact_id = Column(None, ForeignKey(DomainContact.domain_contact_id), index=True)
    # failed_cath_ids = Column(Text) # cath_id_1)underscore)cath_id_2
    cath_id_1 = Column(None, ForeignKey(Domain.cath_id), index=True)
    domain_def_1 = Column(String(255, collation=string_collation))
    alignment_id_1 = Column(String(255, collation=string_collation))
    alignment_score_1 = Column(Integer)
    alignment_identity_1 = Column(Float)

    cath_id_2 = Column(None, ForeignKey(Domain.cath_id), index=True)
    domain_def_2 = Column(String(255, collation=string_collation))
    alignment_id_2 = Column(String(255, collation=string_collation))
    alignment_score_2 = Column(Integer)
    alignment_identity_2 = Column(Float)

    #
    t_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    uniprot_domain_pair = relationship(UniprotDomainPair, uselist=False, cascade='expunge', lazy='joined',
                                       backref=backref('template', uselist=False, cascade='expunge')) # one to one
    domain_1 = relationship(Domain, primaryjoin=cath_id_1==Domain.cath_id, cascade='expunge', lazy='joined')
    domain_2 = relationship(Domain, primaryjoin=cath_id_2==Domain.cath_id, cascade='expunge', lazy='joined')


class UniprotDomainPairModel(Base):
    __tablename__ = 'uniprot_domain_pair_model'
    __table_args__ = ({'schema': 'elaspic'},)

    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPairTemplate.uniprot_domain_pair_id), primary_key=True)
    model_errors = Column(Text)
    model_filename = Column(String(255, collation=binary_collation))
    #
    chain_1 = Column(String(1, collation=string_collation))
    chain_2 = Column(String(1, collation=string_collation))

    norm_dope = Column(Float)
    het_flag_1 = Column(Boolean)
    het_flag_2 = Column(Boolean)
    switch_chain = Column(Boolean)

    interface_area_hydrophobic = Column(Float)
    interface_area_hydrophilic = Column(Float)
    interface_area_total = Column(Float)
    interface_dG = Column(Float)
    interacting_aa_1 = Column(Text)
    interacting_aa_2 = Column(Text)

    #
    m_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    template = relationship(UniprotDomainPairTemplate, uselist=False, cascade='expunge', lazy='joined',
                            backref=backref('model', uselist=False, cascade='expunge')) # one to one


class UniprotDomainPairMutation(Base):
    __tablename__ = 'uniprot_domain_pair_mutation'
    __table_args__ = ({'schema': 'elaspic'},)

    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False, primary_key=True)
    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPairModel.uniprot_domain_pair_id), primary_key=True)
    mutation = Column(String(8, collation=string_collation), nullable=False, primary_key=True)
    mutation_errors = Column(Text)

    model_filename_wt = Column(String(255, collation=binary_collation))
    model_filename_mut = Column(String(255, collation=binary_collation))
    #
    chain_modeller = Column(String(1, collation=string_collation))
    mutation_modeller = Column(String(8, collation=string_collation))
    #
    AnalyseComplex_energy_wt = Column(Text)
    Stability_energy_wt = Column(Text)
    AnalyseComplex_energy_mut = Column(Text)
    Stability_energy_mut = Column(Text)

    physChem_wt = Column(String(255, collation=binary_collation))
    physChem_wt_ownChain = Column(String(255, collation=binary_collation))
    physChem_mut = Column(String(255, collation=binary_collation))
    physChem_mut_ownChain = Column(String(255, collation=binary_collation))

    matrix_score = Column(Float)

    secondary_structure_wt = Column(String(1, collation=string_collation))
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(String(1, collation=string_collation))
    solvent_accessibility_mut = Column(Float)
    contact_distance_wt = Column(Float)
    contact_distance_mut = Column(Float)

    provean_score = Column(Float)

    ddG = Column(Float)

    #
    mut_date_modified = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, nullable=False)

    # Relationships
    model = relationship(UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
                         backref=backref('mutations', cascade='expunge')) # many to one


###############################################################################

@contextmanager
def session_scope():
    """Provide a transactional scope around a series of operations.
    so you can use: 'with session_scope() as session:'
    """
    session = Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.close()


class MyDatabase(object):
    """
    """

    def __init__(self, path_to_sqlite_db='', sql_flavor=sql_flavor, is_immutable=False,
                 path_to_temp='/tmp/', path_to_archive='/home/kimlab1/database_data/elaspic/',
                 clear_database=False, log=None):
        """
        """
        # Based on configuration, use a different database
        # (Defaults to an in-memmory sqlite engine)
        if sql_flavor == 'sqlite':
            autocommit=True
            autoflush=True
            engine = create_engine('sqlite://')
        elif sql_flavor == 'sqlite_file':
            autocommit=True
            autoflush=True
            engine = create_engine('sqlite:///' + path_to_sqlite_db, isolation_level='READ UNCOMMITTED')
        elif sql_flavor == 'mysql':
            autocommit=True
            autoflush=True
            engine = create_engine('mysql://elaspic:elaspic@192.168.6.19:3306/elaspic')
        elif sql_flavor == 'mysql':
            autocommit=False
            autoflush=False
            engine = create_engine('mysql://root:kim630@127.0.0.1:3306/elaspic')
        elif sql_flavor == 'postgresql':
            autocommit=False
            autoflush=False
            engine = create_engine('postgresql://elaspic:elaspic@192.168.6.19:5432/kimlab')

#        if clear_database:
#            Base.metadata.drop_all(engine)
#        Base.metadata.create_all(engine)

        if log is None:
            logger = logging.getLogger(__name__)
            logger.setLevel(logging.DEBUG)
        #    handler = logging.FileHandler(tmp_path + 'templates.log', mode='w', delay=True)
            handler = logging.StreamHandler()
            handler.setLevel(logging.DEBUG)
            logger.addHandler(handler)
            log = logger
        self.log = log

        Session.configure(bind=engine, autocommit=autocommit, autoflush=autoflush)
        self.Session = Session
        self.autocommit = autocommit
#        self.session = Session()
        self.is_immutable = is_immutable
        self.path_to_temp = path_to_temp
        self.path_to_archive = path_to_archive


    ###########################################################################
    def get_uniprot_sequence(self, uniprot_id):
        """
        """
        if uniprot_id in ['A6NF79', 'C9JUS1', 'Q6N045', 'A6NMD1']:
            # these uniprotKBs made problems
            return []
        elif uniprot_id == 'P02735':
            # this sequence got replaced. I don't know right now if I can take
            # replaced sequence so I rather dismiss it.
            return []

        with session_scope() as session:
            uniprot_sequence = session.query(UniprotSequence)\
                                .filter(UniprotSequence.uniprot_id==uniprot_id)\
                                .all()

        if len(uniprot_sequence) == 1:
            uniprot_sequence = uniprot_sequence[0]

        elif len(uniprot_sequence) == 0:
            # the True/False value is used to add new sequences to the database in
            # end. Only usefull if you run one instance at a time otherwise you will
            # get diverging sequence databases.
            username = hf.get_username()
            if username.strip() == 'joan':
                self.log.debug('uniprot sequence not found')
            else:
                self.log.debug('Fetching uniprot sequence', uniprot_id, 'from online server')
                address = 'http://www.uniprot.org/uniprot/' + uniprot_id + '.fasta'
                handle = urllib2.urlopen(address)
                sequence = next(Bio.SeqIO.parse(handle, "fasta"))

                uniprot_sequence = UniprotSequence()
                uniprot_sequence.uniprot_id = uniprot_id
                uniprot_sequence.uniprot_name = sequence.name
                uniprot_sequence.uniprot_description = sequence.description
                uniprot_sequence.uniprot_sequence = str(sequence.seq)
#                self.session.add(uniprot_sequence)

        if uniprot_sequence == []:
            raise error.NoSequenceFound('No sequence found for ' + uniprot_id)

        uniprot_seqio_object = Bio.SeqIO.SeqRecord(seq=Bio.Seq.Seq(str(uniprot_sequence.uniprot_sequence)),
                               id=uniprot_sequence.uniprot_id,
                               name=uniprot_sequence.uniprot_name)

        return uniprot_seqio_object


    def add_uniprot_sequence(self, uniprot_sequence):
        """
        Add the new items (which is a list of tuples) to the database
        """
        with session_scope() as session:
            session.add(uniprot_sequence)

#        if not self.is_immutable:
#            self.session.add(uniprot_sequence)
#            if not self.autocommit:
#                self.session.commit()


    ###########################################################################
    def _get_rows_by_ids(self, row_object, row_object_identifiers, row_object_identifier_values):
        with session_scope() as session:
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


    def _update_rows(self, row_instances):
        with session_scope() as session:
            for row_instance in row_instances:
                session.merge(row_instance)


    def get_domain(self, pfam_names, subdomains=False):
        """
        Contains pdbfam-based definitions of all pfam domains in the pdb
        """
        with session_scope() as session:
            domain_set = set()
            for pfam_name in pfam_names:
                if subdomains:
                    domain = (
                        session.query(Domain).filter(
                            (Domain.pfam_name.like(pfam_name)) |
                            (Domain.pfam_name.like(pfam_name+'+%')) |
                            (Domain.pfam_name.like(pfam_name+'\_%')) | # need an escape character because _ matches any single character
                            (Domain.pfam_name.like('%+'+pfam_name)) |
                            (Domain.pfam_name.like('%+'+pfam_name+'+%')) |
                            (Domain.pfam_name.like('%+'+pfam_name+'\_%')) ) # need an escape character because _ matches any single character
                        .distinct().all() )
                else:
                    domain = (
                        session.query(Domain)
                        .filter(Domain.pfam_name==pfam_name)
                        .distinct().all() )
                domain_set.update(domain)
        if not domain_set:
            self.log.debug('No domain definitions found for pfam: %s' % str(pfam_names))
        return list(domain_set)


    def get_domain_contact(self, pfam_names_1, pfam_names_2, subdomains=False):
        """ Keeps the domain-domain interaction information from pdbfam
        Note that the produced dataframe may not have the same order as the keys
        """
        with session_scope() as session:
            domain_contact_1 = self._get_domain_contact(pfam_names_1, pfam_names_2, session, subdomains)
            domain_contact_2 = self._get_domain_contact(pfam_names_2, pfam_names_1, session, subdomains)

        if not len(domain_contact_1) and not len(domain_contact_2):
            self.log.debug('No domain contact template found for domains %s, %s' % (str(pfam_names_1), str(pfam_names_2),))

        return [domain_contact_1, domain_contact_2]


    def _get_domain_contact(self, pfam_names_1, pfam_names_2, session, subdomains):
        """
        """
        domain_1 = aliased(Domain)
        domain_2 = aliased(Domain)
        domain_contact_set = set()
        for pfam_name_1 in pfam_names_1:
            for pfam_name_2 in pfam_names_2:
                if subdomains:
                    domain_contact = (
                        session.query(DomainContact)
                        .join(domain_1, DomainContact.cath_id_1==domain_1.cath_id)
                        .filter(
                            (domain_1.pfam_name.like(pfam_name_1)) |
                            (domain_1.pfam_name.like(pfam_name_1+'+%')) |
                            (domain_1.pfam_name.like(pfam_name_1+'\_%')) | # need an escape character because _ matches any single character
                            (domain_1.pfam_name.like('%+'+pfam_name_1)) |
                            (domain_1.pfam_name.like('%+'+pfam_name_1+'+%')) |
                            (domain_1.pfam_name.like('%+'+pfam_name_1+'\_%')) ) # need an escape character because _ matches any single character
                        .join(domain_2, DomainContact.cath_id_2==domain_2.cath_id)
                        .filter(
                            (domain_2.pfam_name.like(pfam_name_2)) |
                            (domain_2.pfam_name.like(pfam_name_2+'+%')) |
                            (domain_2.pfam_name.like(pfam_name_2+'\_%')) | # need an escape character because _ matches any single character
                            (domain_2.pfam_name.like('%+'+pfam_name_2)) |
                            (domain_2.pfam_name.like('%+'+pfam_name_2+'+%')) |
                            (domain_2.pfam_name.like('%+'+pfam_name_2+'\_%')) ) # need an escape character because _ matches any single character
                        .distinct().all() )
                else:
                    domain_contact = (
                        session.query(DomainContact)
                        .join(domain_1, DomainContact.cath_id_1==domain_1.cath_id)
                        .filter(domain_1.pfam_name==pfam_name_1)
                        .join(domain_2, DomainContact.cath_id_2==domain_2.cath_id)
                        .filter(domain_2.pfam_name==pfam_name_2)
                        .distinct().all() )
                domain_contact_set.update(domain_contact)
        return list(domain_contact_set)


    ###########################################################################
    def get_uniprot_domain(self, uniprot_id, copy_data=True):
        """
        Initiated using parsed pfam_scan.pl output for the human uniprot (or the entire uniprot)
        The database is updated with information about calculated models
        """
        with session_scope() as session:
            uniprot_definitions = (
                session
                .query(UniprotDomain, UniprotDomainTemplate, UniprotDomainModel)
                .filter(UniprotDomain.uniprot_id==uniprot_id)
                .outerjoin(UniprotDomainTemplate)
                .outerjoin(UniprotDomainModel)
                .all() )

        if len(uniprot_definitions)==0:
            self.log.debug('No domain found in uniprot %s' % uniprot_id)

        if copy_data:
            for d, t, m in uniprot_definitions:
                if d.path_to_data:
                    tmp_save_path = self.path_to_temp + d.path_to_data
                    archive_save_path = self.path_to_archive + d.path_to_data
                    if t:
                        subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
                        if t.alignment_filename:
                            subprocess.check_call('cp -f ' + archive_save_path + t.alignment_filename +
                                                    ' ' + tmp_save_path + t.alignment_filename, shell=True)
                        if t.provean_supset_filename:
                            subprocess.check_call('cp -f ' + archive_save_path + t.provean_supset_filename +
                                                    ' ' + tmp_save_path + t.provean_supset_filename, shell=True)
                    if m:
                        if m.model_filename:
                            subprocess.check_call('cp -f ' + archive_save_path + m.model_filename +
                                                    ' ' + tmp_save_path + m.model_filename, shell=True)
        return uniprot_definitions


    def get_uniprot_domain_pair(self, uniprot_id, copy_data=True):
        """
        Contains known interactions between uniprot proteins
        Checks if the interaction database is already available as a pickled object.
        Generates it from file 'textfile_name' otherwise.
        """
        with session_scope() as session:
            uniprot_domain_pair_1 = self._get_uniprot_domain_pair(uniprot_id, session, reverse=False)
            uniprot_domain_pair_2 = self._get_uniprot_domain_pair(uniprot_id, session, reverse=True)

        if len(uniprot_domain_pair_1)==0 and len(uniprot_domain_pair_2)==0:
            self.log.debug('No known interactions with uniprot %s' % uniprot_id)

        if copy_data:
            for d, t, m in uniprot_domain_pair_1 + uniprot_domain_pair_2:
                if d.path_to_data:
                    tmp_save_path = self.path_to_temp + d.path_to_data
                    archive_save_path = self.path_to_archive + d.path_to_data
                    if t:
                        if t.alignment_filename_1 and t.alignment_filename_2:
                            subprocess.check_call(
                                'mkdir -p ' + tmp_save_path +
                                t.alignment_filename_1[:t.alignment_filename_1.rfind('/')], shell=True)
                            subprocess.check_call(
                                'cp -f ' + archive_save_path + t.alignment_filename_1 +
                                ' ' + tmp_save_path + t.alignment_filename_1, shell=True)
                            subprocess.check_call('cp -f ' + archive_save_path + t.alignment_filename_2 +
                                                    ' ' + tmp_save_path + t.alignment_filename_2, shell=True)
                    if m and (m.model_filename is not None):
                        subprocess.check_call('cp -f ' + archive_save_path + m.model_filename +
                                                ' ' + tmp_save_path + m.model_filename, shell=True)

#        return [uniprot_domain_pair_1, uniprot_domain_pair_2]
        return list(set(uniprot_domain_pair_1 + uniprot_domain_pair_2))


    def _get_uniprot_domain_pair(self, uniprot_id, session, reverse=False):
        """
        """
        if not reverse:
            uniprot_id_of_reference_domain = UniprotDomainPair.uniprot_domain_id_1
        else:
            uniprot_id_of_reference_domain = UniprotDomainPair.uniprot_domain_id_2

        uniprot_domain_pair = session\
            .query(UniprotDomainPair, UniprotDomainPairTemplate, UniprotDomainPairModel)\
            .join(UniprotDomain, UniprotDomain.uniprot_domain_id==uniprot_id_of_reference_domain)\
            .filter(UniprotDomain.uniprot_id == uniprot_id)\
            .outerjoin(UniprotDomainPairTemplate)\
            .outerjoin(UniprotDomainPairModel)\
            .options(joinedload(UniprotDomainPair.uniprot_domain_1).joinedload('template').joinedload('model'))\
            .options(joinedload(UniprotDomainPair.uniprot_domain_2).joinedload('template').joinedload('model'))\
            .all()

        return uniprot_domain_pair


    def get_uniprot_mutation(self, model, uniprot_id, mutation, path_to_data=False):
        """
        """
        with session_scope() as session:
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
        with session_scope() as session:
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
        with session_scope() as session:
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


    def add_uniprot_template(self, t, path_to_data=False):
        """
        """
        t.t_date_modified = datetime.datetime.utcnow()

        # Save a copy of the alignment to the export folder
        if path_to_data:
            tmp_save_path = self.path_to_temp + path_to_data
            archive_save_path = self.path_to_archive + path_to_data
            subprocess.check_call('mkdir -p ' + archive_save_path, shell=True)

            with open(archive_save_path + 'template.json', 'w') as fh:
                json.dump(row2dict(t), fh, indent=4, separators=(',', ': '))

            if self.path_to_temp != self.path_to_archive: # Not running on SciNet
                if isinstance(t, UniprotDomainTemplate):
                    if t.alignment_filename:
                        subprocess.check_call('cp -f ' + tmp_save_path + t.alignment_filename +
                                                ' ' + archive_save_path + t.alignment_filename, shell=True)
                    if t.provean_supset_filename:
                        subprocess.check_call('cp -f ' + tmp_save_path + t.provean_supset_filename +
                                            ' ' + archive_save_path + t.provean_supset_filename, shell=True)
                if isinstance(t, UniprotDomainPairTemplate):
                    if t.alignment_filename_1:
                        subprocess.check_call('cp -f ' + tmp_save_path + t.alignment_filename_1 +
                                                ' ' + archive_save_path + t.alignment_filename_1, shell=True)
                        subprocess.check_call('cp -f ' + tmp_save_path + t.alignment_filename_2 +
                                                ' ' + archive_save_path + t.alignment_filename_2, shell=True)
        if not self.is_immutable:
            with session_scope() as session:
                session.merge(t)


    def add_uniprot_model(self, uniprot_model, path_to_data=False):
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
            with session_scope() as session:
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
            with session_scope() as session:
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
                alignment = Bio.AlignIO.read(tmp_save_path + uniprot_template.alignment_filename, 'clustal')
            elif os.path.isfile(archive_save_path + uniprot_template.alignment_filename):
                alignment = Bio.AlignIO.read(archive_save_path + uniprot_template.alignment_filename, 'clustal')
            else:
                raise error.NoPrecalculatedAlignmentFound(archive_save_path, uniprot_template.alignment_filename)

            return [alignment, None]

        elif isinstance(uniprot_template, UniprotDomainPairTemplate):

            # Read alignment from the temporary folder
            if (os.path.isfile(tmp_save_path + uniprot_template.alignment_filename_1)
            and os.path.isfile(tmp_save_path + uniprot_template.alignment_filename_2)):
                alignment_1 = Bio.AlignIO.read(tmp_save_path + uniprot_template.alignment_filename_1, 'clustal')
                alignment_2 = Bio.AlignIO.read(tmp_save_path + uniprot_template.alignment_filename_2, 'clustal')
            # Read alignment from the export database
            elif (os.path.isfile(archive_save_path + uniprot_template.alignment_filename_1)
            and os.path.isfile(archive_save_path + uniprot_template.alignment_filename_2)):
                alignment_1 = Bio.AlignIO.read(archive_save_path + uniprot_template.alignment_filename_1, 'clustal')
                alignment_2 = Bio.AlignIO.read(archive_save_path + uniprot_template.alignment_filename_2, 'clustal')
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
                self.log.debug(idx)
        self.session.commit()
        self.log.debug('Finished populating table domain')


        # Table `domain_contact`
        names = ['domain_contact_id', 'cath_id_1', 'contact_residues_1', 'cath_id_2', 'contact_residues_2']
        domain_contact_df = pd.read_csv(domain_contact_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
        domain_contact_df = domain_contact_df.dropna() # only a couple of rows are droppeds
        for idx, row in domain_contact_df.iterrows():
            self.session.add(DomainContact(**row.to_dict()))
            if idx % 10000 == 0:
#                self.session.flush()
                self.log.debug(idx)
        self.session.commit()
        self.log.debug('Finished populating table domain_contact')


        # Table `uniprot_sequence`
        names = ['uniprot_id', 'uniprot_name', 'uniprot_description', 'uniprot_sequence']
        uniprot_sequence_df = pd.read_csv(uniprot_sequence_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
        for idx, row in uniprot_sequence_df.iterrows():
            self.session.add(UniprotSequence(**row.to_dict()))
            if idx % 10000 == 0:
#                self.session.flush()
                self.log.debug(idx)
        self.session.commit()
        self.log.debug('Finished populating table uniprot_sequence')


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
                    self.log.debug(idx)
            self.session.commit()
            self.log.debug('Finished populating table uniprot_domain')
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
                    self.log.debug(idx)
            self.session.commit()
            self.log.debug('Finished populating table domain')
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
                    self.log.debug('Error merging %s.\nProbably from an older version of the database. Skipping...' % filename)
                    self.log.debug('\t', e)
                self.log.debug('Merged %s' % filename)
            self.session.commit()
            self.log.debug('Committed changes\n\n\n')


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

