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
from contextlib import contextmanager
from collections import deque

import pandas as pd
import sqlalchemy as sa
import sqlalchemy.ext.declarative as sa_ext_declarative
import sqlalchemy.ext.serializer as sa_ext_serializer

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

DB_TYPE = conf.configs.get('db_type', None)
DB_DATABASE = conf.configs.get('db_database', 'elaspic')
DB_SCHEMA = conf.configs.get('db_schema', 'elaspic')
DB_SCHEMA_UNIPROT = conf.configs.get('db_schema_uniprot', 'elaspic')

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


if DB_TYPE is None:
    print('The `DB_TYPE` has not been set. Do not know what database is being used!')

def get_db_specific_param(key):
    if DB_TYPE is None:
        return
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
            column_names, kwargs = columns
        elif type(columns) == list:
            column_names = columns
            kwargs = {}
            kwargs['unique'] = False
        if 'index_name' in kwargs:
            index_name = kwargs.pop('index_name')
        else:
            index_name = (
                'ix_{table_name}_{column_0_name}'
                .format(table_name=table_name, column_0_name=column_names[0])[:255]
            )
        table_args.append(sa.Index(index_name, *column_names, **kwargs))
    # Other table parameters, such as schemas, etc.
    for db_specific_param in db_specific_params:
        table_args.append(get_db_specific_param(db_specific_param))
    return tuple(table_args)
   


#%%
Base = sa_ext_declarative.declarative_base()
Base.metadata.naming_conventions = naming_convention

#: Get the session that will be used for all future queries.
#: `expire_on_commit` so that you keep all the table objects even after the session closes.
Session = sa.orm.sessionmaker(expire_on_commit=False)
#Session = scoped_session(sa.orm.sessionmaker(expire_on_commit=False))


#%%
class Domain(Base):
    """
    Profs domain definitions for all proteins in the PDB. 

    Columns:
      cath_id
        Unique id identifying each domain in the PDB. Constructed by concatenating the pdb_id,
        pdb_chain, and an index specifying the order of the domain in the chain.

      pdb_id
        The PDB id in which the domain is found.

      pdb_chain
        The PDB chain in which the domain is found.

      pdb_domain_def
        Domain definitions of the domain, in PDB RESNUM coordinates.

      pdb_pdbfam_name
        The Profs name of the domain.

      pdb_pdbfam_idx
        An integer specifying the number of times a domain with domain name ``pdb_pdbfam_name`` has 
        occurred in this chain up to this point. It is used to make every 
        ``(pdb_id, pdb_chain, pdb_pdbfam_name, pdb_pdbfam_idx)`` tuple unique.
      
      domain_errors
        List of errors that occurred when annotating this domain, or when using this domain
        to make structural homology models.
    """
    __tablename__ = 'domain'
    if DB_TYPE != 'mysql':
        # MySQL can't handle long indexes
        _indexes = [
            (['pdb_id', 'pdb_chain', 'pdb_pdbfam_name', 'pdb_pdbfam_idx'], {'unique': True}),
            (['pdb_pdbfam_name'], {'mysql_length': 255})
        ]
    else:
        _indexes = [
            ['pdb_id', 'pdb_chain'],
            (['pdb_pdbfam_name'], {'mysql_length': 255}),
        ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
    cath_id = sa.Column(
        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')), 
        primary_key=True)
    pdb_id = sa.Column(sa.String(SHORT), nullable=False)
    pdb_chain = sa.Column(sa.String(SHORT), nullable=False)
    pdb_domain_def = sa.Column(sa.String(MEDIUM), nullable=False)
    pdb_pdbfam_name = sa.Column(sa.String(LONG), nullable=False)   
    pdb_pdbfam_idx = sa.Column(sa.Integer)
    domain_errors = sa.Column(sa.Text)


class DomainContact(Base):
    u"""
    Interactions between Profs domains in the PDB. Only interactions that were predicted to be
    biologically relevant by `NOXclass`_ are included in this table.
    
    Columns:
      domain_contact_id
        A unique integer identifying each domain pair.

      cath_id_1
        Unique id identifying the first interacting domain in the :ref:`domain` table.

      cath_id_2
        Unique id identifying the second interacting domain in the :ref:`domain` table.

      min_interchain_distance
        The closest that any residue in domain one comes to any residue in domain two.
        
      contact_volume
        The volume covered by contacting residues.
        
      contact_surface_area
        The surface area of the contacting regions of the first and second domains.
        
      atom_count_1
        The number of atoms in the first domain.
        
      atom_count_2
        The number of atoms in the second domain.
        
      number_of_contact_residues_1
        The number of residues in the first domain that come within 5 \u212B of the second domain.
        
      number_of_contact_residues_2
        The number of residues in the second domain that come withing 5 \u212B of the first domain.
      
      contact_residues_1
        A list of all residues in the first domain that come within 5 \u212B of the second domain.
        The residue number corresponds to the position of the residue in the domain.
      
      contact_residues_2
        A list of all residues in the second domain that come within 5 \u212B of the first domain.
        The residue number corresponds to the position of the residue in the domain.

      crystal_packing
        The probability that the interaction is a crystallization artifacts, as defined by `NOXclass`_.
      
      domain_contact_errors
        List of errors that occurred when annotating this domain pair, or when using this domain
        as a template for making structural homology models.

    .. _NOXclass: http://noxclass.bioinf.mpi-inf.mpg.de/
    """
    __tablename__ = 'domain_contact'
    _indexes = [
        (['cath_id_1', 'cath_id_2'], {'unique': True}),
        (['cath_id_2', 'cath_id_1'], {'unique': True}),
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
    
    domain_contact_id = sa.Column(sa.Integer, primary_key=True)
    cath_id_1 = sa.Column(
        None, sa.ForeignKey(Domain.cath_id, onupdate='cascade', ondelete='cascade'), nullable=False)
    cath_id_2 = sa.Column(
        None, sa.ForeignKey(Domain.cath_id, onupdate='cascade', ondelete='cascade'), nullable=False)
#    cath_id_2 = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')), 
#        nullable=False)
    min_interchain_distance = sa.Column(sa.Float)
    contact_volume = sa.Column(sa.Float)
    contact_surface_area = sa.Column(sa.Float)
    atom_count_1 = sa.Column(sa.Integer)
    atom_count_2 = sa.Column(sa.Integer)
    number_of_contact_residues_1 = sa.Column(sa.Integer)
    number_of_contact_residues_2 = sa.Column(sa.Integer)
    contact_residues_1 = sa.Column(sa.Text)
    contact_residues_2 = sa.Column(sa.Text)
    crystal_packing = sa.Column(sa.Float)
    domain_contact_errors = sa.Column(sa.Text)

    # Relationships
    domain_1 = sa.orm.relationship(
        Domain, primaryjoin=cath_id_1==Domain.cath_id, cascade='expunge', lazy='joined')
#    # the second domain may be a ligand or a peptide, and so the foreign key constraint does not work
    domain_2 = sa.orm.relationship(
        Domain, primaryjoin=cath_id_2==Domain.cath_id, cascade='expunge', lazy='joined')


class UniprotSequence(Base):
    """
    Protein sequences from the Uniprot KB, obtained by parsing ``uniprot_sprot_fasta.gz``, 
    ``uniprot_trembl_fasta.gz``, and ``homo_sapiens_variation.txt`` files from the `Uniprot ftp site`_.
    
    Columns:
      db
        The database to which the protein sequence belongs. Possible values are ``sp`` for SwissProt
        and ``tr`` for TrEMBL.

      uniprot_id
        The uniprot id of the protein.

      uniprot_name
        The uniprot name of the protein.

      protein_name
        The protein name.

      organism_name
        Name of the organism in which this protein is found.
        
      gene_name
        Name of the gene that codes for this protein sequence.

      protein_existence
        Evidence for the existence of the protein:

        1. Experimental evidence at protein level
        2. Experimental evidence at transcript level
        3. Protein inferred from homology
        4. Protein predicted
        5. Protein uncertain

      sequence_version
        Version of the protein amino acid sequence.

      uniprot_sequence
        Amino acid sequence of the protein.

    .. _Uniprot ftp site: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/

    """
    __tablename__ = 'uniprot_sequence'
    __table_args__ = get_table_args(__tablename__, [], ['uniprot_kb_schema_tuple'])

    db = sa.Column(sa.String(SHORT), nullable=False)
    uniprot_id = sa.Column(sa.String(SHORT), primary_key=True)
    uniprot_name = sa.Column(sa.String(SHORT), nullable=False)
    protein_name = sa.Column(sa.String(MEDIUM))
    organism_name = sa.Column(sa.String(MEDIUM), index=True)
    gene_name = sa.Column(sa.String(MEDIUM), index=True)
    protein_existence = sa.Column(sa.Integer)
    sequence_version = sa.Column(sa.Integer)
    uniprot_sequence = sa.Column(sa.Text, nullable=False)


class Provean(Base):
    """
    Description of the `Provean`_ supporting set calculated for a protein sequence. The construction
    of a supporting set is the most lengthy step in running Provean. Therefore, the supporting set is
    precalculated and stored for every protein sequence.

    Columns:
      uniprot_id
        The uniprot id of the protein.

      provean_supset_filename
        The filename of the Provean supporting set. The supporting set contains the ids and sequences
        of all proteins in the NCBI nr database that are used by Provean to construct a multiple
        sequence alignment for the given protein.

      provean_supset_length
        The number of sequences in Provean supporting set.

      provean_errors
        List of errors that occurred while the Provean supporting set was being calculated.

      provean_date_modified
        Date and time that this row was last modified.

    .. _provean: http://provean.jcvi.org/downloads.php
    """
    __tablename__ = 'provean'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_id = sa.Column(
        None, sa.ForeignKey(
            UniprotSequence.uniprot_id, 
            onupdate='cascade', ondelete='cascade'), 
        primary_key=True)
    provean_supset_filename = sa.Column(sa.String(MEDIUM))
    provean_supset_length = sa.Column(sa.Integer)
    provean_errors = sa.Column(sa.Text)
    provean_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)

    # Relationships
    uniprot_sequence = sa.orm.relationship(
        UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('provean', uselist=False, cascade='expunge', lazy='joined'))


class UniprotDomain(Base):
    """
    Pfam domain definitions for proteins in the :ref:`uniprot_sequence` table. This table was 
    obtained by downloading Pfam domain definitions for all known proteins from the `SIMAP`_ website,
    and mapping the protein sequence to uniprot using the MD5 hash of each sequence.

    Columns:
      uniprot_domain_id
        Unique id identifying each domain.

      uniprot_id
        The uniprot id of the protein containing the domain.

      pdbfam_name
        The Profs name of the domain. In most cases this will be equivalent to the Pfam name of the domain. 

      pdbfam_idx
        The index of the Profs domain. ``pdbfam_idx`` ranges from 1 to the number of domains with 
        the name ``pdbfam_name`` in the given protein. The ``(pdbfam_name, pdbfam_idx)`` tuple 
        uniquely identifies each domain.

      pfam_clan
        The Pfam clan to which this Profs domain belongs.

      alignment_def
        Alignment domain definitions of the Profs domain. This field is obtained by removing gaps 
        in the ``alignment_subdefs`` column.

      pfam_names
        Pfam names of all Pfam domains that were combined to create the given Profs domain.

      alignment_subdefs
        Comma-separated list of domain definitions for all Pfam domains that were merged to create 
        the given Profs domain. 

      path_to_data
        Location for storing homology models, mutation results, and all other data that are relevant
        to this domain. This path is prefixed by :term:`path_to_archive`.

    .. _SIMAP: http://liferay.csb.univie.ac.at/portal/web/simap
    """
    __tablename__ = 'uniprot_domain'
    
    IS_TRAINING_SCHEMA = 'training' in DB_SCHEMA
    
    uniprot_domain_id = sa.Column(sa.Integer, nullable=False, primary_key=True, autoincrement=True)
    uniprot_id = sa.Column(
        None, sa.ForeignKey(
            UniprotSequence.uniprot_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False)
    pdbfam_name = sa.Column(sa.String(LONG), nullable=False)
    pdbfam_idx = sa.Column(sa.Integer, nullable=False)
    pfam_clan = sa.Column(sa.Text)
    alignment_def = sa.Column(sa.String(MEDIUM))
    pfam_names = sa.Column(sa.String(LONG))
    alignment_subdefs = sa.Column(sa.Text)
    path_to_data = sa.Column(sa.Text)

    if IS_TRAINING_SCHEMA:
        # The database used for storing training data has an extra column `max_seq_identity`,
        # because we want to make homology models at different sequence identities.
        max_seq_identity = sa.Column(sa.Integer, index=True)
        _indexes = [
            (['uniprot_id', 'alignment_def', 'max_seq_identity'], 
             {'unique': True, 'index_name': 'ix_uniprot_id_unique'}),
            (['pdbfam_name'], {'mysql_length': 255}), 
        ]
    else:
        _indexes = [
            (['pdbfam_name'], {'mysql_length': 255}),   
        ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
    
    # Relationships
    uniprot_sequence = sa.orm.relationship(
        UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('uniprot_domain', cascade='expunge')) # many to one



class UniprotDomainPair(Base):
    """
    Potentially-interacting pairs of domains for proteins that are known to interact, according to 
    `Hippie`_, `IRefIndex`_, and `Rolland et al. 2014`_.

    Columns:
      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.
      
      uniprot_domain_id_1
        Unique id of the first domain.
      
      uniprot_domain_id_2
        Unique id of the second domain.
      
      rigids
        Phased out.
      
      domain_contact_ids
        List of unique ids identifying all domain-domain pairs in the PDB, where one domain
        belongs to the protein containing ``uniprot_domain_id_1`` and the other domain
        belongs to the protein containing ``uniprot_domain_id_2``. This was used as crystallographic 
        evidence that the two proteins interact.
      
      path_to_data
        Location for storing homology models, mutation results, and all other data that is relevant
        to this domain pair. This path is prefixed by :term:`path_to_archive`.
        
    .. _Hippie: http://cbdm.mdc-berlin.de/tools/hippie/
    .. _IRefIndex: http://irefindex.org
    .. _Rolland et al. 2014: http://dx.doi.org/10.1016/j.cell.2014.10.050
    """
    __tablename__ = 'uniprot_domain_pair'
    _indexes = [
            (['uniprot_domain_id_1', 'uniprot_domain_id_2'], {'unique': True}),
            (['uniprot_domain_id_2', 'uniprot_domain_id_1'], {'unique': True}),
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_domain_pair_id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    uniprot_domain_id_1 = sa.Column(
        None, sa.ForeignKey(
            UniprotDomain.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    uniprot_domain_id_2 = sa.Column(
        None, sa.ForeignKey(
            UniprotDomain.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    rigids = sa.Column(sa.Text) # Interaction references from iRefsa.Index
    domain_contact_ids = sa.Column(sa.Text) # interaction references from the PDB
    path_to_data = sa.Column(sa.Text)

    # Relationships
    uniprot_domain_1 = sa.orm.relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_1==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one
    uniprot_domain_2 = sa.orm.relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_2==UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined') # many to one


class UniprotDomainTemplate(Base):
    """
    Structural templates for domains in the :ref:`uniprot_domain` table. Lists PDB crystal structures 
    that will be used for making homology models. 
    
    Columns:
      uniprot_domain_id
        An integer which uniquely identifies each uniprot domain in the :ref:`uniprot_domain` table.
      
      template_errors
        List of errors that occurred during the process for finding the template.
      
      cath_id
        The unique id identifying the structural template of the domain. 
      
      domain_start
        The Uniprot position of the first amino acid of the Profs domain.
      
      domain_end
        The Uniprot position of the last amino acid of the Profs domain.
      
      domain_def
        Profs domain definitions for domains with structural templates. Domain definitions in this 
        column are different from domain definitions in the ``alignment_def`` column of the 
        :ref:`uniprot_domain` table in that they have been expanded to match domain boundaries of the
        Profs structural template, identified by the ``cath_id``.
      
      alignment_identity
        Percent identity of the domain to its structural template.
      
      alignment_coverage
        Percent coverage of the domain to its structural template.

      alignment_score
        A score obtained by combining ``alignment_identity`` (:math:`SeqId`) and ``alignment_coverage``
        (:math:`Cov`) using the following equation, as described by `Mosca et al.`_:
      
        .. math:: 
           :label: score_function

           Score = 0.95 \\cdot \\frac{SeqId}{100} \\cdot \\frac{Cov}{100} + 0.05 \\cdot \\frac{Cov}{100}

      t_date_modified
        The date and time when this row was last modified.
    
    .. _Mosca et al.: http://doi.org/10.1038/nmeth.2289
    """
    __tablename__ = 'uniprot_domain_template'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_domain_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomain.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    template_errors = sa.Column(sa.Text)
    cath_id = sa.Column(
        None, sa.ForeignKey(
            Domain.cath_id, 
            onupdate='cascade', ondelete='cascade'), 
        index=True, nullable=False)
    domain_start = sa.Column(sa.Integer, index=True)
    domain_end = sa.Column(sa.Integer, index=True)
    domain_def = sa.Column(sa.String(MEDIUM))
    alignment_identity = sa.Column(sa.Float)
    alignment_coverage = sa.Column(sa.Float)
    alignment_score = sa.Column(sa.Float)
    t_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)
        
    # Relationships
    uniprot_domain = sa.orm.relationship(
        UniprotDomain, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('template', uselist=False, cascade='expunge', lazy='joined')) # one to one
    domain = sa.orm.relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('uniprot_domain', cascade='expunge')) # many to one



class UniprotDomainModel(Base):
    """
    Homology models for templates in the :ref:`uniprot_domain_template` table.

    Columns:
      uniprot_domain_id
        An integer which uniquely identifies each uniprot domain in the :ref:`uniprot_domain` table.

      model_errors
        List of errors that occurred when making the homology model.

      alignment_filename
        The name of the alignment that was given to Modeller when making the homology model.

      model_filename
        The name of the homology model that was produced by Modeller.
      
      chain
        The chain that contains the domain in question in the homology (this is now set to 'A' 
        in all models).

      norm_dope
        Normalized DOPE score of the model (lower is better).

      sasa_score
        Comma-separated list of the percent solvent-accessible surface area for each residue.

      m_date_modified
        The date and time when this row was last modified.

      model_domain_def
        Domain definitions for the region of the domain that is covered by the structural template.
        
        In most cases, this field is identical to the ``domain_def`` field in the
        :ref:`uniprot_domain_template` table. However, it sometimes happens that the best 
        Profs structural template only covers a fraction of the Pfam domain. In that case, the
        ``alignment_def`` column in the :ref:`uniprot_domain` table, and the ``domain_def`` column
        in the :ref:`uniprot_domain_template` table, will contain the original Pfam domain definitions,
        and the ``model_domain_def`` column will contain domain definitions for only the region that 
        is covered by the structural template.

    """
    __tablename__ = 'uniprot_domain_model'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_domain_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainTemplate.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    model_errors = sa.Column(sa.Text)
    alignment_filename = sa.Column(sa.String(MEDIUM))
    model_filename = sa.Column(sa.String(MEDIUM))
    chain = sa.Column(sa.String(SHORT))
    norm_dope = sa.Column(sa.Float)
    sasa_score = sa.Column(sa.Text)
    m_date_modified = sa.Column(sa.DateTime, default=datetime.datetime.utcnow,
                             onupdate=datetime.datetime.utcnow, nullable=False)
    model_domain_def = sa.Column(sa.String(MEDIUM))
                             
    # Relationships
    template = sa.orm.relationship(
        UniprotDomainTemplate, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('model', uselist=False, cascade='expunge', lazy='joined')) # one to one



class UniprotDomainMutation(Base):
    """
    Characterization of mutations introduced into structures in the :ref:`uniprot_domain_model` table.
    
    Columns:
      uniprot_id
        Uniprot ID of the protein that was mutated. 

      uniprot_domain_id 
        Unique id which identifies the Profs domain that was mutated in the :ref:`uniprot_domain` table.

      mutation
        Mutation that was introduced into the protein, in Uniprot coordinates.

      mutation_errors
        List of errors that occured while evaluating the mutation.

      model_filename_wt
        The name of the file which contains the homology model of the domain after the model was
        relaxed with FoldX but before the mutation was introduced.

      model_filename_mut
        The name of the file which contains the homology model of the domain after the model was
        relaxed with FoldX and after the mutation was introduced.

      chain_modeller
        The chain which contains the domain that was mutated in the ``model_filename_wt`` and the
        ``model_filename_mut`` structures.

      mutation_modeller
        The mutation that was introduced into the protein, in PDB RESNUM coordinates.
        This identifies the mutated residue in the ``model_filename_wt`` and the
        ``model_filename_mut`` structures.

      stability_energy_wt
        Comma-separated list of scores returned by FoldX for the wildtype protein.
        The comma-separated list can be converted into a DataFrame with each column clearly labelled 
        using the :func:`elaspic.domain_mutation.format_mutation_features`.

      stability_energy_mut
        Comma-separated list of scores returned by FoldX for the mutant protein.

      physchem_wt
        Physicochemical properties describing the interaction of the wildtype residue with residues
        on the opposite chain.

      physchem_wt_ownchain
        Physicochemical properties describing the interaction of the wildtype residue with residues
        on the same chain.

      physchem_mut
        Physicochemical properties describing the interaction of the mutant residue with residues
        on the opposite chain.

      physchem_mut_ownchain
        Physicochemical properties describing the interaction of the mutant residue with residues
        on the same chain.

      matrix_score
        Score assigned to the wt -> mut transition by the BLOSUM substitution matrix.

      secondary_structure_wt
        Secondary structure of the wildtype residue predicted by `stride`_.

      solvent_accessibility_wt
        Percent solvent accessible surface area of the wildtype residue, predicted by `msms`_.

      secondary_structure_mut
        Secondary structure of the mutated residue predicted by `stride`_.

      solvent_accessibility_mut
        Percent solvent accessible surface area of the mutated residue, predicted by `msms`_.

      provean_score
        Score produced by `Provean`_ for this mutation.

      ddg
        Change in the Gibbs free energy of folding that our classifier predicts for this mutation.

      mut_date_modified
        Date and time that this row was last modified.

    .. _stride: http://webclu.bio.wzw.tum.de/stride/
    .. _msms: http://mgltools.scripps.edu/
    """
    __tablename__ = 'uniprot_domain_mutation'
    _indexes = [
        ['uniprot_id', 'mutation'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_id = sa.Column(
        None, sa.ForeignKey(
            UniprotSequence.uniprot_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    uniprot_domain_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainModel.uniprot_domain_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True, index=True)
    mutation = sa.Column(sa.String(SHORT), index=True, nullable=False, primary_key=True)
    mutation_errors = sa.Column(sa.Text)
    model_filename_wt = sa.Column(sa.String(MEDIUM))
    model_filename_mut = sa.Column(sa.String(MEDIUM))
    chain_modeller = sa.Column(sa.String(SHORT))
    mutation_modeller = sa.Column(sa.String(SHORT))
    stability_energy_wt = sa.Column(sa.Text)
    stability_energy_mut = sa.Column(sa.Text)
    physchem_wt = sa.Column(sa.Text)
    physchem_wt_ownchain = sa.Column(sa.Text)
    physchem_mut = sa.Column(sa.Text)
    physchem_mut_ownchain = sa.Column(sa.Text)
    matrix_score = sa.Column(sa.Float)
    secondary_structure_wt = sa.Column(sa.Text)
    solvent_accessibility_wt = sa.Column(sa.Float)
    secondary_structure_mut = sa.Column(sa.Text)
    solvent_accessibility_mut = sa.Column(sa.Float)
    provean_score = sa.Column(sa.Float)
    ddg = sa.Column(sa.Float, index=True)
    mut_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)
        
    # Relationships
    model = sa.orm.relationship(
        UniprotDomainModel, cascade='expunge', uselist=False, lazy='joined',
        backref=sa.orm.backref('mutations', cascade='expunge')) # many to one



class UniprotDomainPairTemplate(Base):
    """
    Structural templates for pairs of domains in the :ref:`uniprot_domain_pair` table.

    Columns:
      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.

      domain_contact_id

      cath_id_1
        Unique id of the structural template for the first domain.
      
      cath_id_2
        Unique id of the structural template for the second domain.

      identical_1
        Fraction of residues in the Blast alignment of the first domain to its template that are 
        *identical*.

      conserved_1
        Fraction of residues in the Blast alignment of the first domain to its template that are *conserved*.

      coverage_1
        Fraction of the first domain that is covered by the blast alignment.

      score_1
        Score obtained by multiplying ``identical_1`` by ``coverage_1``.

      identical_if_1
        Fraction of interface residues [#f1]_ that are *identical* in the Blast alignment of the first domain.

      conserved_if_1
        Fraction of interface residues [#f1]_ that are *conserved* in the Blast alignment of the first domain.

      coverage_if_1
        Fraction of interface residues [#f1]_ that are *covered* by the Blast alignment of the first domain.

      score_if_1
        Score obtained by combining ``identical_if_1`` and ``coverage_if_1`` using :eq:`score_function`.
    
      identical_2
        Fraction of residues in the Blast alignment of the second domain to its template that are 
        *identical*.
      
      conserved_2
        Fraction of residues in the Blast alignment of the second domain to its template that are *conserved*.
      
      coverage_2
        Fraction of the second domain that is covered by the blast alignment.
      
      score_2
        Score obtained by multiplying ``identical_2`` by ``coverage_2``.
      
      identical_if_2
        Fraction of interface residues [#f1]_ that are *identical* in the Blast alignment of the second domain.

      conserved_if_2
        Fraction of interface residues [#f1]_ that are *conserved* in the Blast alignment of the second domain.

      coverage_if_2 
        Fraction of interface residues [#f1]_ that are *covered* by the Blast alignment of the second domain.

      score_if_2 
        Score obtained by combining ``identical_if_2`` and ``coverage_if_2`` using :eq:`score_function`.

      score_total
        The product of ``score_1`` and ``score_2``.
      
      score_if_total 
        The product of ``score_if_1`` and ``score_if_2``.

      score_overall 
        The product of ``score_total`` and ``score_if_total``. This is the score that was used to 
        select the best Profs domain pair to be used as a template.
        
      t_date_modified 
        The date and time when this row was last updated.

      template_errors 
        List of errors that occured while looking for the structural template.


    .. [#f1] Interface residues are defined as residues that are within 5 \u212B of the partner domain.
    """
    __tablename__ = 'uniprot_domain_pair_template'
    _indexes = [
        ['cath_id_1', 'cath_id_2'],
        ['cath_id_2', 'cath_id_1'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_domain_pair_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainPair.uniprot_domain_pair_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    domain_contact_id = sa.Column(
        None, sa.ForeignKey(
            DomainContact.domain_contact_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False)
    cath_id_1 = sa.Column(
        None, sa.ForeignKey(
            Domain.cath_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    cath_id_2 = sa.Column(
        None, sa.ForeignKey(
            Domain.cath_id, 
            onupdate='cascade', ondelete='cascade'),
        nullable=False)

    identical_1 = sa.Column(sa.Float)
    conserved_1 = sa.Column(sa.Float)
    coverage_1 = sa.Column(sa.Float)
    score_1 = sa.Column(sa.Float)

    identical_if_1 = sa.Column(sa.Float)
    conserved_if_1 = sa.Column(sa.Float)
    coverage_if_1 = sa.Column(sa.Float)
    score_if_1 = sa.Column(sa.Float)

    identical_2 = sa.Column(sa.Float)
    conserved_2 = sa.Column(sa.Float)
    coverage_2 = sa.Column(sa.Float)
    score_2 = sa.Column(sa.Float)

    identical_if_2 = sa.Column(sa.Float)
    conserved_if_2 = sa.Column(sa.Float)
    coverage_if_2 = sa.Column(sa.Float)
    score_if_2 = sa.Column(sa.Float)

    score_total = sa.Column(sa.Float)
    score_if_total = sa.Column(sa.Float)
    score_overall = sa.Column(sa.Float)

    t_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow,
        onupdate=datetime.datetime.utcnow, nullable=False)
    template_errors = sa.Column(sa.Text)

    # Relationships
    domain_pair = sa.orm.relationship(
        UniprotDomainPair, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('template', uselist=False, cascade='expunge', lazy='joined')) # one to one
    domain_contact = sa.orm.relationship(
        DomainContact, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('uniprot', cascade='expunge')) # one to one
    domain_1 = sa.orm.relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        primaryjoin=(cath_id_1==Domain.cath_id)) # many to one
    domain_2 = sa.orm.relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        primaryjoin=(cath_id_2==Domain.cath_id)) # many to one



class UniprotDomainPairModel(Base):
    """
    Structural models of interactions between pairs of domains in the :ref:`uniprot_domain_pair`
    table.
    
    Columns:
      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.

      model_errors 
        List of errors that occured while making the homology model.

      alignment_filename_1
        Name of the file containing the alignment of the first domain with its structural template.

      alignment_filename_2
        Name of the file containing the alignment of the second domain with its structural template.

      model_filename
        Name of the file containing the homology model of the domain-domain interaction 
        created by Modeller.

      chain_1 
        Chain containing the first domain in the model specified by ``model_filename``.

      chain_2
        Chain containing the second domain in the model specified by ``model_filename``.

      norm_dope
        The normalized DOPE score of the model.

      interface_area_hydrophobic
        Hydrophobic surface area of the interface, calculated using `POPS`_.

      interface_area_hydrophilic
        Hydrophilic surface area of the interface, calculated using `POPS`_.

      interface_area_total
        Total surface area of the interface, calculated using `POPS`_.

      interface_dg
        Gibbs free energy of binding for this domain-domain interaction, predicted using `FoldX`_.
        Not implemented yet!

      interacting_aa_1
        List of amino acid positions in the first domain that are within 5 \u212B of the second domain. 
        Positions are specified using uniprot coordinates.

      interacting_aa_2
        List of amino acids in the second domain that are within 5 \u212B of the first domain.
        Position are specified using uniprot coordinates.

      m_date_modified 
        Date and time that this row was last modified.

      model_domain_def_1
        Domain boundaries of the first domain that are covered by the Profs structural template.
      
      model_domain_def_2
        Domain boundaries of the second domain that are covered by the Profs structural template.


    .. _POPS: http://mathbio.nimr.mrc.ac.uk/wiki/POPS
    .. _FoldX: http://foldx.crg.es/
    """
    __tablename__ = 'uniprot_domain_pair_model'
    __table_args__ = get_table_args(__tablename__, [], ['schema_version_tuple'])

    uniprot_domain_pair_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainPairTemplate.uniprot_domain_pair_id, 
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    model_errors = sa.Column(sa.Text)
    alignment_filename_1 = sa.Column(sa.String(MEDIUM))
    alignment_filename_2 = sa.Column(sa.String(MEDIUM))
    model_filename = sa.Column(sa.String(MEDIUM))
    chain_1 = sa.Column(sa.String(SHORT))
    chain_2 = sa.Column(sa.String(SHORT))
    norm_dope = sa.Column(sa.Float)
    interface_area_hydrophobic = sa.Column(sa.Float)
    interface_area_hydrophilic = sa.Column(sa.Float)
    interface_area_total = sa.Column(sa.Float)
    interface_dg = sa.Column(sa.Float)
    interacting_aa_1 = sa.Column(sa.Text)
    interacting_aa_2 = sa.Column(sa.Text)
    m_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow, 
        nullable=False)
    model_domain_def_1 = sa.Column(sa.String(MEDIUM))
    model_domain_def_2 = sa.Column(sa.String(MEDIUM))
        
    # Relationships
    template = sa.orm.relationship(
        UniprotDomainPairTemplate, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('model', uselist=False, cascade='expunge', lazy='joined')) # one to one
        


class UniprotDomainPairMutation(Base):
    u"""
    Characterization of interface mutations introduced into structures in the :ref:`uniprot_domain_pair_model` table.

    Columns:
      uniprot_id 
        Uniprot ID of the protein that is being mutated.

      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.

      mutation
        Mutation for which the :math:`\Delta \Delta G` score is being predicted, specified in Uniprot coordinates.

      mutation_errors
        List of errors obtained when evaluating the impact of the mutation.

      model_filename_wt
        Filename of the homology model relaxed by FoldX but containing the wildtype residue.

      model_filename_mut
        Filename of the homology model relaxed by FoldX and containing the mutated residue.

      chain_modeller
        Chain containing the domain that was mutated, in homology models specified by
        ``model_filename_wt`` and ``model_filename_mut``.

      mutation_modeller
        Mutation for which the :math:`\Delta \Delta G` score is being predicted, specified in PDB RESNUM coordinates.

      analyse_complex_energy_wt
        Comma-separated list of FoldX scores describing the effect of the wildtype residue on 
        the stability of the protein domain.

      stability_energy_wt
        Comma-separated list of FoldX scores describing the effect of the wildtype residue on 
        protein-protein interaction interface.

      analyse_complex_energy_mut
        Comma-separated list of FoldX scores describing the effect of the mutated residue on 
        the stability of the protein domain.

      stability_energy_mut
        Comma-separated list of FoldX scores describing the effect of the mutated residue on 
        protein-protein interaction interface.

      physchem_wt
        Comma-separated list of physicochemical properties describing the interaction between 
        the wildtype residue and other residues on the opposite chain.

      physchem_wt_ownchain
        Comma-separated list of physicochemical properties describing the interaction between 
        the wildtype residue and other residues on the same chain.

      physchem_mut
        Comma-separated list of physicochemical properties describing the interaction between 
        the mutated residue and other residues on the opposite chain.

      physchem_mut_ownchain
        Comma-separated list of physicochemical properties describing the interaction between 
        the mutated residue and other residues on the same chain.

      matrix_score
        Score assigned to the wt -> mut transition by the BLOSUM substitution matrix.

      secondary_structure_wt
        Secondary structure of the wildtype residue, predicted by `stride`_.
      
      solvent_accessibility_wt
        Percent solvent accessible surface area of the wildtype residue, predicted by `msms`_.

      secondary_structure_mut
        Secondary structure of the mutated residue, predicted by `stride`_.

      solvent_accessibility_mut
        Percent solvent accessible surface area of the mutated residue, predicted by `msms`_.

      contact_distance_wt
        Shortest distance between the wildtype residue and a residue on the opposite chain.
      
      contact_distance_mut
        Shortest distance between the mutated reside and a residue on the opposite chain.

      provean_score
        `Provean`_ score for this mutation.

      ddg 
        Predicted change in Gibbs free energy of binding caused by this mutation.

      mut_date_modified
        Date and time when this row was last modified.
    """
    __tablename__ = 'uniprot_domain_pair_mutation'
    _indexes = [
        ['uniprot_id', 'mutation'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])

    uniprot_id = sa.Column(None, sa.ForeignKey(
        UniprotSequence.uniprot_id, onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    uniprot_domain_pair_id = sa.Column(None, sa.ForeignKey(
        UniprotDomainPairModel.uniprot_domain_pair_id, onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    mutation = sa.Column(sa.String(SHORT),
        nullable=False, primary_key=True)
    mutation_errors = sa.Column(sa.Text)
    model_filename_wt = sa.Column(sa.String(MEDIUM))
    model_filename_mut = sa.Column(sa.String(MEDIUM))
    chain_modeller = sa.Column(sa.String(SHORT))
    mutation_modeller = sa.Column(sa.String(SHORT))
    analyse_complex_energy_wt = sa.Column(sa.Text)
    stability_energy_wt = sa.Column(sa.Text)
    analyse_complex_energy_mut = sa.Column(sa.Text)
    stability_energy_mut = sa.Column(sa.Text)
    physchem_wt = sa.Column(sa.Text)
    physchem_wt_ownchain = sa.Column(sa.Text)
    physchem_mut = sa.Column(sa.Text)
    physchem_mut_ownchain = sa.Column(sa.Text)
    matrix_score = sa.Column(sa.Float)
    secondary_structure_wt = sa.Column(sa.Text)
    solvent_accessibility_wt = sa.Column(sa.Float)
    secondary_structure_mut = sa.Column(sa.Text)
    solvent_accessibility_mut = sa.Column(sa.Float)
    contact_distance_wt = sa.Column(sa.Float)
    contact_distance_mut = sa.Column(sa.Float)
    provean_score = sa.Column(sa.Float)
    ddg = sa.Column(sa.Float, index=False)
    mut_date_modified = sa.Column(sa.DateTime, default=datetime.datetime.utcnow,
                               onupdate=datetime.datetime.utcnow, nullable=False)
    # Relationships
    model = sa.orm.relationship(
        UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('mutations', cascade='expunge')) # many to one

         
         
#%%
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
        elif self.configs['db_type'] in ['postgresql', 'mysql']:
            autocommit = False
            autoflush = False
        Session.configure(bind=self.engine, autocommit=autocommit, autoflush=autoflush)
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


    def get_uniprot_domain_pair(self, uniprot_id, copy_data=False):
        """
        """
        with self.session_scope() as session:
            uniprot_domain_pairs = (
                session.query(UniprotDomainPair)
                .filter(sa.or_(
                    sa.text("uniprot_domain_1.uniprot_id='{}'".format(uniprot_id)),
                    sa.text("uniprot_domain_2.uniprot_id='{}'".format(uniprot_id))))
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
            except subprocess.CalledProcessError as e:
                self.logger.error(e)
                uniprot_mutation.model_filename_wt = None
        return uniprot_mutation


    def _copy_mutation_data(self, mutation, path_to_data, path_to_archive):
        if mutation and mutation.model_filename_wt:
            tmp_save_path = self.configs['temp_path'] + path_to_data
            archive_save_path = path_to_archive + path_to_data
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
