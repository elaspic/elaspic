# -*- coding: utf-8 -*-
"""
Naming convention for classes includes an underscore since peewee converts 
database table names to lowercase.
"""
import os
import pandas as pd
import urllib2
import subprocess
import json
from string import uppercase

from sqlalchemy import create_engine, and_, or_
from sqlalchemy import Column, Index, UniqueConstraint
from sqlalchemy import Integer, Float, String, Boolean, Text
from sqlalchemy import ForeignKey
from sqlalchemy.orm import sessionmaker, relationship, backref, aliased
from sqlalchemy.ext.declarative import declarative_base

import Bio

import parse_pfamscan
import class_error as error


###############################################################################

# Not the best place to define this, bot collation changes depending on the 
# database type...

sql_flavor = 'sqlite_file'
#sql_flavor = 'postgresql'

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
Session = sessionmaker()

###############################################################################
# Helper functions for dealing with sql objects

def decode_domain(domains, merge=True):
    """ Unlike split_domain(), this function returns a tuple of tuples of strings,
    preserving letter numbering (e.g. 10B)
    """
    if domains[-1] == ',':
        domains = domains[:-1]
    x = domains
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

    return d

#############################################################################

class Domain(Base):
    __tablename__ = 'domain'
    
    cath_id = Column(String(15, collation=binary_collation), primary_key=True)

    pdb_id = Column(String(4, collation=string_collation), nullable=False)
    pdb_type = Column(String(31, collation=string_collation), nullable=True)
    pdb_resolution = Column(Float, nullable=True)
    pdb_chain = Column(String(1, collation=string_collation), nullable=False)
    pdb_domain_def = Column(String(255, collation=string_collation), nullable=False)
    
    pfam_autopfam = Column(Integer, nullable=False)
    pfam_name = Column(String(255, collation=string_collation), nullable=False)


class DomainContact(Base):
    __tablename__ = 'domain_contact'
    __table_args__ = ({'sqlite_autoincrement': True},)
    # Columns
    domain_contact_id = Column(Integer, primary_key=True)
    
    cath_id_1 = Column(None, ForeignKey('domain.cath_id'), index=True, nullable=False)
    pdb_contact_residues_1 = Column(Text, nullable=False)

    cath_id_2 = Column(None, ForeignKey('domain.cath_id'), index=True, nullable=False)
    pdb_contact_residues_2 = Column(Text, nullable=False)
    
    Index('cath_id_1_2', 'cath_id_1', 'cath_id_1', unique=True)
    # Relationships
    domain_1 = relationship("Domain", primaryjoin=cath_id_1==Domain.cath_id)
    domain_2 = relationship("Domain", primaryjoin=cath_id_2==Domain.cath_id)
    
###############################################################################

class UniprotSequence(Base):
    __tablename__ = 'uniprot_sequence'
    
    uniprot_id = Column(String(50, collation=string_collation), primary_key=True, nullable=False)
    uniprot_name = Column(String(255, collation=string_collation))
    uniprot_description = Column(String(255, collation=string_collation))
    uniprot_sequence = Column(Text, nullable=False)


class UniprotDomain(Base):
    __tablename__ = 'uniprot_domain'
    __table_args__ = (UniqueConstraint('uniprot_id', 'pfam_name', 'alignment_def', name='uniprot_domain_tuple'),)
    
    uniprot_domain_id = Column(Integer, primary_key=True)
    
    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False)
    pfam_name = Column(String(255, collation=string_collation), index=True, nullable=False) #seq_id
    alignment_def = Column(String(255, collation=string_collation), nullable=False)
    
    alignment_defs = Column(Text)
    hmm_accession = Column(Text)
    hmm_name = Column(Text)
    hmm_type = Column(Text)
    hmm_clan_accession = Column(Text)
    hmm_clan_name = Column(Text)
    
    path_to_data = Column(Text)
    
    # Relationships
    uniprot_sequence = relationship(UniprotSequence, backref='domain') # many to one


class UniprotDomainPair(Base):
    __tablename__ = 'uniprot_domain_pair'
    __table_args__ = (Index('uniprot_domain_id_1_2', "uniprot_domain_id_1", "uniprot_domain_id_2", unique=True),
                      {'sqlite_autoincrement': True},)
                      
    uniprot_domain_pair_id = Column(Integer, primary_key=True)
    
    uniprot_domain_id_1 = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), index=True, nullable=False)
    uniprot_domain_id_2 = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), index=True, nullable=False)
    
    path_to_data = Column(Text)
    
    # Relationships
    uniprot_domain_1 = relationship(UniprotDomain, primaryjoin=uniprot_domain_id_1==UniprotDomain.uniprot_domain_id) # many to one
    uniprot_domain_2 = relationship(UniprotDomain, primaryjoin=uniprot_domain_id_2==UniprotDomain.uniprot_domain_id) # many to one
    

###############################################################################

class UniprotDomainTemplate(Base):
    __tablename__ = 'uniprot_domain_template'
    
    uniprot_domain_id = Column(None, ForeignKey(UniprotDomain.uniprot_domain_id), primary_key=True)
    template_errors = Column(Text)
    alignment_filename = Column(String(255, collation=string_collation))
    
    #
    cath_id = Column(None, ForeignKey(Domain.cath_id), index=True)
    domain_def = Column(String(255, collation=string_collation))
    alignment_id = Column(String(255, collation=string_collation))
    alignment_score = Column(Float)
    
    # Relationships
    uniprot_domain = relationship(UniprotDomain, uselist=False, backref='template') # one to one
    domain = relationship(Domain)


class UniprotDomainModel(Base):
    __tablename__ = 'uniprot_domain_model'
    
    uniprot_domain_id = Column(None, ForeignKey(UniprotDomainTemplate.uniprot_domain_id), primary_key=True)
    model_errors = Column(Text)
    model_filename = Column(String(255, collation=string_collation))
    
    chain = Column(String(1, collation=string_collation))
    norm_dope = Column(Float)
    het_flag = Column(Boolean)
    switch_chain = Column(Boolean)
    surface_aa = Column(Text)
    
    # Relationships
    template = relationship(UniprotDomainTemplate, uselist=False, backref='model') # one to one


class UniprotDomainMutation(Base):
    __tablename__ = 'uniprot_domain_mutation'
    
    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False, primary_key=True)
    uniprot_domain_id = Column(None, ForeignKey(UniprotDomainModel.uniprot_domain_id), primary_key=True)
    mutation = Column(String(11, collation=string_collation), nullable=False, primary_key=True)
    mutation_errors = Column(Text)
    
    model_filename_wt = Column(String(255, collation=string_collation))
    model_filename_mut = Column(String(255, collation=string_collation))
    #
    mutation_pdb = Column(Integer)
    mutation_modeller = Column(Integer)
    
    chain_pdb = Column(String(1, collation=string_collation))
    chain_modeller = Column(String(1, collation=string_collation))
    
    AnalyseComplex_energy_wt = Column(Float)
    Stability_energy_wt = Column(Float)
    AnalyseComplex_energy_mut = Column(Float)
    Stability_energy_mut = Column(Float)
    
    interface_size = Column(Float)
    
    physChem_wt = Column(Float)
    physChem_wt_ownChain = Column(Float)
    physChem_mut = Column(Float)
    physChem_mut_ownChain = Column(Float)

    matrix_score = Column(Float)
    
    secondary_structure_wt = Column(Float)
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(Float)
    solvent_accessibility_mut = Column(Float)
           
    ddG = Column(Float)
    
    # Relationships
    model = relationship(UniprotDomainModel, backref='mutations') # many to one

    
###############################################################################

class UniprotDomainPairTemplate(Base):
    __tablename__ = 'uniprot_domain_pair_template'
    
    # Columns
    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPair.uniprot_domain_pair_id), primary_key=True)
    template_errors = Column(Text)
    alignment_filename_1 = Column(String(255, collation=string_collation))
    alignment_filename_2 = Column(String(255, collation=string_collation))
    
#    domain_contact_id = Column(None, ForeignKey(DomainContact.domain_contact_id), index=True)        
    cath_id_1 = Column(None, ForeignKey(Domain.cath_id), index=True)
    domain_def_1 = Column(String(255, collation=string_collation))
    alignment_id_1 = Column(String(255, collation=string_collation))
    alignment_score_1 = Column(Float)

    cath_id_2 = Column(None, ForeignKey(Domain.cath_id), index=True)
    domain_def_2 = Column(String(255, collation=string_collation))
    alignment_id_2 = Column(String(255, collation=string_collation))
    alignment_score_2 = Column(Float)
    
    # Relationships
    uniprot_domain_pair = relationship(UniprotDomainPair, uselist=False, backref='template') # one to one
    domain_1 = relationship(Domain, primaryjoin=cath_id_1==Domain.cath_id)
    domain_2 = relationship(Domain, primaryjoin=cath_id_2==Domain.cath_id)
    

class UniprotDomainPairModel(Base):
    __tablename__ = 'uniprot_domain_pair_model'
    
    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPairTemplate.uniprot_domain_pair_id), primary_key=True)
    model_errors = Column(Text)
    model_filename = Column(String(255, collation=string_collation))
    
    #
    chain_1 = Column(String(1, collation=string_collation))
    chain_2 = Column(String(1, collation=string_collation))

    norm_dope = Column(Float)
    het_flag = Column(Boolean)
    switch_chain = Column(Boolean)
    
    interface_aa_1 = Column(Text)
    interface_aa_2 = Column(Text)
    
    # Relationships
    template = relationship(UniprotDomainPairTemplate, uselist=False, backref='model') # one to one


class UniprotDomainPairMutation(Base):
    __tablename__ = 'uniprot_domain_pair_mutation'
    
    uniprot_id = Column(None, ForeignKey(UniprotSequence.uniprot_id), index=True, nullable=False, primary_key=True)
    uniprot_domain_pair_id = Column(None, ForeignKey(UniprotDomainPairModel.uniprot_domain_pair_id), primary_key=True)
    mutation = Column(String(11, collation=string_collation), nullable=False, primary_key=True)
    mutation_errors = Column(Text)

    model_filename_wt = Column(String(255, collation=string_collation))
    model_filename_mut = Column(String(255, collation=string_collation))    
    #
    mutation_position_pdb = Column(Integer)
    mutation_position_modeller = Column(Integer)
    
    chain_pdb = Column(String(1, collation=string_collation))
    chain_modeller = Column(String(1, collation=string_collation))
    
    AnalyseComplex_energy_wt = Column(Float)
    Stability_energy_wt = Column(Float)
    AnalyseComplex_energy_mut = Column(Float)
    Stability_energy_mut = Column(Float)
    
    interface_size = Column(Float)
    
    physChem_wt = Column(Float)
    physChem_wt_ownChain = Column(Float)
    physChem_mut = Column(Float)
    physChem_mut_ownChain = Column(Float)

    matrix_score = Column(Float)
    
    secondary_structure_wt = Column(Float)
    solvent_accessibility_wt = Column(Float)
    secondary_structure_mut = Column(Float)
    solvent_accessibility_mut = Column(Float)
    
    ddG = Column(Float)
    
    # Relationships
    model = relationship(UniprotDomainPairModel, backref='mutations') # many to one

###############################################################################
class MyDatabase(object):
    """
    """
    
    def __init__(self, path_to_sqlite_db='', sql_flavor='sqlite_file', is_immutable=False,
                 path_to_temp='/tmp/', path_to_archive='/tmp/human/', clear_database=False):
        """
        """
        # Based on configuration, use a different database
        # (Defaults to an in-memmory sqlite engine)
        if sql_flavor == 'sqlite':
            engine = create_engine('sqlite://')
            if clear_database:
                Base.metadata.drop_all(engine)
            Base.metadata.create_all(engine)
        elif sql_flavor == 'sqlite_file':
            engine = create_engine('sqlite:///' + path_to_sqlite_db)
            if clear_database:
                Base.metadata.drop_all(engine)            
            Base.metadata.create_all(engine)
        elif sql_flavor == 'mysql':
            engine = create_engine('mysql://elaspic:elaspic@192.168.6.19:3306/elaspic')
            if clear_database:
                Base.metadata.drop_all(engine)            
            Base.metadata.create_all(engine)
        elif sql_flavor == 'mysql':
            engine = create_engine('mysql://root:kim630@127.0.0.1:3306/elaspic')
            if clear_database:
                Base.metadata.drop_all(engine)            
            Base.metadata.create_all(engine)
        elif sql_flavor == 'postgresql':
            engine = create_engine('postgresql://elaspic:elaspic@192.168.6.19:5432/elaspic')
            if clear_database:
                Base.metadata.drop_all(engine)
            Base.metadata.create_all(engine)
            
        Session.configure(bind=engine)
        self.session = Session()
        self.is_immutable = is_immutable
        self.path_to_temp = path_to_temp
        self.path_to_archive = path_to_archive
    
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
#                session.flush()
                print idx
        self.session.commit()
        print 'Finished populating table domain'
        
        
        # Table `domain_contact`
        names = ['domain_contact_id', 'cath_id_1', 'pdb_contact_residues_1', 'cath_id_2', 'pdb_contact_residues_2']
        domain_contact_df = pd.read_csv(domain_contact_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
        domain_contact_df = domain_contact_df.dropna() # only a couple of rows are droppeds
        for idx, row in domain_contact_df.iterrows():
            self.session.add(DomainContact(**row.to_dict()))
            if idx % 10000 == 0:
#                session.flush()
                print idx
        self.session.commit()
        print 'Finished populating table domain_contact'
        
        
        # Table `uniprot_sequence`
        names = ['uniprot_id', 'uniprot_name', 'uniprot_description', 'uniprot_sequence']
        uniprot_sequence_df = pd.read_csv(uniprot_sequence_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
        for idx, row in uniprot_sequence_df.iterrows():
            self.session.add(UniprotSequence(**row.to_dict()))
            if idx % 10000 == 0:
#                session.flush()
                print idx
        self.session.commit()
        print 'Finished populating table uniprot_sequence'     
        
        
        # Table `uniprot_domain`
        if os.path.isfile(uniprot_domain_infile):
#            names = ['uniprot_domain_id', 'uniprot_id', 'pfam_name', 'alignment_def']
            uniprot_domain_df_with_id = pd.read_csv(uniprot_domain_infile, sep='\t', na_values='\N', index_col=False)
            uniprot_domain_df_with_id['alignment_defs'] = uniprot_domain_df_with_id['alignment_def']
            uniprot_domain_df_with_id['alignment_def'] = uniprot_domain_df_with_id['alignment_defs'].apply(lambda x: encode_domain(decode_domain(x)))
            uniprot_domain_df_with_id['path_to_data'] = (uniprot_domain_df_with_id['uniprot_id'] + '/' +
                                                        uniprot_domain_df_with_id['pfam_name'] + '*' +
                                                        uniprot_domain_df_with_id['alignment_def'].apply(lambda x: x.replace(':','-')) + '/')
            for idx, row in uniprot_domain_df_with_id.iterrows():
                self.session.add(UniprotDomain(**row.to_dict()))
                if idx % 10000 == 0:
#                    session.flush()
                    print idx
            self.session.commit()
            print 'Finished populating table uniprot_domain'
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
            temp['path_to_data'] = (temp['uniprot_id_1'] + '/' +
                                    temp['pfam_name_1'] + '*' +
                                    temp['alignment_def_1'].apply(lambda x: x.replace(':','-')) + '/' +
                                    temp['pfam_name_2'] + '*' +
                                    temp['alignment_def_2'].apply(lambda x: x.replace(':','-')) + '/' +
                                    temp['uniprot_id_2'] + '/')
            temp = temp[['uniprot_domain_pair_id', 'path_to_data']]
            uniprot_domain_pair_df_with_id = uniprot_domain_pair_df_with_id.merge(temp, how='left')

#            uniprot_domain_pair_df_with_id.to_sql('uniprot_domain_pair', conn, flavor=sql_flavor, if_exists='append')
            for idx, row in uniprot_domain_pair_df_with_id.iterrows():
                self.session.add(UniprotDomainPair(**row.to_dict()))
                if idx % 10000 == 0:
#                   self.session.flush()
                    print idx
            self.session.commit()
            print 'Finished populating table domain'  
        else:
            pass
#            pfam_parser = parse_pfamscan.make_uniprot_domain_pair_database(domain_df, domain_contact_df, uniprot_domain_df_with_id,
#                                infile='/home/kimlab1/strokach/working/databases/biogrid/pairs_of_interacting_uniprots_human.tsv')
#            uniprot_domain_pair_df = pfam_parser.get_dataframe()
##            uniprot_domain_pair_df.to_sql('uniprot_domain_pair', conn, flavor=sql_flavor, if_exists='append')
##            uniprot_domain_pair_df_with_id = pd.read_sql('SELECT * from uniprot_domain_pair', conn)
#            uniprot_domain_pair_df_with_id.to_csv(uniprot_domain_pair_infile, sep='\t', na_values='\N', index=False)


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
        
        uniprot_sequence = self.session.query(UniprotSequence)\
                            .filter(UniprotSequence.uniprot_id==uniprot_id)\
                            .all()
        
        if len(uniprot_sequence) == 1:
            uniprot_sequence = uniprot_sequence[0]        
        
        elif len(uniprot_sequence) == 0:
            # the True/False value is used to add new sequences to the database in
            # end. Only usefull if you run one instance at a time otherwise you will
            # get diverging sequence databases.
            childProcess = subprocess.Popen('whoami', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            whoami, e = childProcess.communicate()
            if whoami.strip() == 'joan':
                print 'uniprot sequence not found'
            else:
                print 'Fetching uniprot sequence', uniprot_id, 'from server'
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
        
        uniprot_seqio_oject = Bio.SeqIO.SeqRecord(seq=Bio.Seq.Seq(str(uniprot_sequence.uniprot_sequence)), 
                               id=uniprot_sequence.uniprot_id,
                               name=uniprot_sequence.uniprot_name,
                               description=uniprot_sequence.uniprot_description)
        
        return uniprot_seqio_oject
    
    
    def add_uniprot_sequence(self, uniprot_sequence):
        """ Add the new items (which is a list of tuples) to the database
        """
        if not self.is_immutable:
            self.session.add(uniprot_sequence)
            self.session.commit()
        
    
    ###########################################################################
    def get_domain(self, pfam_name):
        """ Contains pdbfam-based definitions of all pfam domains in the pdb
        """
        
        domain = self.session.query(Domain)\
                    .filter(Domain.pfam_name==pfam_name)\
                    .all()

        if domain == []:
            print 'No domain template found for domain %s' % pfam_name
            
        return domain


    def get_domain_contact(self, pfam_name_1, pfam_name_2):
        """ Keeps the domain-domain interaction information from pdbfam
        Note that the produced dataframe may not have the same order as the keys
        """
        domain_contact_1 = self._get_domain_contact(pfam_name_1, pfam_name_2, reverse=False)
        domain_contact_2 = self._get_domain_contact(pfam_name_1, pfam_name_2, reverse=True)
            
        if len(domain_contact_1)==0 and len(domain_contact_2)==0:
            print 'No domain contact template found for domains %s, %s' % (pfam_name_1, pfam_name_2,)
            
        return [domain_contact_1, domain_contact_2]
        
        
    def _get_domain_contact(self, pfam_name_1, pfam_name_2, reverse=False):
        """        
        """  
        if reverse:
            pfam_name_1, pfam_name_2 = pfam_name_2, pfam_name_1

        domain_1 = aliased(Domain)
        domain_2 = aliased(Domain)
        
        domain_contact = self.session\
            .query(DomainContact)\
            .join(domain_1, DomainContact.cath_id_1==domain_1.cath_id)\
            .filter(domain_1.pfam_name==pfam_name_1)\
            .join(domain_2, DomainContact.cath_id_2==domain_2.cath_id)\
            .filter(domain_2.pfam_name==pfam_name_2)\
            .all()
        
        return domain_contact


    ###########################################################################
    def get_uniprot_domain(self, uniprot_id):
        """ Initiated using parsed pfam_scan.pl output for the human uniprot (or the entire uniprot)
        The database is updated with information about calculated models
        """
        uniprot_definitions = self.session\
            .query(UniprotDomain, UniprotDomainTemplate, UniprotDomainModel)\
            .filter(UniprotDomain.uniprot_id==uniprot_id)\
            .outerjoin(UniprotDomainTemplate)\
            .outerjoin(UniprotDomainModel)\
            .all()
        
        if len(uniprot_definitions)==0:
            print 'No domain found in uniprot %s' % uniprot_id
        
        for uniprot_domain, uniprot_template, uniprot_model in uniprot_definitions:
            tmp_save_path = self.path_to_temp + uniprot_domain.path_to_data
            archive_save_path = self.path_to_archive + uniprot_domain.path_to_data      
            if uniprot_template:
                subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
                subprocess.check_call('cp ' + archive_save_path + uniprot_template.alignment_filename +
                                        ' ' + tmp_save_path + uniprot_template.alignment_filename, shell=True)
            if uniprot_model:
                subprocess.check_call('cp ' +  + archive_save_path + uniprot_model.model_filename +
                                        ' ' + tmp_save_path + uniprot_model.model_filename, shell=True)
        return uniprot_definitions
        
    
    def get_uniprot_domain_mutation(self, uniprot_domain_id, mutation, path_to_data=False):
        """
        """
        uniprot_mutation = self.session\
            .query(UniprotDomainMutation)\
            .filter(UniprotDomainMutation.uniprot_domain_id==uniprot_domain_id)\
            .filter(UniprotDomainMutation.mutation==mutation)\
            .all()
            
        if len(uniprot_mutation) == 0:
            print 'No precalculated mutation %s for uniprot domain number %s' % (mutation, uniprot_domain_id)
        
#        if path_to_data:
#            tmp_save_path = self.path_to_temp + path_to_data
#            archive_save_path = self.path_to_archive + path_to_data
        
        return uniprot_mutation
        

    ###########################################################################
    def get_uniprot_domain_pair(self, uniprot_id):
        """ Contains known interactions between uniprot proteins
        Checks if the interaction database is already available as a pickled object.
        Generates it from file 'textfile_name' otherwise.
        """
        
        uniprot_domain_pair_1 = self._get_uniprot_domain_pair(uniprot_id, reverse=False)
        uniprot_domain_pair_2 = self._get_uniprot_domain_pair(uniprot_id, reverse=True)
        
        if len(uniprot_domain_pair_1)==0 and len(uniprot_domain_pair_2)==0:
            print 'No known interactions with uniprot %s' % uniprot_id
        
        for uniprot_domain, uniprot_template, uniprot_model in uniprot_domain_pair_1 + uniprot_domain_pair_2:
            tmp_save_path = self.path_to_temp + uniprot_domain.path_to_data
            archive_save_path = self.path_to_archive + uniprot_domain.path_to_data
            if uniprot_template:
                subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
                subprocess.check_call('cp ' + archive_save_path + uniprot_template.alignment_filename_1 +
                                        ' ' + tmp_save_path + uniprot_template.alignment_filename_1, shell=True)                    
                subprocess.check_call('cp ' + archive_save_path + uniprot_template.alignment_filename_2 +
                                        ' ' + tmp_save_path + uniprot_template.alignment_filename_2, shell=True)
            if uniprot_model:
                subprocess.check_call('cp ' + archive_save_path + uniprot_model.model_filename +
                                        ' ' + tmp_save_path + uniprot_model.model_filename, shell=True)                    
        
        return [uniprot_domain_pair_1, uniprot_domain_pair_2]


    def _get_uniprot_domain_pair(self, uniprot_id_1, reverse=False):
        
        domain_id_1 = UniprotDomainPair.uniprot_domain_id_1
        domain_id_2 = UniprotDomainPair.uniprot_domain_id_2
        if reverse:
            domain_id_1, domain_id_2 = domain_id_2, domain_id_1
            
        uniprot_domain_1 = aliased(UniprotDomain)
        uniprot_domain_2 = aliased(UniprotDomain)
        
        uniprot_domain_pair = self.session\
            .query(UniprotDomainPair, UniprotDomainPairTemplate, UniprotDomainPairModel)\
            .join(uniprot_domain_1, domain_id_1==uniprot_domain_1.uniprot_domain_id)\
            .filter(uniprot_domain_1.uniprot_id==uniprot_id_1)\
            .join(uniprot_domain_2, domain_id_2==uniprot_domain_2.uniprot_domain_id)\
            .outerjoin(UniprotDomainPairTemplate)\
            .outerjoin(UniprotDomainPairModel)\
            .all()
        
        return uniprot_domain_pair
    
    
    def get_uniprot_domain_pair_mutation(self, uniprot_domain_pair_id, mutation):
        """
        """
        position = int(mutation[1:-1])
        from_aa = mutation[0]
        to_aa = mutation[-1]

        uniprot_mutation = self.session\
            .query(UniprotDomainPairMutation)\
            .filter(UniprotDomainMutation.uniprot_domain_pair_id==uniprot_domain_pair_id)\
            .filter(UniprotDomainMutation.mutation_position_uniprot==position)\
            .filter(UniprotDomainMutation.mutation_from==from_aa)\
            .filter(UniprotDomainMutation.mutation_to==to_aa)\
            .all()
            
        if len(uniprot_mutation) == 0:
            print 'No precalculated mutation %s for uniprot domain pair number %s' % (mutation, uniprot_domain_pair_id)
            
#        if path_to_data:
#            tmp_save_path = self.path_to_temp + path_to_data
#            archive_save_path = self.path_to_archive + path_to_data
            
        return uniprot_mutation

        
    ###########################################################################
    def add_uniprot_template(self, uniprot_template, path_to_data=False):
        
        # Save a copy of the alignment to the export folder
        if path_to_data:
            tmp_save_path = self.path_to_temp + path_to_data 
            archive_save_path = self.path_to_archive + path_to_data
            subprocess.check_call('mkdir -p ' + archive_save_path, shell=True)
            if type(uniprot_template) == UniprotDomainTemplate:
                subprocess.check_call('cp ' + tmp_save_path + uniprot_template.alignment_filename +
                                        ' ' + archive_save_path + uniprot_template.alignment_filename, shell=True)
            elif type(uniprot_template) == UniprotDomainPairTemplate:
                subprocess.check_call('cp ' + tmp_save_path + uniprot_template.alignment_filename_1 +
                                        ' ' + archive_save_path + uniprot_template.alignment_filename_1, shell=True)
                subprocess.check_call('cp ' + tmp_save_path + uniprot_template.alignment_filename_2 +
                                        ' ' + archive_save_path + uniprot_template.alignment_filename_2, shell=True)
            else:
                raise Exception('Wrong uniprot_template type!')
            with open(archive_save_path + 'template.json', 'w') as fh:
                json.dump(row2dict(uniprot_template), fh, indent=4, separators=(',', ': '))
        
        if not self.is_immutable:
            self.session.add(uniprot_template)
            self.session.commit()
    
    
    def add_uniprot_model(self, uniprot_model, path_to_data=False):
        
        # Save a copy of the alignment to the export folder
        if path_to_data:
            tmp_save_path = self.path_to_temp + path_to_data 
            archive_save_path = self.path_to_archive + path_to_data
            subprocess.check_call('mkdir -p ' + archive_save_path, shell=True)
            subprocess.check_call('cp ' + tmp_save_path + uniprot_model.model_filename +
                                    ' ' + archive_save_path + uniprot_model.model_filename, shell=True)
            with open(archive_save_path + 'model.json', 'w') as fh:
                json.dump(row2dict(uniprot_model), fh, indent=4, separators=(',', ': '))
        
        if not self.is_immutable:
            self.session.add(uniprot_model)
            self.session.commit()
    
    
    def add_uniprot_mutation(self, uniprot_mutation, path_to_data=False):
        
        if path_to_data:
            tmp_save_path = self.path_to_temp + path_to_data
            archive_save_path = self.path_to_archive + path_to_data
            archive_save_subpath = uniprot_mutation.model_filename_wt.split('/')[0] + '/'
            subprocess.check_call('mkdir -p ' + archive_save_path + archive_save_subpath, shell=True)
            subprocess.check_call('cp ' + tmp_save_path + uniprot_mutation.model_filename_wt +
                                    ' ' + archive_save_path + uniprot_mutation.model_filename_wt, shell=True)
            subprocess.check_call('cp ' + tmp_save_path + uniprot_mutation.model_filename_mut +
                                    ' ' + archive_save_path + uniprot_mutation.model_filename_mut, shell=True)
            with open(archive_save_path + archive_save_subpath + 'mutation.json', 'w') as fh:
                json.dump(row2dict(uniprot_mutation), fh, indent=4, separators=(',', ': '))
        
        if not self.is_immutable:
            self.session.add(uniprot_mutation)
            self.session.commit()
 
    
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
    
    
    def close(self):
        self.session.commit()
        self.session.close()

    ###########################################################################


if __name__ == '__main__':
#    return
    # run to generate an initial state database (with no precalculatios)
    print 0/0
    print sql_flavor
    db = MyDatabase('/home/kimlab1/strokach/working/pipeline/db/', 
                    sql_flavor=sql_flavor,
                    clear_database=True)
    db.load_db_from_csv()
    db.session.close()

###############################################################################

#
#class UniprotSequence(object):
#    
#    def __init__(self, database):
#        self.database = database
#        
#        self.database_name = 'pipeline_uniprot_sequences.pickle'
#        
#        if isfile(self.database_name):
#            print 'Loading uniprot sequence database'
#            f = open(self.database_name, 'rb')    
#            self.uniprot_data = pickle.load(f)
#            print 'it contains', len(self.uniprot_data), 'sequences'
#            f.close()
#        else:
#            self.uniprot_data = dict()
#
#
#            
##        #######################################################
##        ### one time: add sequences from a files            ###
##        ### do that if you need to expand the database      ###
##        ### fetching new sequences does not work on scient! ###
##
##        fileNames = ['/home/niklas/playground/mutations_clasified_recep/uniprot_list_interactions.fasta', \
##                     '/home/niklas/playground/mutations_clasified_recep/uniprot_list.fasta' ]
##        fileNames = ['/home/niklas/playground/mutations_clasified_recep/hapmap_uniprot_sequences.fasta', ]
##        for fileName in fileNames:
##            for seq_record in SeqIO.parse(fileName, "fasta"):
##                seq_record.id = seq_record.id.split('|')[1]
##                self.uniprot_data[seq_record.id] = seq_record
##        # save the newly added sequences
##        self.close()
##
##        ###           end                                   ###
##        #######################################################
#
#        
#    def __call__(self, uniprot_id):
#        """
#        returns the uniprot sequence. If the sequence is not in the database it
#        tries to retrieve it from the uniprot website.
#        
#        Note: retrieval from the website does not work on Scinet!
#        
#        """
#        if uniprot_id in ['A6NF79', 'C9JUS1', 'Q6N045', 'A6NMD1']:
#            # these uniprotKBs made problems
#            return 'no sequences'
#        elif uniprot_id == 'P02735':
#            # this sequence got replaced. I don't know right now if I can take
#            # replaced sequence so I rather dismiss it.
#            return 'no sequences'
#        
#        # the True/False value is used to add new sequences to the database in
#        # end. Only usefull if you run one instance at a time otherwise you will
#        # get diverging sequence databases.
#        try:
#            return self.uniprot_data[uniprot_id]
#        except KeyError:
#            childProcess = subprocess.Popen('whoami', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#            whoami, error = childProcess.communicate()
#            if whoami.strip() == 'joan':
#                print 'uniprot sequence not found'
#                return None
#            else:
#                print 'Fetching uniprot sequence', uniprot_id, 'from server'
#                address = 'http://www.uniprot.org/uniprot/' + uniprot_id + '.fasta'
#                handle = urllib2.urlopen(address)
#                sequence = next(SeqIO.parse(handle, "fasta"))
#                sequence.id = uniprot_id
#                self.uniprot_data[uniprot_id] = sequence
#                return sequence
#
#    
#    def add(self, items):
#        """
#        add the new items (which is a list of tuples) to the database
#        """
#        for key, value in items:
#            if self.uniprot_data.has_key(key):
#                print 'Strange..., the uniprot database seems to already have the entry', key, value
#                print 'I am overwriting it'
#            print 'adding', key, 'to the uniprot database'
#            self.uniprot_data[key] = value
#
#    
#    def close(self):
#        """
#        save the database back to disk. Do it at the end if
#        new sequences where added.
#        """
#        print 'saving the uniprot database'
#        f = open(self.database_name, 'wb')    
#        pickle.dump(self.uniprot_data, f)
#        f.close()
#        


#class PDBResolution(object):
#    """
#    In the file pdbtosp.txt from the pdb website one finds the measurement type
#    and resolution of each structure. This information is used for the selection
#    of the best interface template.
#    """
#    def __init__(self, pdbtosp):
#        self.pdbResolution_xray = dict()
#        self.pdbResolution_nmr = dict()
#        self.pdbResolution_rest = dict()
#        
#        with open(pdbtosp, 'r') as f:
#            for x in range(25):
#                f.readline()
#            for l in f:
#                if l.strip() == '':
#                    break
#                if l[0] == ' ':
#                    continue
#                line = [ item.strip() for item in l.split(' ') if item != '']
#                try:
#                    if line[1] == 'X-ray' or line[1] == 'Neutron':
#                        self.pdbResolution_xray[line[0]] = float(line[2])
#                    elif line[1] == 'NMR':
#                        self.pdbResolution_nmr[line[0]] = float(line[2])
#                    elif line[1] in ['Model', 'Other', 'IR']:
#                        continue
#                    elif line[1] in ['EM', 'Fiber']:
#                        self.pdbResolution_rest[line[0]] = float(line[2])
#                    else:
#                        print 'Could not associate the method for', line[1]
#                except:
#                    continue
#    
#    def __call__(self, pdbID):
#        """
#        the first return value is used to indicate the measurement type,
#        the second return value is the resolution.
#        
#        0 means X-ray, 2 NMR, 3 other
#        """
#        if self.pdbResolution_xray.has_key(pdbID):
#            return 0, self.pdbResolution_xray[pdbID]
#        elif self.pdbResolution_nmr.has_key(pdbID):
#            return 1, self.pdbResolution_nmr[pdbID]
#        elif self.pdbResolution_rest.has_key(pdbID):
#            return 2, self.pdbResolution_rest[pdbID]
#        else:
#            return 'None', '-'

