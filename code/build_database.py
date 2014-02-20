# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 22:02:14 2014

@author: alexey
"""
import os
import subprocess
import json 
import argparse 

import pandas as pd
from sqlalchemy import create_engine

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Alphabet import IUPAC

import class_sql as sql
import class_pdbTemplate



sql_flavor = sql.sql_flavor
path_to_sqlite_db = '/home/kimlab1/database_data/elaspic/temp.db'
clear_database = False
Base = sql.Base
Session = sql.Session


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
    engine = create_engine('postgresql://elaspic:elaspic@192.168.6.19:5432/elaspic')

if clear_database:
    Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)


Session.configure(bind=engine, autocommit=autocommit, autoflush=autoflush)
session = Session()






class CD_HIT(object):
    
    def __init__(self, bin_path, pdb_path, output_path):
        self.bin_path = bin_path
        self.pdb_path = pdb_path
        self.output_path = output_path
        if not os.path.isdir(output_path):
            print 'mkdir -p %s' % output_path
            subprocess.check_call('mkdir -p %s' % output_path, shell=True)
            subprocess.check_call('cp %s %s' % (self.bin_path + 'cd-hit', output_path), shell=True)




    def run(self, domains, identity_cutoff=0.9):
               
        seqrecords = []
        for domain in domains:
            try:
                pdb_structure = class_pdbTemplate.get_PDB(domain.pdb_id, self.pdb_path)
                model = pdb_structure[0]
                chain = model[domain.pdb_chain]
                # setting aa_only to False Selenomethionines are reported in the
                # sequence as well
                # see: http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html
                sequence = Seq('', IUPAC.protein)
                for pb in PPBuilder().build_peptides(chain, aa_only=False):
                    sequence += pb.get_sequence()
                if str(sequence) == '':
                    continue
                seqrecord = SeqRecord(sequence, id=domain.cath_id, name=domain.cath_id, description=domain.pdb_id + '\t' + domain.pdb_chain + '\t' + domain.pfam_name)
            except Exception as e:
                with open(domain.cath_id + '_error.log', 'w') as fh:
                    fh.writelines(e)
                continue
            seqrecords.append(seqrecord)
        
        SeqIO.write(seqrecords, self.output_path + domain.pfam_name + '.fasta', 'fasta')
        
        system_command = 'cd-hit -i %s -o %s -c %f -M 12000 -T 7' \
            % (self.output_path + domain.pfam_name + '.fasta', self.output_path + domain.pfam_name + '.clusters.fasta ', identity_cutoff)
        childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, error = childProcess.communicate()
        if childProcess.returncode != 0:
            raise Exception(error)
        
        with open(self.output_path + domain.pfam_name + '.clusters.fasta.clstr', 'r') as fh:
            domain_cluster_assignments = dict()
            have_cluster = False
            for line in fh:
                if line[0] == '>':
                    if have_cluster:
                        cluster_data.sort(key=lambda x: x[1], reverse=True)
                        for cluster_idx, (cath_id, seq_length, seq_identity) in enumerate(cluster_data):
                            domain_cluster_assignments[cath_id] = (cluster_id, cluster_idx, seq_length, seq_identity)
                    cluster_id = int(line.split()[-1])
                    cluster_data = []
                    have_cluster = True
                    continue
                row = line.split()
                seq_length, cath_id = int(row[1][:-3]), row[2][1:-3]
                if row[-1] == '*':
                    seq_identity = None
                else:
                    seq_identity = float(row[-1][:-1])
                cluster_data.append([cath_id, seq_length, seq_identity])
            if have_cluster:
                cluster_data.sort(key=lambda x: x[1], reverse=True)
                for cluster_idx, (cath_id, seq_length, seq_identity) in enumerate(cluster_data):
                    domain_cluster_assignments[cath_id] = (cluster_id, cluster_idx, seq_length, seq_identity)
        return domain_cluster_assignments


def load_db_from_csv():
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
        session.add(sql.Domain(**row.to_dict()))
        if idx % 10000 == 0:
#                self.session.flush()
            print idx
    session.commit()
    print 'Finished populating table domain'
    
    
    # Table `domain_contact`
    names = ['domain_contact_id', 'cath_id_1', 'pdb_contact_residues_1', 'cath_id_2', 'pdb_contact_residues_2']
    domain_contact_df = pd.read_csv(domain_contact_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
    domain_contact_df = domain_contact_df.dropna() # only a couple of rows are droppeds
    for idx, row in domain_contact_df.iterrows():
        session.add(sql.DomainContact(**row.to_dict()))
        if idx % 10000 == 0:
#                self.session.flush()
            print idx
    session.commit()
    print 'Finished populating table domain_contact'
    
    
    # Table `uniprot_sequence`
    names = ['uniprot_id', 'uniprot_name', 'uniprot_description', 'uniprot_sequence']
    uniprot_sequence_df = pd.read_csv(uniprot_sequence_infile, sep='\t', quoting=1, na_values='\N', names=names, header=None)
    for idx, row in uniprot_sequence_df.iterrows():
        session.add(sql.UniprotSequence(**row.to_dict()))
        if idx % 10000 == 0:
#                self.session.flush()
            print idx
    session.commit()
    print 'Finished populating table uniprot_sequence'     
    
    
    # Table `uniprot_domain`
    if os.path.isfile(uniprot_domain_infile):
#            names = ['uniprot_domain_id', 'uniprot_id', 'pfam_name', 'alignment_def']
        uniprot_domain_df_with_id = pd.read_csv(uniprot_domain_infile, sep='\t', na_values='\N', index_col=False)
        uniprot_domain_df_with_id['alignment_defs'] = uniprot_domain_df_with_id['alignment_def']
        uniprot_domain_df_with_id['alignment_def'] = uniprot_domain_df_with_id['alignment_defs'].apply(lambda x: sql.encode_domain(sql.decode_domain(x)))
        
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
            session.add(sql.UniprotDomain(**row.to_dict()))
            if idx % 10000 == 0:
#                    self.session.flush()
                print idx
        session.commit()
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
        
        temp = uniprot_domain_pair_df_with_id.merge(uniprot_domain_df_with_id, how='left', left_on='uniprot_domain_id_1', right_on='uniprot_domain_id').merge(uniprot_domain_df_with_id, how='left', left_on='uniprot_domain_id_2', right_on='uniprot_domain_id', suffixes=('_1', '_2'))
        temp['path_to_data'] = (temp['path_to_data_1'] + temp['pfam_name_2'] + '*' + temp['alignment_def_2'].apply(lambda x: x.replace(':','-')) + '/' + temp['uniprot_id_2'] + '/')
        temp = temp[['uniprot_domain_pair_id', 'path_to_data']]
        uniprot_domain_pair_df_with_id = uniprot_domain_pair_df_with_id.merge(temp, how='left')

#            uniprot_domain_pair_df_with_id.to_sql('uniprot_domain_pair', conn, flavor=sql_flavor, if_exists='append')
        for idx, row in uniprot_domain_pair_df_with_id.iterrows():
            session.add(sql.UniprotDomainPair(**row.to_dict()))
            if idx % 10000 == 0:
#                    self.session.flush()
                print idx
        session.commit()
        print 'Finished populating table domain'  
    else:
        pass
#            pfam_parser = parse_pfamscan.make_uniprot_domain_pair_database(domain_df, domain_contact_df, uniprot_domain_df_with_id,
#                                infile='/home/kimlab1/strokach/working/databases/biogrid/pairs_of_interacting_uniprots_human.tsv')
#            uniprot_domain_pair_df = pfam_parser.get_dataframe()
##            uniprot_domain_pair_df.to_sql('uniprot_domain_pair', conn, flavor=sql_flavor, if_exists='append')
##            uniprot_domain_pair_df_with_id = pd.read_sql('SELECT * from uniprot_domain_pair', conn)
#            uniprot_domain_pair_df_with_id.to_csv(uniprot_domain_pair_infile, sep='\t', na_values='\N', index=False)


def load_db_from_archive():
    """
    """
    data = [
        ['human/*/*/*/*/template.json', sql.UniprotDomainTemplate],
        ['human/*/*/*/*/model.json', sql.UniprotDomainModel],
        ['human/*/*/*/*/*/mutation.json', sql.UniprotDomainMutation],
        ['human/*/*/*/*/*/*/template.json', sql.UniprotDomainPairTemplate],
        ['human/*/*/*/*/*/*/model.json', sql.UniprotDomainPairModel],
        ['human/*/*/*/*/*/*/*/mutation.json', sql.UniprotDomainPairMutation],
    ]
    
    for d in data:
        command = 'ls ' + path_to_archive + d[0]
        childProcess = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, __ = childProcess.communicate() 
        filenames = [fname for fname in result.split('\n') if fname != '']
        for filename in filenames:
            with open(filename, 'r') as fh:
                row = json.load(fh)
            try:
                session.merge(d[1](**row))
            except TypeError as e:
                print 'Error merging %s.\nProbably from an older version of the database. Skipping...' % filename
                print '\t', e
            print 'Merged %s' % filename
        session.commit()
        print 'Committed changes\n\n\n'
    



if __name__ == '__main__':
            
    bin_path = '/home/kimlab1/strokach/working/pipeline/bin/'
    pdb_path = '/home/kimlab1/database_data/pdb/data/data/structures/divided/pdb/'
    output_path = '/tmp/run_cdhit/'
    db = sql.MyDatabase()
    cd_hit = CD_HIT(bin_path, pdb_path, output_path)
    
    # read which configFile to use    
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pfam_name')
    arguments = parser.parse_args()

    domains = []
    for domain in db.get_domain(arguments.pfam_name):
        if not domain.cdhit_cluster_length:
            domains.append(domain)

    if domains:
        domain_cluster_assignment = cd_hit.run(domains)
        print domain_cluster_assignment
        for domain in domains:
            domain.cdhit_cluster, domain.cdhit_cluster_idx, \
            domain.cdhit_cluster_length, domain.cdhit_cluster_identity = \
                domain_cluster_assignment.get(domain.cath_id, (-1,None,None,None))[:4]
            db.session.add(domain)
        db.session.commit()
        db.session.close()

    