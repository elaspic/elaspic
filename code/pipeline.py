# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:09:11 2013

@author: niklas
"""
import numpy as np
import pandas as pd

from sets import ImmutableSet

from ConfigParser import SafeConfigParser
import subprocess
import os
import sys

import multiprocessing
import optparse
from string import uppercase

import logging

from Bio.SubsMat import MatrixInfo
from Bio import SeqIO

from class_multi import Consumer
from class_multi import Task
from class_error import DataError, ConfigError
from class_logging import MultiProcessingLog

from scinetCleanup import scinetCleanup
from os.path import isfile
import cPickle as pickle
import urllib2

import parse_pfamscan

from collections import OrderedDict


class NullDevice():
    def write(self, s):
        pass

# credit goes to here:
# http://pymotw.com/2/multiprocessing/communication.html#controlling-access-to-resources
class ActivePool(object):
    """
    Used to control how many parallel T_Coffee calls can be run
    Originally implemented because I had problems with running T_Coffee in
    parallel it can now be used if there are memory limitations
    
    The problem I had with T_Coffee is that I made the programm call from the
    same directory. Switching to a unique directory for every call solved
    the issue.
    """
    def __init__(self):
        super(ActivePool, self).__init__()
        self.mgr = multiprocessing.Manager()
        self.active = self.mgr.list()
        self.lock = multiprocessing.Lock()
    def makeActive(self, name):
        with self.lock:
            self.active.append(name)
    def makeInactive(self, name):
        with self.lock:
            self.active.remove(name)
    def __str__(self):
        with self.lock:
            return str(self.active)



class myDatabases(object):
    
    def __init__(self):
        pass
    
    def split_domain(self, domain):
        """ 
        Takes a string of two domain boundaries and returns a list with int
        The separator is '-' and it can happen that both or one boundary is
        negative, i.e.
        
            -150-200,   meaning from -150 to 200
            -150--100,  meaning from -150 to -100
            etc.
        
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
        
        
    def split_domain_semicolon(self, domains):
        """ Unlike split_domain(), this function returns a tuple of tuples of strings,
        preserving letter numbering (e.g. 10B)
        """
        x = domains
        return tuple([ tuple([r.strip() for r in ro.split(':')]) for ro in x.split(',') ])

        
    def split_interface_aa(self, interface_aa):
        """
        """
        if (interface_aa is not np.nan) and (interface_aa != '') and (interface_aa != 'NULL'):
            if interface_aa[-1] == ',':
                interface_aa = interface_aa[:-1]
        
            x  = interface_aa
            return_tuple = tuple([int(r.strip()) for r in x.split(',')])
            
        else:
            return_tuple = np.nan
            
        return return_tuple



class UniprotSequence(object):
    
    def __init__(self, database):
        self.database = database
        
        self.database_name = 'pipeline_uniprot_sequences.pickle'
        
        if isfile(self.database_name):
            print 'Loading uniprot sequence database'
            f = open(self.database_name, 'rb')    
            self.uniprot_data = pickle.load(f)
            print 'it contains', len(self.uniprot_data), 'sequences'
            f.close()
        else:
            self.uniprot_data = dict()
            
#        #######################################################
#        ### one time: add sequences from a files            ###
#        ### do that if you need to expand the database      ###
#        ### fetching new sequences does not work on scient! ###
#
#        fileNames = ['/home/niklas/playground/mutations_clasified_recep/uniprot_list_interactions.fasta', \
#                     '/home/niklas/playground/mutations_clasified_recep/uniprot_list.fasta' ]
#        fileNames = ['/home/niklas/playground/mutations_clasified_recep/hapmap_uniprot_sequences.fasta', ]
#        for fileName in fileNames:
#            for seq_record in SeqIO.parse(fileName, "fasta"):
#                seq_record.id = seq_record.id.split('|')[1]
#                self.uniprot_data[seq_record.id] = seq_record
#        # save the newly added sequences
#        self.close()
#
#        ###           end                                   ###
#        #######################################################
        
    def __call__(self, uniprotKB):
        """
        returns the uniprot sequence. If the sequence is not in the database it
        tries to retrieve it from the uniprot website.
        
        Note: retrieval from the website does not work on Scinet!
        
        """
        if uniprotKB in ['A6NF79', 'C9JUS1', 'Q6N045', 'A6NMD1']:
            # these uniprotKBs made problems
            return 'no sequences'
        elif uniprotKB == 'P02735':
            # this sequence got replaced. I don't know right now if I can take
            # replaced sequence so I rather dismiss it.
            return 'no sequences'
        
        # the True/False value is used to add new sequences to the database in
        # end. Only usefull if you run one instance at a time otherwise you will
        # get diverging sequence databases.
        try:
            return self.uniprot_data[uniprotKB], False
        except KeyError:
            childProcess = subprocess.Popen('whoami', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            whoami, error = childProcess.communicate()
            if whoami.strip() == 'joan':
                print 'uniprot sequence not found'
                return 'no sequences'
            else:
                print 'Fetching uniprot sequence', uniprotKB, 'from server'
                address = 'http://www.uniprot.org/uniprot/' + uniprotKB + '.fasta'
                handle = urllib2.urlopen(address)
                sequence = next(SeqIO.parse(handle, "fasta"))
                sequence.id = uniprotKB
                self.uniprot_data[uniprotKB] = sequence
                return sequence, True

    
    def add(self, items):
        """
        add the new items (which is a list of tuples) to the database
        """
        for key, value in items:
            if self.uniprot_data.has_key(key):
                print 'Strange..., the uniprot database seems to already have the entry', key, value
                print 'I am overwriting it'
            print 'adding', key, 'to the uniprot database'
            self.uniprot_data[key] = value

    
    def close(self):
        """
        save the database back to disk. Do it at the end if
        new sequences where added.
        """
        print 'saving the uniprot database'
        f = open(self.database_name, 'wb')    
        pickle.dump(self.uniprot_data, f)
        f.close()
        
        

class DomainDefinition(myDatabases):
    """ Contains pdbfam-based definitions of all pfam domains in the pdb
    """
    
    def __init__(self, textfile_name):
        """
        """
        db_name = 'domain_definition_database.pickle'
        
        if isfile(db_name):
            # A precompiled database exists
            print 'Loading domain definition database...'
            
            self.db = pd.read_pickle(db_name)
            self.grouped_db = self.db.groupby('pfam_name') 

        else:        
            # No precompiled database exists
            print 'Generating domain definition database...'
            
            names = ['pdb_id', 'autoPFAMA', 'pfam_name', 'pdb_chain', 'pdb_definition', 'cath_id']
            
            db = pd.read_csv(textfile_name, sep='\t', header=False, names=names)
            
            db['pdb_definition'] = db['pdb_definition'].apply(self.split_domain_semicolon)
            
            self.db = db
            self.db.to_pickle(db_name)
            self.grouped_db = db.groupby('pfam_name')
                    
                    
    def __call__(self, pfam_name):
        if self.grouped_db.groups.has_key(pfam_name):
            return self.grouped_db.get_group(pfam_name)
        else:
            print "No entries for the given pfam_name:", pfam_name
            return None



class DomainInteraction(myDatabases):
    """ Keeps the domain-domain interaction information from pdbfam
    """
    
    def __init__(self, textfile_name):
        """ Create the internal dictionary using a tsv file exported from mysql
        """
        db_name = 'domain_interaction_database.pickle'
        
        # Check if a pickled database exists
        if isfile(db_name):
            # A precompiled database exists
            print 'Loading domain interaction database...'
            
            self.db = pd.read_pickle(db_name)
            self.grouped_db = self.db.groupby(['pfam_name_1','pfam_name_2'])
            
        else:        
            # No precompiled database exists
            print 'Generating domain interaction database...'

            names = [
            'pdb_id', 'pfam_name_1', 'pdb_chain_1', 'pdb_definition_1', 'pdb_interface_aa_1', 'cath_id_1', 
            'pfam_name_2', 'pdb_chain_2', 'pdb_definition_2', 'pdb_interface_aa_2', 'cath_id_2', 'domain_contact_id',]
            
            db = pd.read_csv(textfile_name, sep='\t', header=False, names=names)
            
            db['pdb_definition_1'] = db['pdb_definition_1'].apply(self.split_domain_semicolon)
            db['pdb_definition_2'] = db['pdb_definition_2'].apply(self.split_domain_semicolon)
            
            db['pdb_interface_aa_1'] = db['pdb_interface_aa_1'].apply(self.split_interface_aa)
            db['pdb_interface_aa_2'] = db['pdb_interface_aa_2'].apply(self.split_interface_aa)
            
            self.db = db
            self.grouped_db = self.db.groupby(['pfam_name_1', 'pfam_name_2'])
            
            self.db.to_pickle(db_name)


    def __call__(self, pfam_name_1, pfam_name_2):
        """
        Note that the produced dataframe may not have the same order as the keys
        """
        key1 = tuple([pfam_name_1, pfam_name_2])
        key2 = tuple([pfam_name_2, pfam_name_1])
            
        if self.grouped_db.groups.has_key(key1):
            return_df_1 = self.grouped_db.get_group(key1)
        else:
            return_df_1 = None
            
        if self.grouped_db.groups.has_key(key2):
            return_df_2 = self.grouped_db.get_group(key2)
        else:
            return_df_2 = None
            
        if not return_df_1 and not return_df_2:
            print "No entries for the given domain combination:", pfam_name_1, pfam_name_2
            
        return [return_df_1, return_df_2]



class PDBResolution(object):
    """
    In the file pdbtosp.txt from the pdb website one finds the measurement type
    and resolution of each structure. This information is used for the selection
    of the best interface template.
    """
    def __init__(self, pdbtosp):
        self.pdbResolution_xray = dict()
        self.pdbResolution_nmr = dict()
        self.pdbResolution_rest = dict()
        
        with open(pdbtosp, 'r') as f:
            for x in range(25):
                f.readline()
            for l in f:
                if l.strip() == '':
                    break
                if l[0] == ' ':
                    continue
                line = [ item.strip() for item in l.split(' ') if item != '']
                try:
                    if line[1] == 'X-ray' or line[1] == 'Neutron':
                        self.pdbResolution_xray[line[0]] = float(line[2])
                    elif line[1] == 'NMR':
                        self.pdbResolution_nmr[line[0]] = float(line[2])
                    elif line[1] in ['Model', 'Other', 'IR']:
                        continue
                    elif line[1] in ['EM', 'Fiber']:
                        self.pdbResolution_rest[line[0]] = float(line[2])
                    else:
                        print 'Could not associate the method for', line[1]
                except:
                    continue
    
    def __call__(self, pdbID):
        """
        the first return value is used to indicate the measurement type,
        the second return value is the resolution.
        
        0 means X-ray, 2 NMR, 3 other
        """
        if self.pdbResolution_xray.has_key(pdbID):
            return 0, self.pdbResolution_xray[pdbID]
        elif self.pdbResolution_nmr.has_key(pdbID):
            return 1, self.pdbResolution_nmr[pdbID]
        elif self.pdbResolution_rest.has_key(pdbID):
            return 2, self.pdbResolution_rest[pdbID]
        else:
            return 'None', '-'



class ProteinDefinition(myDatabases):
    """ Initiated using parsed pfam_scan.pl output for the human uniprot (or the entire uniprot)
    The database is updated with information about calculated models
    """
    
    def __init__(self):
        
        # Check if a pickled database exists
        db_name = 'protein_definition_database.pickle'
        
        if isfile(db_name):
            # A precompiled database exists
            print 'Loading protein definition database...'
            self.db = pd.read_pickle(db_name)
            self.grouped_db = self.db.groupby('uniprot_id')
        
        else:
            # Create a database from a file containing fasta sequences, the pfam_scan.pl output from that file
            # Also needs a list of clans, a list os repeating domains and a list of superdomains
            # Filenames are hardcoded at this points
            print 'Generating protein definition database...'

            pfam_parser = parse_pfamscan.make_uniprot_pfam_database()
            pfam_parser.run()
            db = pfam_parser.get_dataframe()
            
            self.db = db
            print self.db
            self.db.to_pickle(db_name)
            self.grouped_db = self.db.groupby('uniprot_id')


    def __call__(self, uniprot_id):
        """ 
        """
        # using (uniprot_id, uniprot_id,) because in the future we might want to incorporate
        # splicing and other variants, which will be identified by the second uniprot_id
        if self.grouped_db.groups.has_key(uniprot_id):
            return self.grouped_db.get_group(uniprot_id)
        else:
            # In the future, make it fetch the uniprot sequence, analyse it with pfam_scan.pl,
            # and save the output to the database
            print "No protein definition entries for the given uniprot_id:", uniprot_id
            return None
            
            
    def update(self, modified_df):
        """ Updates the internal data frame with the modifications present in the
        copied dataframe
        """
        self.db.update(modified_df)
        


class ProteinInteraction(object):
    """ Contains known interactions between uniprot proteins
    """
    
    def __init__(self, textfile_name, ProteinDefinition, DomainInteraction):
        """ 
        Checks if the interaction database is already available as a pickled object.
        Generates it from file 'textfile_name' otherwise.
        """
        db_name = 'protein_interaction_database.pickle'
        
        if isfile(db_name):
            # A precompiled database exists
            print 'Loading protein interaction database...'
            
            self.db = pd.read_pickle(db_name)
            self.grouped_db_1 = self.db.groupby(['uniprot_id_1'])
            self.grouped_db_2 = self.db.groupby(['uniprot_id_2'])
        
        else:            
            # No precompiled database exists
            print 'Generating protein interaction database...'
            
            # Read in a list of unteracting uniprot pairs obtained using BioGrid
            # (and biomart to convert enesmble gene accession to uniprot)
            set_of_interactions = set()
            with open(textfile_name, 'r') as fh:
                next(fh)
                for line in fh:
                    row = [l.strip() for l in line.split('\t')]
                    set_of_interactions.add(ImmutableSet([row[0], row[1]]))
            
            
            # For every unique uniprot pair, try to find a template in the 
            # DomainInteraction database (which is based on pfam...)
            
            # Turn off print statement output (otherwise too many domains not found)
            original_stdout = sys.stdout
            sys.stdout = NullDevice()
            number_of_missing_pfams = 0 # for debugging
            
            # Lists to keep the obtained information
            list_of_uniprot_ids_1 = []
            list_of_uniprot_ids_2 = []
            list_of_contact_ids_forward = []
            list_of_contact_ids_reversed = []
            
            for idx, interaction in enumerate(set_of_interactions):
                
                interaction = list(interaction)
                # Get a list of pfam domains for each protein in a pair
                if len(interaction) == 1:
                    uniprot_id_1 = interaction[0]
                    uniprot_id_2 = uniprot_id_1
                    
                    pfam_names_1 = ProteinDefinition(uniprot_id_1)
                    if pfam_names_1:
                        pfam_names_1 = list(pfam_names_1['seq_id'])
                        pfam_names_2 = pfam_names_1
                    else:
                        number_of_missing_pfams += 1
                        continue
                elif len(interaction) == 2:
                    uniprot_id_1 = interaction[0]
                    uniprot_id_2 = interaction[1]
                    
                    pfam_names_1 = ProteinDefinition(uniprot_id_1)
                    if pfam_names_1:
                        pfam_names_1 = list(pfam_names_1['seq_id'])
                    else:
                        number_of_missing_pfams += 1
                        continue
                    
                    pfam_names_2 = ProteinDefinition(uniprot_id_2)
                    if pfam_names_2:
                        pfam_names_2 = list(pfam_names_2['seq_id'])
                    else:
                        number_of_missing_pfams += 1
                        continue
                else:
                    raise Exception
                
                # Go over each possible pair of pfam domains and see if we can 
                # find a pdb template
                for pfam_name_1 in pfam_names_1:
                    for pfam_name_2 in pfam_names_2:
                        domain_interaction_forward, domain_interaction_reverse = DomainInteraction(pfam_name_1,pfam_name_2)
                        
                        contact_ids_forward = tuple()
                        if domain_interaction_forward:
                            contact_ids_forward = tuple(domain_interaction_forward['domain_contact_id'])
                        
                        contact_ids_reverse = tuple()
                        if domain_interaction_reverse:
                            contact_ids_reverse = tuple(domain_interaction_reverse['domain_contact_id'])
                
                        if len(contact_ids_forward) > 0 or len(contact_ids_reverse) > 0:
                            list_of_uniprot_ids_1.append(uniprot_id_1)
                            list_of_uniprot_ids_2.append(uniprot_id_2)
                            list_of_contact_ids_forward.append(contact_ids_forward)
                            list_of_contact_ids_reversed.append(contact_ids_reverse)
            
            # Turn print statement output back on
            sys.stdout = original_stdout
            print number_of_missing_pfams
            
            # Create a dataframe with the compiled results
            db = pd.DataFrame(index=range(len(list_of_uniprot_ids_1)))
            db['uniprot_id_1'] = list_of_uniprot_ids_1
            db['uniprot_id_2'] = list_of_uniprot_ids_2
            db['contact_ids_forward'] = list_of_contact_ids_forward
            db['contact_ids_reversed'] = list_of_contact_ids_reversed
            print db
            db.drop_duplicates(inplace=True)
            print db
            
            self.db = db
            self.db.to_pickle(db_name)
            self.grouped_db_1 = self.db.groupby(['uniprot_id_1'])
            self.grouped_db_2 = self.db.groupby(['uniprot_id_2'])
            
    
    def __call__(self, uniprot_id):
        """
        """
        if self.grouped_db_1.groups.has_key(uniprot_id):
            result_df_1 = self.grouped_db_1.get_group(uniprot_id)
        else:
            result_df_1 = None
            
        if self.grouped_db_2.groups.has_key(uniprot_id):
            result_df_1 = self.grouped_db_2.get_group(uniprot_id)
        else:
            result_df_2 = None
            
        if not result_df_1 and not result_df_2:
            print "No entries for the given uniprot_id:", uniprot_id
           
        return [result_df_1, result_df_2]
    
    
    def upldate(self, modified_df):
        """ Add new interaction models to the database
        """
        pass



#class get3DID(myDatabases):
#    """
#    To improve the 3DID file the domain boundaries are corrected by
#    the improved boundaries generated by Sebastian
#    """
#    def __init__(self, threeDidFile, domainTableFile):
#        # I hardcoded the filenames
#        database_name = 'pipeline_3DID_database.pickle'
#        
#        # check if the database was already created
#        if isfile(database_name):
#            # load the database
#            print 'Loading the corrected 3DID table'
#            f = open(database_name, 'rb')
#            self.database = pickle.load(f)
#            f.close()
#            
#        else:            
#            print 'Generating the correceted 3DID table'
#            # read the 3Did file and store the information in a dict()
#            self.database_3did = self.__make_3did_table(threeDidFile)
#            
#            # read the corrected domain boundaries and store them in a dict()
#            self.domainTable = self.__make_domain_boundary_table(domainTableFile)
#            
#            # create the database
#            self.database = self.__correct_3did_domain_boundaries()
#            
#            # save the database
#            f = open(database_name, 'wb')
#            pickle.dump(self.database, f)
#            f.close()
#          
#    
#    def __call__(self, PfamID_firstGuy, PfamID_secondGuy):
#        if self.database.has_key((PfamID_firstGuy, PfamID_secondGuy)):
#            return self.database[PfamID_firstGuy, PfamID_secondGuy]
#        else:
#            return 'no entry'
#
#    
#    def __make_3did_table(self, threeDidFile):
#        """
#        read the 3DID file and store the information in a dictionary
#        """
#        database = dict()
#        with open(threeDidFile, 'r') as f:
#            for line in f:
#                # go to the line corresponding to one family pair
#                if line[:4] == '#=ID':
#                    # get the two family ids
#                    Pfam1, Pfam2 = line.split('\t')[1], line.split('\t')[2]
#                
#                # read the entries for the above set pfam family
#                if line[:4] == '#=3D':                    
#                    pdb = line.split('\t')[1]
#                    chain1 = line.split('\t')[2]
#                    chain2 = line.split('\t')[3]
#                    
#                    # there often are interactions with one chain
#                    # we are only interested in interactions between different chains
#                    if chain1.split(':')[0] !=  chain2.split(':')[0]:
#                        chain1, domain_pdb1 = chain1.split(':')
#                        chain2, domain_pdb2 = chain2.split(':')
#                        add = [pdb.upper(), chain1, chain2, domain_pdb1, domain_pdb2]
#                        database.setdefault((Pfam1, Pfam2), []).append(add)
#
#        return database
#
#    
#    def __make_domain_boundary_table(self, domainTableFile):
#        """
#        creates the domain boundary table Sebastian provided to correct the
#        domain boundaries of the 3did_flat file
#        """
#        domainTable = dict()
#        with open(domainTableFile, 'r') as f:
#            f.readline()
#            for line in f:
#                pdb, AutoPfamA, pfam_id, chain, domain_boundary = line.split('\t')
#                # it happend that Sebastian was "glueing" two Pfam domain together
#                # forming kind of a super domain (to get a "real", i.e. meaningful,
#                # domain). In that case he added them with '+'. I split them and
#                # add them in a list. Later it is checked wether the pfam ID is in
#                # this list, i.e. if it is either one or the other
#                pfam_id = pfam_id.strip('+').split('+')
#                domain_boundary = domain_boundary.strip().split(',')
#                for item in domain_boundary:
#                    # it can happen that the pdb numbering is negativ
#                    # see split_domain() function explanation
#                    domain = self.split_domain(item)
#
#                domainTable.setdefault(pdb + chain, []).append([pfam_id, domain])
#        
#        return domainTable
#
#    
#    def __correct_3did_domain_boundaries(self):
#        """
#        The 3DID file contains all pdb files that have one specific interaction
#        given a family pair. Sebastian improved the domain boundaries for the
#        whole pdb and his corrected boundaries are implemented here.
#        """
#        database = dict()
#        
#        for key in self.database_3did:
#            Pfam1_ID, Pfam2_ID = key
#            for interaction in self.database_3did[key]:
#                pdb, chain1, chain2, domain1_pdb, domain2_pdb = interaction
#                
#                domain1_pdb = self.split_domain(domain1_pdb)
#                domain2_pdb = self.split_domain(domain2_pdb)
#                
#                domain1 = self.__correct_domain_boundary(Pfam1_ID, pdb, chain1, domain1_pdb)
#                domain2 = self.__correct_domain_boundary(Pfam2_ID, pdb, chain2, domain2_pdb)
#                
#                add = [pdb, chain1, chain2, domain1, domain2]
#                database.setdefault(key, []).append(add)
#                
#        return database
#                
#    
#    def __correct_domain_boundary(self, Pfam_ID, pdb, chain, domain_pdb):
#        """
#        correct the domain boundary for chain if there is information
#        for it in the database with the corrected domain information
#        
#        You can have more than one domain of the same family type in
#        one chain so one has to
#        make sure that the domain boundaries one wants to correct belong
#        to the domain that is checked. This is done by calculating the
#        overlap between the two domains.
#        see http://stackoverflow.com/a/5095171
#        to see how the overlap is calculated
#        """
#        is_weird, domain_pdb_multiset = self.__check_weird_domain_numbering(pdb, chain)
#        if not is_weird:
#            domain_pdb_multiset = collections.Counter(range(domain_pdb[0],domain_pdb[1]))
#        else:
#            domain_pdb_multiset = collections.Counter(domain_pdb_multiset)
#                    
#        set_new = False
#        if self.domainTable.has_key(pdb + chain):
#            # there might be more than one domain (of the same family), thus use a loop
#            for item in self.domainTable[pdb + chain]:
#                pfam_id, domain_new = item
#                
#                # look for the correct family
#                if Pfam_ID not in pfam_id:
#                    # that means one is not looking at the correct domain
#                    # in the chain
#                    continue
#                else:
#                    # check if the domains overlap, i.e. if one is looking
#                    # at the correct domain
#                    is_weird, domain_multiset = self.__check_weird_domain_numbering(pdb, chain)
#
#                    # see __check_weird_domain_numbering() for explanation
#                    # the function is necessary to create correct multisets
#                    if not is_weird:
#                        domain_multiset = collections.Counter(range(int(domain_new[0]),int(domain_new[1])))
#                    else:
#                        domain_multiset = collections.Counter(domain_multiset)
#
#                    # get the overlap of the two domains
#                    overlap = list( (domain_pdb_multiset & domain_multiset).elements() )
#                    # if they are not overlapping, check the next domain
#                    # in the chain
#                    if len(overlap) <= 2:
#                        continue
#                    else:
#                        domain = domain_new
#                        set_new = True
#                        
#        # if the domain was not corrected, set it to the old value
#        if set_new == False:
#            domain = domain_pdb
#        
#        return domain
#    
#    
#    def __check_weird_domain_numbering(self, pdb, chain):
#        """
#        I manually filtered some weird domain numberings. They are handled here.
#        Returns True/False as first value depending if the numbering is "weird"
#        """
#        if pdb == '6INS':
#            return True, range(1,30)
#        elif pdb + chain == '2GEDB':
#            return True, range(37, 245)
#        elif pdb + chain == '3ALXB':
#            r = range(32,127)
#            r.extend([600,601,602,603,604,605,606])
#            return True, r
#        elif pdb + chain == '3ALXA':
#            r = range(32,127)
#            r.extend([606])
#            return True, r
#        elif pdb + chain == '3ERRB':
#            r = range(3266,3427)
#            r.extend([94,95,96])
#            return True, r
#        elif pdb + chain == '1JR3B':
#            r = range(2041,2070)
#            r.extend(range(70,179))
#            return True, r
#        elif pdb + chain == '1BMV':
#            r = range(3001,3183)
#            r.extend(range(2001,2191))
#            return True, r
#        else:
#            return False, []
  
           
        
#class core_Templates(myDatabases):
#    """
#    Sebastian provided the templates for the core mutations. Here his file
#    is read and the information is stored in a dict()
#    """
#    def __init__(self, database):
#        database_filename = 'pipeline_core_templates.pickle'
#        if isfile(database_filename):
#            print 'Loading the core template database'
#            f = open(database_filename, 'rb')    
#            self.template_information = pickle.load(f)
#            f.close()
#        else:
#            print 'Reading core templates from', database
#            self.template_information = self.__read_templates(database)
#            
#            f = open(database_filename, 'wb')    
#            pickle.dump(self.template_information, f)
#            f.close()
#    
#    def __call__(self, uniprotKB):
#        if self.template_information.has_key(uniprotKB):
#            return self.template_information[uniprotKB]
#        else:
#            return []
#
#        
#    def __read_templates(self, inFile):
#        """
#        create a dict containing the template information for fast lookup
#        make sure that:
#        domain_boundary_uniprot is a list
#        amino_acids_core is a list
#        """
#        template_information = dict()
#        
#        with open(inFile, 'r') as f:
#            f.readline() # skip the header
#            for line in f:
#                uniprotKB, family_name, domain_boundary_uniprot, amino_acids_core, \
#                pdb_template, pdb_chain, pdb_domain_boundary = line.split('\t')
#                
#                # if no core is given, skip
#                if amino_acids_core == 'NULL' or pdb_domain_boundary == 'NULL':
#                    continue
#                
#                # the uniprot numbering is not weird as can be the pdb numbering
#                # so for the uniprot domain boundaries one can easily split them
#                domain_boundary_uniprot = [int(domain_boundary_uniprot.split('-')[0]), int(domain_boundary_uniprot.split('-')[1])]
#                
#                # it can happen that a domain is not continues but consists of
#                # several parts. In that case the maximum is taken to be sure
#                # to cover everything. This is one point that could be improved
#                # for that one would have to split the domain into the relevant
#                # parts and make sure that only those parts are modelled!
#                # To improve this one would have to do quite a lot of work
#                # changing the whole pipeline
#                pdb_domain_boundary_list = list()
#                for item in pdb_domain_boundary.split(','):
#                    pdb_domain_boundary_list.extend(self.split_domain(item))
#                pdb_domain_boundary = [pdb_domain_boundary_list[0], pdb_domain_boundary_list[-1]]
#                
#                # convert the amino_acids_core to a list
#                amino_acids_core = [ int(i) for i in amino_acids_core.split(',') ]
#                
#                template_information.setdefault(uniprotKB, list()).append( [family_name, \
#                                                                            domain_boundary_uniprot, \
#                                                                            amino_acids_core, \
#                                                                            pdb_template, \
#                                                                            pdb_chain, \
#                                                                            pdb_domain_boundary
#                                                                            ] )
#
#        return template_information
    
#    def close(self):
#        pass
##        print 'Saving the core database'
##        f = open(self.database, 'wb')    
##        pickle.dump(self.uniprot_data, f)
##        f.close()




class pipeline():
    
    def __init__(self, configFile):
        
        ###################################################        
        # read the configuration file and set the variables
        
        configParser = SafeConfigParser(
                                defaults={'tmpPath':'/tmp/pipeline/',
                                          'HETATM':True,
                                          'outputPath':'results/',
                                          'DEBUG':False,
                                          'saveTo':'$SCRATCH/niklas-pipeline/',
                                          'saveScinet':False,
                                          'webServer': False
                                          }
                                )

        configParser.read(configFile)
        defaults = configParser.defaults()
        ### read all the values
        ## from [DEFAULT]
        # tmpPath
        if configParser.has_option('DEFAULT', 'tmpPath'):
            tmpPath = configParser.get('DEFAULT', 'tmpPath')
        else:
            tmpPath = defaults['tmpPath']
        # DEBUG
        if configParser.has_option('DEFAULT', 'DEBUG'):
            self.DEBUG = configParser.getboolean('DEFAULT', 'DEBUG')
        else:
            self.DEBUG = defaults['DEBUG']
        # HETATM
        if configParser.has_option('DEFAULT', 'HETATM'):
            self.HETATM = configParser.getboolean('DEFAULT', 'HETATM')
        else:
            self.HETATM = defaults['HETATM']
        # saveTo
        if configParser.has_option('DEFAULT', 'saveTo'):
            self.saveTo = configParser.get('DEFAULT', 'saveTo')
        else:
            self.saveTo = defaults['saveTo']
        # saveScinet
        if configParser.has_option('DEFAULT', 'saveScinet'):
            self.saveScinet = configParser.getboolean('DEFAULT', 'saveScinet')
        else:
            self.saveScinet = defaults['saveScinet']
        # 3DID or PDBFam interaction file
        if configParser.has_option('DEFAULT', 'use_pdbfam_database'):
            self.use_pdbfam_database = configParser.getboolean('DEFAULT', 'use_pdbfam_database')
        else:
            self.use_pdbfam_database = False
        
        # web-server
        if configParser.has_option('DEFAULT', 'webServer'):
            self.webServer = configParser.get('DEFAULT', 'webServer')
        else:
            self.webServer = defaults['webServer']
            
        ## from [SETTINGS]
        # name
        if configParser.has_option('SETTINGS', 'name'):
            self.name = configParser.get('SETTINGS', 'name')
        else:
            self.name = None
        # numConsumers
        if configParser.has_option('SETTINGS', 'numConsumers'):
            self.num_consumers = configParser.getint('SETTINGS', 'numConsumers')
        else:
            self.num_consumers = multiprocessing.cpu_count()
        # tcoffee_parallel_runs
        if configParser.has_option('SETTINGS', 'tcoffee_parallel_runs'):
            self.tcoffee_parallel_runs = configParser.getint('SETTINGS', 'tcoffee_parallel_runs')
        else:
            self.tcoffee_parallel_runs = multiprocessing.cpu_count()
        # pdbPath
        if configParser.has_option('SETTINGS', 'pdbPath'):
            self.pdbPath = configParser.get('SETTINGS', 'pdbPath')
        else:
            raise ConfigError('pdbPath')
        # matrix (currently hardcoded, can be changed if needed)
        # in that case also change the gap_start and gap_extend options
        # they are used in conjuntion with the matrix to determine the
        # sequence similariy (could be used to determine the sequence similarity
        # between to interfaces. Might help to improve the modelling)
#        if configParser.has_option('SETTINGS', 'matrix'):
#            self.matrix_option = configParser.get('SETTINGS', 'matrix')
        self.matrix_option = 'blosum80'
        # gap_start
#        if configParser.has_option('SETTINGS', 'gap_start'):
#            self.gap_start = configParser.getint('SETTINGS', 'gap_start')
#        else:
#            raise ConfigError('gap_start')
        self.gap_start = -16
        # gap_extend
#        if configParser.has_option('SETTINGS', 'gap_extend'):
#            self.gap_extend = configParser.getint('SETTINGS', 'gap_extend')
#        else:
#            raise ConfigError('gap_extend')
        self.gap_extend = -4
        
        
#        # crystalPDB
#        if configParser.has_option('SETTINGS', 'crystalPDB'):
#            self.crystalPDB = configParser.getboolean('SETTINGS', 'crystalPDB')
#        else:
#            raise ConfigError('crystalPDB')
        # outputPath
        if configParser.has_option('SETTINGS', 'outputPath'):
            self.outputPath = configParser.get('SETTINGS', 'outputPath')
            # adjust the savePDB folder if the outputPath changed
            if configParser.has_option('SETTINGS', 'savePDB'):
                self.savePDB = configParser.get('SETTINGS', 'savePDB')
            else:
                self.savePDB = self.outputPath + 'pdbFiles/'
        else:
            self.outputPath = defaults['outputPath']
            self.savePDB = self.outputPath + 'pdbfiles/'
        # runTime
        if configParser.has_option('SETTINGS', 'runTime'):
            self.runTime = configParser.get('SETTINGS', 'runTime')
        else:
            self.runTime = 'INFINITE'
        # bin
        if configParser.has_option('SETTINGS', 'bin'):
            self.executables = configParser.get('SETTINGS', 'bin')
        else:
            raise ConfigError('bin')
            
        ## from [INPUT]
        # file
        if configParser.has_option('INPUT', 'file'):
            self.inputFile = configParser.get('INPUT', 'file')
            if not os.path.isfile(self.inputFile):
                raise DataError(self.inputFile)
        else:
            raise ConfigError('file')
        # mutation_uniprot
        if configParser.has_option('INPUT', 'mutation_uniprot'):
            self.mutation_uniprot = configParser.getboolean('INPUT', 'mutation_uniprot')
        else:
            raise ConfigError('mutation_uniprot')
        
        ## from [MODELLER]
        # MODELLER
#        if configParser.has_option('MODELLER', 'MODELLER'):
#            self.MODELLER = configParser.getboolean('MODELLER', 'MODELLER')
#        else:
#            raise ConfigError('MODELLER')
        # modeller_runs
        if configParser.has_option('MODELLER', 'modeller_runs'):
            self.modeller_runs = configParser.getint('MODELLER', 'modeller_runs')
        else:
            raise ConfigError('modeller_runs')
        
        ## from [FOLDX]
        # WATER
        if configParser.has_option('FOLDX', 'WATER'):
            self.foldX_WATER = configParser.get('FOLDX', 'WATER')
        else:
            raise ConfigError('WATER')
        # buildModel_runs
        if configParser.has_option('FOLDX', 'buildModel_runs'):
            self.buildModel_runs = configParser.get('FOLDX', 'buildModel_runs')
        else:
            raise ConfigError('buildModel_runs')
#        # MUTATE
#        if configParser.has_option('FOLDX', 'MUTATE'):
#            self.MUTATE = configParser.getboolean('FOLDX', 'MUTATE')
#        else:
#            raise ConfigError('MUTATE')
        
        ## from [DATABASES]
        # protein_interaction_file
        if configParser.has_option('DATABASES', 'protein_interaction_file'):
            self.protein_interaction_file = configParser.get('DATABASES', 'protein_interaction_file')
        else:
            raise ConfigError('protein_interaction_file')
        # threeDidFile
        if configParser.has_option('DATABASES', 'threeDidFile'):
            self.threeDidFile = configParser.get('DATABASES', 'threeDidFile')
        else:
            raise ConfigError('threeDidFile')
        # PfamBoundaryCorrection
        if configParser.has_option('DATABASES', 'PfamBoundaryCorrection'):
            self.PfamBoundaryCorrection = configParser.get('DATABASES', 'PfamBoundaryCorrection')
        else:
            raise ConfigError('threeDidFile')
        # uniprotSequenceDatabase
        if configParser.has_option('DATABASES', 'uniprotSequenceDatabase'):
            self.uniprotSequenceDatabase = configParser.get('DATABASES', 'uniprotSequenceDatabase')
        else:
            raise ConfigError('uniprotSequenceDatabase')
        # pdb_resolution_file
        if configParser.has_option('DATABASES', 'pdb_resolution_file'):
            self.pdb_resolution_file = configParser.get('DATABASES', 'pdb_resolution_file')
        else:
            raise ConfigError('pdb_resolution_file')
        # coreTemplatesFile
        if configParser.has_option('DATABASES', 'coreTemplatesFile'):
            self.coreTemplatesFile = configParser.get('DATABASES', 'coreTemplatesFile')
        else:
            raise ConfigError('coreTemplatesFile')
        # domain_definition_file
        if configParser.has_option('DATABASES', 'domain_definition_file'):
            self.domain_definition_file = configParser.get('DATABASES', 'domain_definition_file')
        else:
            raise ConfigError('domain_definition_file')            
        # domain_interaction_file
        if configParser.has_option('DATABASES', 'domain_interaction_file'):
            self.domain_interaction_file = configParser.get('DATABASES', 'domain_interaction_file')
        else:
            raise ConfigError('domain_interaction_file')
            
        # check the TMPDIR
        # if a TMPDIR is given as environment variable the tmp directory
        # is created relative to that. This is useful when running on banting
        # (the cluster in the ccbr) and also on Scinet (I might have set the
        # environment variable on Scinet myself though..). Make sure that it
        # points to '/dev/shm/' on Scinet
        childProcess = subprocess.Popen('echo $TMPDIR', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, error = childProcess.communicate()
        TMPDIR_CLUSTER = result.strip()
        try:
            if TMPDIR_CLUSTER[-1] == '/':
                # the tmpPath might be given as an absolute Path
                # thus the last '/' has to be removed
                self.TMPDIR_CLUSTER = TMPDIR_CLUSTER[:-1]
            else:
                self.TMPDIR_CLUSTER = TMPDIR_CLUSTER
        except IndexError:
            self.TMPDIR_CLUSTER = TMPDIR_CLUSTER
        
        if tmpPath[0] == '/':
            # i.e. tmpPath is given as an absolute Path
            self.tmpPath = self.TMPDIR_CLUSTER + tmpPath
        else:
            self.tmpPath = self.TMPDIR_CLUSTER + '/' + tmpPath


        # create the tmp directories and copy the binary files
        self.__prepareTMP()
        self.__prepareOutputPaths()
        
        # load the databases
#        self.threeDID_database = get3DID(self.threeDidFile, self.PfamBoundaryCorrection)
#        self.core_template_database = core_Templates(self.coreTemplatesFile)
        
        self.UniprotSequence = UniprotSequence(self.uniprotSequenceDatabase)
        self.DomainDefinition = DomainDefinition(self.domain_definition_file)
        self.DomainInteraction = DomainInteraction(self.domain_interaction_file)
        self.PDBResolution = PDBResolution(self.pdb_resolution_file)
        
        self.ProteinDefinition = ProteinDefinition()
        self.ProteinInteraction = ProteinInteraction(self.protein_interaction_file,
                                                     self.ProteinDefinition,
                                                     self.DomainInteraction)

        # set the matrix for the substitution score
        self.matrix = self.__selectMatrix(self.matrix_option)
        
        # if running on the cluster copy the database to the tmp DIR for
        # speedup and to avoid killing the network. BLAST is very IO intensive
        # and you don't want that to be run over the network!
        #
        # I can distinguish where the pipeline is running by checking the username
        # you will have to adjust that!
        # my usernames are:
        # local: niklas
        # banting: nberliner
        # Scinet: joan
        childProcess = subprocess.Popen('whoami', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        whoami, error = childProcess.communicate()
        if whoami.strip() == 'strokach':
            # when running the pipeline on beagle or banting
            system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                'cp -r /home/kimlab1/strokach/ncbi-blast-2.2.28+/pdbaa_db ' + \
                                self.tmpPath + 'blast/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, error = childProcess.communicate()
            assert childProcess.returncode == 0
            if childProcess.returncode != 0:
                print 'couldnt copy the blast database!!\n\n\n\n'
        if whoami.strip() == 'alexey':
            # when running the pipeline locally there is no need to copy the database
            # a symlink is enough
            system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                'cd ' + self.tmpPath + 'blast && ' + \
                                'ln -s /home/kimlab1/strokach/ncbi-blast-2.2.28+/pdbaa_db'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, error = childProcess.communicate()
            assert childProcess.returncode == 0
        if whoami.strip() == 'joan':
            # for scinet, blast is already installed, but the database needs to be copied
            system_command = 'mkdir -p ' + self.tmpPath + 'blast && ' + \
                                'cp -r $HOME/niklas-pipeline/blastdb/pdbaa_db ' + \
                                self.tmpPath + 'blast/'
            childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            result, error = childProcess.communicate()
            assert childProcess.returncode == 0
            
    
    
    def __prepareTMP(self):
        # create the basic tmp directory
        # delete its content if it exists
        if not os.path.isdir(self.tmpPath):
            subprocess.check_call('mkdir -p ' + self.tmpPath, shell=True)
        else:
            if not self.tmpPath[-1] == '/':
                self.tmpPath = self.tmpPath +'/'
            subprocess.check_call('rm -r ' + self.tmpPath + '*', shell=True)
        
        
        for i in range(1, self.num_consumers + 1):
            # the consumers
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i)):
                subprocess.check_call('mkdir ' + self.tmpPath + 'Consumer-' + str(i), shell=True)
                
            # tcoffee
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/tcoffee'):
                # create tmp for tcoffee
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee/tmp && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee/lck && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/tcoffee/cache'
                subprocess.check_call(mkdir_command, shell=True)
        
            # FoldX
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/FoldX'):
                # make the directories
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/FoldX'
                # copy the executables
                cp_command = 'cp ' + self.executables + 'FoldX.linux64 ' + self.tmpPath + 'Consumer-' + str(i) + '/FoldX/ && ' + \
                           'cp ' + self.executables + 'rotabase.txt ' + self.tmpPath + 'Consumer-' + str(i) + '/FoldX/'
                # call the command
                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
            
            # modeller
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/modeller'):
                # create workingfolder for modeller
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/modeller'
                subprocess.check_call(mkdir_command, shell=True)
            
            # create tmp for KNOT
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/KNOT'):
                # make the directories
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/KNOT'
                # copy the executables
                cp_command = 'cp ' + self.executables + 'topol ' + self.tmpPath + 'Consumer-' + str(i) + '/KNOT'
                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
            
            # create tmp for pops
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/pops'):
                # make the directories
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/pops'
                # copy the executables
                cp_command = 'cp ' + self.executables + 'pops ' + self.tmpPath + 'Consumer-' + str(i) + '/pops'
                subprocess.check_call(mkdir_command + ' && ' + cp_command, shell=True)
            
            # create tmp for output
            if not os.path.isdir(self.tmpPath + 'Consumer-' + str(i) + '/output'):
                # make the directories
                mkdir_command = 'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/output && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/output/alignments && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/output/bestModels && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/output/pdbFiles && ' + \
                                'mkdir ' + self.tmpPath + 'Consumer-' + str(i) + '/output/pickled'
                subprocess.check_call(mkdir_command, shell=True)


    def __prepareOutputPaths(self):
        if not os.path.isdir(self.outputPath):
            subprocess.check_call('mkdir ' + self.outputPath, shell=True)
        # Files are compressed into individual tar archives, do not need these paths anymore
#        if not os.path.isdir(self.outputPath + 'alignments/'):
#            subprocess.check_call('mkdir ' + self.outputPath + 'alignments/', shell=True)
#        if not os.path.isdir(self.outputPath + 'bestModels/'):
#            subprocess.check_call('mkdir ' + self.outputPath + 'bestModels/', shell=True)
#        if not os.path.isdir(self.outputPath + 'pdbFiles/'):
#            subprocess.check_call('mkdir ' +  self.outputPath + 'pdbFiles/', shell=True)

        
    def __selectMatrix(self, matrix):
        if matrix == 'blosum80':
            return MatrixInfo.blosum80
        if matrix == 'blosum60':
            return MatrixInfo.blosum60
        else:
            print 'specified matrix not found!'
        
        

    def run(self):
                    
        ## Part for multiprocessing ##
        # see: http://doughellmann.com/2009/04/pymotw-multiprocessing-part-1.html
        # Establish communication queues
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        
        # create a logger instance
        # I started to play with the logger a little bit but didn't have the time
        # to fully implement it to be really usefull. It works, one just has
        # to set the desired logging with the information where needed
        logger = MultiProcessingLog(self.outputPath + self.name + '.log', mode='w', maxsize=0, rotate=1)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        logger.setFormatter(formatter)
        
        log = logging.getLogger(__name__)
        log.addHandler(logger)
        if self.DEBUG:
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.INFO)
        
        # create the pool to control the number of t_coffee instances
        # that are allowed to run in parallel
        pool = ActivePool()
        s = multiprocessing.Semaphore(self.tcoffee_parallel_runs)
        
        # Start consumers
        print 'Creating %d consumers' % self.num_consumers
        proc_name = [ 'Consumer-' + str(i) for i in range(1, self.num_consumers + 1) ]
        consumers = [ Consumer(proc_name[i-1], tasks, results, self.runTime, pool, s, self.DEBUG, self.outputPath, logger, webServer=self.webServer)
                      for i in range(1, self.num_consumers + 1) ]

        for w in consumers:
            w.start()
    
        # check if a skip file is present
        # I normally didn't use it but the idea is that if you run the pipeline
        # and it reached the maximum runtime, you could easily re-run and skip
        # the already calculated entries.
        if os.path.isfile('processed.log'):
            skipFile = open('processed.log', 'r')
            skip = list()
            for line in skipFile:
                skip.append(line.strip())
            skipFile.close()
            SKIP = True
        else:
            SKIP = False
        
        num_jobs = 0
        with open(self.inputFile, 'r') as f:
            for l in f:
                if l[0][0] == ' ' or l[0][0] == '\t':
                    continue
                
                line = [ ll.strip() for ll in l.split('\t') ]
                
                # AS: Mutation does not necessarily have to be specified
                if len(line) > 1:
                    uniprotKB, mutation = line[0], line[1]
                elif len(line) == 1:
                    uniprotKB = line[0]
                    mutation = ''
                    
                # check if some entries should be skipped
                if SKIP:
                    if uniprotKB + mutation in skip:
                        print 'skipping ', uniprotKB, mutation
                        continue
                
                log.debug('Added to queue uniprot %s with mutation %s' % (uniprotKB, mutation) )
                
                # Enqueue jobs                
                num_jobs += 1
                tasks.put( Task(self.mutation_uniprot,
                                uniprotKB, 
                                mutation,
                                self.savePDB,
                                self.tmpPath,
                                self.outputPath,
                                self.pdbPath,
                                self.matrix,
                                self.gap_start,
                                self.gap_extend,
                                self.modeller_runs,
                                self.buildModel_runs,
                                self.foldX_WATER,
                                self.DomainDefinition,
                                self.DomainInteraction,
                                self.UniprotSequence,
                                self.PDBResolution,
                                self.ProteinDefinition,
                                self.ProteinInteraction,
                                )
                          )
                          

        
        # Add a poison pill for each consumer
        for i in range(1, self.num_consumers+1):
            tasks.put( None )
       
        # Wait for all of the tasks to finish
        tasks.join()

        # process the result
        res_wt = open(self.outputPath + 'result_wt.log', 'w')
        res_mut = open(self.outputPath + 'result_mut.log', 'w')
        res_if = open(self.outputPath + 'result_additional_information.log', 'w')
        skiplog = open(self.outputPath + 'not_processed.log', 'w')
        
        # Write header
        id_labels = 'uniprotIDs\t' + 'pfamIDs\t' + 'domain_defs\t' + 'mutation\t' + 'wt_or_mut\t'
                    
        value_labels = 'normDOPE\t' + \
                    'intraclashes_energy1\t' + 'intraclashes_energy2\t' + \
                    'interaction_energy\t' + 'backbone_hbond\t' + \
                    'sidechain_hbond\t' + 'van_der_waals\t' + \
                    'electrostatics\t' + 'solvation_polar\t' + \
                    'solvation_hydrophobic\t' + 'Van_der_Waals_clashes\t' + \
                    'entropy_sidechain\t' + 'entropy_mainchain\t' + \
                    'sloop_entropy\t' + 'mloop_entropy\t' + 'cis_bond\t' + \
                    'torsional_clash\t' + 'backbone_clash\t' + 'helix_dipole\t' + \
                    'water_bridge\t' + 'disulfide\t' + 'electrostatic_kon\t' + \
                    'partial_covalent_bonds\t' + 'energy_ionisation\t' + \
                    'entropy_complex\t' + 'number_of_residues\t' + \
                    'stability_energy\t' + 'stability_backbone_hbond\t' + \
                    'stability_sidechain_hbond\t' + 'stability_Van_der_Waals\t' + \
                    'stability_electrostatics\t' + 'stability_solvation_polar\t' + \
                    'stability_solvation_hydrophobic\t' + 'stability_Van_der_Waals_clashes\t' + \
                    'stability_entropy_sidechain\t' + 'stability_entropy_mainchain\t' + \
                    'stability_sloop_entropy\t' + 'stability_mloop_entropy\t' + \
                    'stability_cis_bond\t' + 'stability_torsional_clash\t' + \
                    'stability_backbone_clash\t' + 'stability_helix_dipole\t' + \
                    'stability_water_bridge\t' + 'stability_disulfide\t' + \
                    'stability_electrostatic_kon\t' + 'stability_partial_covalent_bonds\t' + \
                    'stability_energy_ionisation\t' + 'stability_entropy_complex\t' + \
                    'stability_number_of_residues\n'
                    
        value_labels_extra = 'core_or_interface\t' + 'seq_id_avg\t' + \
                    'seq_id_chain1\t' + 'seq_id_chain2\t' + \
                    'matrix_score\t' + 'if_hydrophobic\t' + \
                    'if_hydrophilic\t' + 'if_total\t' + \
                    'contactVector_wt_ownChain\t' + 'contactVector_wt\t' + \
                    'contactVector_mut_ownChain\t' + 'contactVector_mut\t' + \
                    'secondary_structure_wt\t' + 'solvent_accessibility_wt\t' + \
                    'secondary_structure_mut\t' + 'solvent_accessibility_mut\n'
        
        res_wt.write(id_labels + value_labels)
        res_mut.write(id_labels + value_labels)
        res_if.write(id_labels + value_labels_extra)

        # Start printing results
        while num_jobs:
            num_jobs -= 1
            # if the timeout was reached and the calculation stopped before
            # every task was calculated the queue is empty before num_jobs is 0
            
            # check if the que is empty
            if results.empty():
                break
            
            output_data = results.get()
            log.debug('output_data: \n%s' % output_data)            
            for output_dict in output_data:
                if isinstance(output_dict, list) and (output_dict[0] == None):
                    skiplog.write(output_dict[1] + '\t' + output_dict[2] + '\n')
                    continue
                elif isinstance(output_dict, list) and (output_dict[0] == 'no template found'):
                    skiplog.write(output_dict[1] + '\t' + 'no template found' + '\n')
                    continue
                else:
                    # Add sequences that were not present already in the database
                    if output_dict.has_key('new_sequences'):
                        self.get_uniprot_sequence.add(output_dict['new_sequences'])
                        
                    # Unique identifier for each line
                    id_data = ('-'.join(output_dict['uniprotIDs']) + '\t' + 
                             '-'.join(output_dict['pfamIDs']) + '\t' + 
                             '_'.join(['-'.join([str(i) for i in x]) for x in output_dict['domain_defs']]) + '\t' + 
                             output_dict['mutation'] + '\t')
#
#                    # Unique identifier for each line
#                    id_data =   output_dict['uniprotIDs'][0] + '\t' + \
#                                output_dict['uniprotIDs'][1] + '\t' + \
#                                output_dict['pfamIDs'][0] + '\t' + \
#                                output_dict['pfamIDs'][1] + '\t' + \
#                                str(output_dict['domain_defs'][0][0]) + '-' + str(output_dict['domain_defs'][0][1]) + '\t' + \
#                                str(output_dict['domain_def2'][1][0]) + '-' + str(output_dict['domain_def2'][1][1]) + '\t' + \
#                                output_dict['mutation'] + '\t'
                                          
                    # Make line for wildtype file
                    resForWrite_wt =  id_data + 'wt\t' + \
                                        str(output_dict['normDOPE_wt'])[:6] + '\t' +\
                                        '\t'.join(output_dict['AnalyseComplex_energy_wt']) + '\t' +\
                                        '\t'.join(output_dict['Stability_energy_wt']) + \
                                        '\n'
                    
                    # Make line for mutant file
                    resForWrite_mut = id_data + 'mut\t' + \
                                        '-' + '\t' + \ # mutant structure has no normDOPE score (same as wildtype)
                                        '\t'.join(output_dict['AnalyseComplex_energy_mut']) + '\t' + \
                                        '\t'.join(output_dict['Stability_energy_mut']) + \
                                        '\n'
                    
                    # Make line for additional information file
                    resForWrite_if = id_data + '-\t' + \
                                str(output_dict['is_in_core']) + '\t' + \
                                '\t'.join([str(i) for i in output_dict['alignment_scores']]) + '\t' + \
                                 str(output_dict['matrix_score']) + '\t' + \
                                 '\t'.join(output_dict['interface_size']) + '\t' + \
                                 ','.join(output_dict['physChem_wt_ownChain']) + '\t' + \
                                 ','.join(output_dict['physChem_wt']) + '\t' + \
                                 ','.join(output_dict['physChem_mut_ownChain']) + '\t' + \
                                 ','.join(output_dict['physChem_mut']) + '\t' + \
                                 str(output_dict['secondary_structure_wt'])  + '\t' + \
                                 str(output_dict['solvent_accessibility_wt']) + '\t' + \
                                 str(output_dict['secondary_structure_mut'])  + '\t' + \
                                 str(output_dict['solvent_accessibility_mut']) + \
                                 '\n'

                    # Make another file to keep precalculated values?
                    # resForWrite_precalc = []
                    # output_dict['interactingAA']
                    # output_dict['surfaceAA']
                                 
                    # Write output lines             
                    res_wt.write(str(resForWrite_wt))
                    res_mut.write(str(resForWrite_mut))
                    res_if.write(str(resForWrite_if))                 
        
        # save the database (needed if new sequences were added)
        self.get_uniprot_sequence.close()
        
        # close and flush the output files        
        res_wt.close()
        res_mut.close()
        res_if.close()
        skiplog.close()
        logger.close()
        
        # save the results from ramdisk
        if self.saveScinet:
            scinetCleanup(self.outputPath, self.saveTo, self.name)



if __name__ == '__main__':
    # read which configFile to use    
    optParser = optparse.OptionParser()
    optParser.add_option('-c', '--config', action="store", dest='configFile')
    options, args = optParser.parse_args()
       
    configFile = options.configFile
    
    if not os.path.isfile(configFile):
        print 'Error: configFile not found!'
        print 'exiting'
    else:
        try:
            p = pipeline(configFile)
            p.run()
        except DataError, e:
            print 'Error: input file',  e.inputFile, ' not found!'
            print 'exiting...'
        except ConfigError, e:
            print 'Error: option', e.option, ' not found!'
            print 'exiting...'