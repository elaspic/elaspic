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

        
    def __call__(self, uniprot_id):
        """
        returns the uniprot sequence. If the sequence is not in the database it
        tries to retrieve it from the uniprot website.
        
        Note: retrieval from the website does not work on Scinet!
        
        """
        if uniprot_id in ['A6NF79', 'C9JUS1', 'Q6N045', 'A6NMD1']:
            # these uniprotKBs made problems
            return 'no sequences'
        elif uniprot_id == 'P02735':
            # this sequence got replaced. I don't know right now if I can take
            # replaced sequence so I rather dismiss it.
            return 'no sequences'
        
        # the True/False value is used to add new sequences to the database in
        # end. Only usefull if you run one instance at a time otherwise you will
        # get diverging sequence databases.
        try:
            return self.uniprot_data[uniprot_id]
        except KeyError:
            childProcess = subprocess.Popen('whoami', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            whoami, error = childProcess.communicate()
            if whoami.strip() == 'joan':
                print 'uniprot sequence not found'
                return None
            else:
                print 'Fetching uniprot sequence', uniprot_id, 'from server'
                address = 'http://www.uniprot.org/uniprot/' + uniprot_id + '.fasta'
                handle = urllib2.urlopen(address)
                sequence = next(SeqIO.parse(handle, "fasta"))
                sequence.id = uniprot_id
                self.uniprot_data[uniprot_id] = sequence
                return sequence

    
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


class MyDatabase(object):
    
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


    ###########################################################################

    def get_dicts(self, key):
        """
        """
        query_df = self.__call__(key)
        
        if not query_df or query_df == [None, None]:
            return None
        
        if type(query_df) is pd.DataFrame:
            list_of_dicts = []
            for idx, df in query_df.iterrows():
                row_dict = df.to_dict()
                row_dict.update({'idx': idx})
                list_of_dicts.append(row_dict)
            return list_of_dicts
        
        
        elif type(query_df) is list and len(query_df) == 2:
            list_of_dicts = []
            query_df_forward, query_df_reversed = query_df
            if query_df_forward:
                for idx, df in query_df_forward.iterrows():
                    row_dict = df.to_dict()
                    row_dict.update({'idx': idx, 'reversed': False})
                    list_of_dicts.append(row_dict)
            if query_df_reversed:
                for idx, df in query_df_reversed.iterrows():
                    row_dict = df.to_dict()
                    row_dict.update({'idx': idx, 'reversed': True})
                    list_of_dicts.append(row_dict)
            return list_of_dicts
            
        else:
            raise Exception()
            
    
    
    def update_from_df(self, modified_df):
        """ Updates the internal data frame with the modifications present in the
        copied dataframe
        """
        self.db.update(modified_df)
        
        

    def update_from_dicts(self, list_of_dicts):
        """
        """
        index = []
        for row_dict in list_of_dicts:
            index = row_dict.pop('idx')
            if row_dict.has_key('reversed'):
                del row_dict['reversed']
        
        modified_df = pd.DataFrame(list_of_dicts, index=index)
        
        self.update_from_df(modified_df)      


    ###########################################################################

    def load_db(self):
        if isfile(self.db_name + '.pickle') and isfile(self.db_name + '_updates.pickle'):
            print 'Loading %s...' % ' '.join(self.db_name.split('_'))
            self.db = pd.read_pickle(self.db_name + '.pickle')
            self.db_updates = pd.read_pickle(self.db_name + '_updates.pickle')
            return True
        else:
            return False
    
    
    def save_db(self):
        self.db.to_pickle(self.db_name + '.pickle')
        self.db_updates.to_pickle(self.db_name + '_updates.pickle')



class DomainDefinition(MyDatabase):
    """ Contains pdbfam-based definitions of all pfam domains in the pdb
    """
    
    def __init__(self, textfile_name):
        """
        """
        self.db_name = 'domain_definition'
        
        if not self.load_db():
            # No precompiled database exists
            print 'Generating %s database...' % ' '.join(self.db_name.split('_'))
            
            names = ['pdb_id', 'autoPFAMA', 'pfam_name', 'pdb_chain', 'pdb_definition', 'cath_id']
            
            db = pd.read_csv(textfile_name, sep='\t', header=False, names=names)
            
            db['pdb_definition'] = db['pdb_definition'].apply(self.split_domain_semicolon)
            
            self.db = db
            self.db_updates = pd.DataFrame()
            self.save_db()
            
#        self.grouped_db = self.db.groupby('pfam_name')
                    
                    
    def __call__(self, pfam_name):
        return self.db[self.db['pfam_name'] == pfam_name]
#        if self.grouped_db.groups.has_key(pfam_name):
#            return self.grouped_db.get_group(pfam_name)
#        else:
#            print "No entries for the given pfam_name:", pfam_name
#            return None



class DomainInteraction(MyDatabase):
    """ Keeps the domain-domain interaction information from pdbfam
    """
    
    def __init__(self, textfile_name):
        """ Create the internal dictionary using a tsv file exported from mysql
        """
        self.db_name = 'domain_interaction'
        
        # Check if a pickled database exists
        if not self.load_db():    
            # No precompiled database exists
            print 'Generating %s database...' % ' '.join(self.db_name.split('_'))

            names = [
            'pdb_id', 'pfam_name_1', 'pdb_chain_1', 'pdb_definition_1', 'pdb_interface_aa_1', 'cath_id_1', 
            'pfam_name_2', 'pdb_chain_2', 'pdb_definition_2', 'pdb_interface_aa_2', 'cath_id_2', 'domain_contact_id',]
            
            db = pd.read_csv(textfile_name, sep='\t', header=False, names=names)
            
            db['pdb_definition_1'] = db['pdb_definition_1'].apply(self.split_domain_semicolon)
            db['pdb_definition_2'] = db['pdb_definition_2'].apply(self.split_domain_semicolon)
            
            db['pdb_interface_aa_1'] = db['pdb_interface_aa_1'].apply(self.split_interface_aa)
            db['pdb_interface_aa_2'] = db['pdb_interface_aa_2'].apply(self.split_interface_aa)
            
            self.db = db
            self.db_updates = pd.DataFrame()
            self.save_db()
            
#        self.grouped_db = self.db.groupby(['pfam_name_1', 'pfam_name_2'])
        

    def __call__(self, pfam_names):
        """
        Note that the produced dataframe may not have the same order as the keys
        """
        key1 = tuple([pfam_names[0], pfam_names[1]])
        key2 = tuple([pfam_names[1], pfam_names[0]])
            
        if self.grouped_db.groups.has_key(key1):
            return_df_1 = self.grouped_db.get_group(key1)
        else:
            return_df_1 = None
            
        if self.grouped_db.groups.has_key(key2):
            return_df_2 = self.grouped_db.get_group(key2)
        else:
            return_df_2 = None
            
        if not return_df_1 and not return_df_2:
            print "No entries for the given domain combination:", pfam_names[0], pfam_names[1]
            
        return [return_df_1, return_df_2]



class ProteinDefinition(MyDatabase):
    """ Initiated using parsed pfam_scan.pl output for the human uniprot (or the entire uniprot)
    The database is updated with information about calculated models
    """
    
    def __init__(self):
        
        # Check if a pickled database exists
        self.db_name = 'protein_definition'
        
        if not self.load_db():
            # Create a database from a file containing fasta sequences, the pfam_scan.pl output from that file
            # Also needs a list of clans, a list os repeating domains and a list of superdomains
            # Filenames are hardcoded at this points
            print 'Generating %s database...' % ' '.join(self.db_name.split('_'))

            pfam_parser = parse_pfamscan.make_uniprot_pfam_database()
            pfam_parser.run()
            db = pfam_parser.get_dataframe()
            
            db['cath_id'] = None
            db['template_errors'] = None
            db['model_id'] = None
            db['modelling_errors'] = None
            
            self.db = db
            print self.db
            self.db_updates = pd.DataFrame()
            self.save_db()
            
#        self.grouped_db = self.db.groupby('uniprot_id')


    def __call__(self, uniprot_id):
        """ 
        """
        # using (uniprot_id, uniprot_id,) because in the future we might want to incorporate
        # splicing and other variants, which will be identified by the second uniprot_id
        
        return self.db[self.db['uniprot_id'] == uniprot_id]
        
#        if self.grouped_db.groups.has_key(uniprot_id):
#            return self.grouped_db.get_group(uniprot_id)
#        else:
#            # In the future, make it fetch the uniprot sequence, analyse it with pfam_scan.pl,
#            # and save the output to the database
#            print "No protein definition entries for the given uniprot_id:", uniprot_id
#            return None



class ProteinInteraction(MyDatabase):
    """ Contains known interactions between uniprot proteins
    """
    
    def __init__(self, textfile_name, ProteinDefinition):
        """ 
        Checks if the interaction database is already available as a pickled object.
        Generates it from file 'textfile_name' otherwise.
        """
        self.db_name = 'protein_interaction'
        
        if not self.load_db():
            # No precompiled database exists
            print 'Generating %s database...' % ' '.join(self.db_name.split('_'))
            
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
            
            
            column_uniprot_id_1 = []
            column_uniprot_id_2 = []
            column_pfam_name_1 = []
            column_pfam_name_2 = []
            column_alignment_def_1 = []
            column_alignment_def_2 = []
            
            # Go over each unique protein-protein interaction
            for interaction in set_of_interactions:
                
                interaction = list(interaction)
                # Get a list of pfam domains for each protein in a pair
                
                if len(interaction) == 1:
                    # Homodimer
                    uniprot_id_1 = interaction[0]
                    uniprot_id_2 = uniprot_id_1
                    
                    protein_domains = ProteinDefinition(uniprot_id_1)
                    if protein_domains:
                        pfam_names_1 = list(protein_domains['seq_id'])
                        alignment_defs_1 = list(protein_domains['alignment_def'])
                        pfam_names_2 = pfam_names_1
                        alignment_defs_2 = alignment_defs_1
                    else:
                        number_of_missing_pfams += 1
                        continue
                    
                elif len(interaction) == 2:
                    # Heterodimer
                    uniprot_id_1 = interaction[0]
                    uniprot_id_2 = interaction[1]
                    
                    protein_domains_1 = ProteinDefinition(uniprot_id_1)
                    if protein_domains_1:
                        pfam_names_1 = list(protein_domains_1['seq_id'])
                        alignment_defs_1 = list(protein_domains_1['alignment_def'])
                    else:
                        number_of_missing_pfams += 1
                        continue
                    
                    protein_domains_2 = ProteinDefinition(uniprot_id_2)
                    if protein_domains_2:
                        pfam_names_2 = list(protein_domains_2['seq_id'])
                        alignment_defs_2 = list(protein_domains_2['alignment_def'])
                    else:
                        number_of_missing_pfams += 1
                        continue
                else:
                    raise Exception
                        
                for pfam_name_1, alignment_def_1 in zip(pfam_names_1, alignment_defs_1):
                    for pfam_name_2, alignment_def_2 in zip(pfam_names_2, alignment_defs_2):
                        column_uniprot_id_1.append(uniprot_id_1)
                        column_uniprot_id_2.append(uniprot_id_2)
                        column_pfam_name_1.append(pfam_name_1)
                        column_pfam_name_2.append(pfam_name_2)
                        column_alignment_def_1.append(alignment_def_1)
                        column_alignment_def_2.append(alignment_def_2)
            
            # Turn print statement output back on
            sys.stdout = original_stdout
            print number_of_missing_pfams
            
            # Create a dataframe with the compiled results
            db = pd.DataFrame(index=range(len(column_uniprot_id_1)))
            db['uniprot_id_1'] = column_uniprot_id_1
            db['uniprot_id_2'] = column_uniprot_id_2
            db['domain_name_1'] = column_pfam_name_1
            db['domain_name_2'] = column_pfam_name_2
            db['alignment_def_1'] = column_alignment_def_1
            db['alignment_def_2'] = column_alignment_def_2
            db['cath_id_1'] = None
            db['cath_id_2'] = None
            db['domain_contact_id'] = None
            db['template_errors'] = None
            db['model_path'] = None
            db['modelling_errors'] = None
            
            print db
            db.drop_duplicates(inplace=True)
            print db

            self.db = db
            self.db_updates = pd.DataFrame()
            self.save_db()
        
#        self.grouped_db_1 = self.db.groupby(['uniprot_id_1'])
#        self.grouped_db_2 = self.db.groupby(['uniprot_id_2'])
    
    
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



class ProteinCoreMutation(MyDatabase):
    """ Stores information about mutations in the protein core
    Unimplemented functionality
    """
    def __init__(self):
        self.db_name = 'protein_core_mutation'



class ProteinInterfaceMutation(MyDatabase):
    """ Stores information about mutations in the protein interface
    Unimplemented functionality
    """
    def __init__(self):
        self.db_name = 'protein_interface_mutation'




class pipeline():
    
    def __init__(self, configFile):
        
        #######################################################################  
        # read the configuration file and set the variables
        
        configParser = SafeConfigParser(
                                defaults={'tmpPath':'/tmp/pipeline/',
                                          'HETATM':True,
                                          'outputPath':'results/',
                                          'DEBUG':False,
                                          'saveTo':'$SCRATCH/niklas-pipeline/',
                                          'saveScinet':False,
                                          'path_to_archive': '/home/kimlab1/database_data/elaspic/',
                                          'webServer': False,
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
        # path_to_archive
        if configParser.has_option('DEFAULT', 'path_to_archive'):
            self.path_to_archive = configParser.get('DEFAULT', 'path_to_archive')
        else:
            self.path_to_archive = defaults['path_to_archive']        
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
        # mutation_uniprot ### should be renamed template_finding
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
        
        self.uniprot_sequence_database = UniprotSequence(self.uniprotSequenceDatabase)
        self.domain_definition_database = DomainDefinition(self.domain_definition_file)
        self.domain_interaction_database = DomainInteraction(self.domain_interaction_file)
        self.pdb_resolution_database = PDBResolution(self.pdb_resolution_file)
        
        self.protein_definition_database = ProteinDefinition()
        self.protein_interaction_database = ProteinInteraction(self.protein_interaction_file,
                                                     self.protein_definition_database)

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
                                self.uniprot_sequence_database,
                                self.domain_definition_database,
                                self.domain_interaction_database,
                                self.pdb_resolution_database,
                                self.protein_definition_database,
                                self.protein_interaction_database,
                                self.path_to_archive
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
                    resForWrite_wt = (id_data + 'wt\t' +
                                        str(output_dict['normDOPE_wt'])[:6] + '\t' +
                                        '\t'.join(output_dict['AnalyseComplex_energy_wt']) + '\t' +
                                        '\t'.join(output_dict['Stability_energy_wt']) + 
                                        '\n')
                    
                    # Make line for mutant file
                    resForWrite_mut = (id_data + 'mut\t' + 
                                        '-' + '\t' + # mutant structure has no normDOPE score (same as wildtype)
                                        '\t'.join(output_dict['AnalyseComplex_energy_mut']) + '\t' +
                                        '\t'.join(output_dict['Stability_energy_mut']) +
                                        '\n')
                    
                    # Make line for additional information file
                    resForWrite_if = (id_data + '-\t' +
                                str(output_dict['is_in_core']) + '\t' +
                                '\t'.join([str(i) for i in output_dict['alignment_scores']]) + '\t' +
                                 str(output_dict['matrix_score']) + '\t' +
                                 '\t'.join(output_dict['interface_size']) + '\t' +
                                 ','.join(output_dict['physChem_wt_ownChain']) + '\t' +
                                 ','.join(output_dict['physChem_wt']) + '\t' +
                                 ','.join(output_dict['physChem_mut_ownChain']) + '\t' +
                                 ','.join(output_dict['physChem_mut']) + '\t' +
                                 str(output_dict['secondary_structure_wt'])  + '\t' +
                                 str(output_dict['solvent_accessibility_wt']) + '\t' +
                                 str(output_dict['secondary_structure_mut'])  + '\t' +
                                 str(output_dict['solvent_accessibility_mut']) +
                                 '\n')

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