# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 22:16:38 2013

@author: alexey
"""

from Bio import SeqIO, AlignIO

import class_callTcoffee as callTcoffee
import class_multi as multi
from Bio.SubsMat import MatrixInfo
from pipeline_splicing import getInteractions, get3DID, getUniprotSequence, pdb_resolution, core_Templates

import pandas as pd
import numpy as np
import pylab as pl
import cPickle as pickle
from os.path import isfile
import urllib2

import csv


class splicingSequenceDatabaseByUniprotKB:
    
    def __init__(self):     
        
        database_name = 'splicing_sequences.pickle'
        
        if isfile(database_name):
            print 'Loading splicing sequence database'
            handle = open(database_name, 'rb')
            splicing_db = pickle.load(handle)
            handle.close()
            print 'it contains', len(splicing_db), 'uniprot IDs'
        else:
            splicing_db = dict()

        
    def __call__(self, uniprot_accession):
        
        if splicing_db.has_key(uniprot_accession):
            return splicing_db[uniprot_accession]
        else:
            print 'No such key found!'

    
    def add(self, uniprot_accession, seqio_objects):
        """
        add the new items (which is a list of tuples) to the database
        """
        if not splicing_db.has_key(uniprot_accession):
            try:
                address = 'http://www.uniprot.org/uniprot/' + uniprot_accession + '.fasta'
                handle = urllib2.urlopen(address)
                seqio_object = next(SeqIO.parse(handle, "fasta"))
                splicing_db[uniprot_accession] = [seqio_object,]
                print 'Created canonical entry for ', uniprot_accession
            except:
                print 'Could not create canonical entry for ', uniprot_accession
                return uniprot_accession
        for seqio_object in seqio_objects:
            if seqio_object not in splicing_db[uniprot_accession]:
                splicing_db[uniprot_accession].append(seqio_object)
                print 'Added splice sequence to ', uniprot_accession, ' database'
            else:
                print 'Splice sequence already present'
                return seqio_object
        
    def export_sequences(self):
        for dict_key, dict_values in splicing_db.iteritems():
            output_handle = open('/home/alexey/working/pipeline-splicing/lib/splicing-sequences-by-uniprotKB/' +  dict_key + '.fasta', 'w')
            SeqIO.write(dict_values, output_handle, 'fasta')
            output_handle.close()
    
    def allign_in_tcoffee(self, seqIDs):
        """
        Copied from 'class_get_uniprot_template_core_and_interface.py'
        
        """
        if seqIDs[0] == 'all':
            seqIDs = [dict_key for dict_key in splicing_db]
        
        # use the try statement to ensure that the sempaphore is released
        # in the end
        tcoffee = callTcoffee.tcoffee_alignment('/tmp/pipeline-splicing/', # tmpPath, 
                                       'Consumer-1', # unique
                                       '/home/alexey/working/pipeline-splicing/lib/splicing-alignments-by-uniprotKB/', # saveAlignments = HOME + outputPath + 'alignments/'
                                       ['/home/alexey/working/pipeline-splicing/lib/splicing-sequences-by-uniprotKB/' + seqID + '.fasta' for seqID in seqIDs], # seqfiles
                                       seqIDs,  # seqids
                                       )
        alignments = tcoffee.align()
    
    def close(self):
        """
        save the database back to disk. Do it at the end if
        new sequences where added.
        """
        print 'saving the splicing database'
        f = open(database_name, 'wb')    
        pickle.dump(splicing_db, f)
        f.close()




class splicingSequenceDatabaseBySplicingID:
    
    def __init__(self):     
        
        self.database_name = 'splicing_sequence_alignmnets_by_splicingID.pickle'
        self.tmpPath = '/tmp/pipeline-splicing/'
        self.unique = 'Consumer-1'
        self.saveSequencesPath = '/home/alexey/working/pipeline-splicing/lib/splicing-sequences-by-splicingID/'
        self.saveAlignmentsPath = '/home/alexey/working/pipeline-splicing/lib/splicing-alignments-by-splicingID/'
        
        if isfile(self.database_name):
            print 'Loading splicing sequence alignments database'
            handle = open(self.database_name, 'rb')
            self.splicing_db = pickle.load(handle)
            handle.close()
            print 'it contains', len(self.splicing_db), 'uniprot IDs'
        else:
            self.splicing_db = dict()

        
    def __call__(self, splicing_accession):
        
        if self.splicing_db.has_key(splicing_accession):
            return self.splicing_db[splicing_accession]
        else:
            print 'No such key found!'

    
    def add(self, (uniprot_accession, splicing_accession,), seqio_object2):
        """
        add the new items (which is a list of tuples) to the database
        """
        if self.splicing_db.has_key((uniprot_accession, splicing_accession,)):
            return splicing_accession
        else:
            try:
                self.splicing_db[(uniprot_accession, splicing_accession,)] = []
                
                address = 'http://www.uniprot.org/uniprot/' + uniprot_accession + '.fasta'
                handle = urllib2.urlopen(address)
                seqio_object1 = next(SeqIO.parse(handle, "fasta"))
                self.splicing_db[(uniprot_accession, splicing_accession,)].append(seqio_object1)
                print 'Created canonical entry for ', uniprot_accession,  ' ', splicing_accession    

                self.splicing_db[(uniprot_accession, splicing_accession,)].append(seqio_object2)
                print 'Added splice sequence to ', uniprot_accession, ' ', splicing_accession, ' database'
                
            except:
                print 'Could not create canonical entry for ', uniprot_accession,  ' ', splicing_accession
                return splicing_accession

        
    def export_sequences(self):
        for dict_key, dict_values in self.splicing_db.iteritems():
            output_handle = open(self.saveSequencesPath +  dict_key[0] + '_' + dict_key[1] + '.fasta', 'w')
            SeqIO.write(dict_values, output_handle, 'fasta')
            output_handle.close()

    
    def allign_in_tcoffee(self):
        """
        Copied from 'class_get_uniprot_template_core_and_interface.py'
        """
        # use the try statement to ensure that the sempaphore is released
        # in the end
        for dict_key in self.splicing_db:
            tcoffee = callTcoffee.tcoffee_alignment(self.tmpPath, # tmpPath, 
                                           self.unique, # unique
                                           self.saveAlignmentsPath, # saveAlignments = HOME + outputPath + 'alignments/'
                                           [self.saveSequencesPath + dict_key[0] + '_' + dict_key[1] + '.fasta', ], # seqfiles
                                           [dict_key, ],  # seqids
                                           )
            alignments = tcoffee.align()                                      
            self.splicing_db[dict_key].append(alignments[0]) # made after the initial run


    def analyse_alignments(self):
        for dict_key, dict_value in self.splicing_db.iteritems():
            seq1, seq2 = list(self.splicing_db[dict_key][2][0])
            assert len(seq1) == len(seq2)
            resID = 0
            splicing_mutation_positions = set()
            for i in range(len(seq1)):
                if seq1[i] != '-':
                    resID += 1
                if seq1[i] != seq2[i]:
                    if seq1[i] != '-':
                        splicing_mutation_positions.add(resID)
                    else:
                        splicing_mutation_positions.add(resID)
                        splicing_mutation_positions.add(resID+1)
            self.splicing_db[dict_key].append(splicing_mutation_positions)

    
    def export_pickle(self):
        """
        save the database back to disk. Do it at the end if
        new sequences where added.
        """
        print 'saving the splicing database'
        f = open(self.database_name, 'wb')    
        pickle.dump(self.splicing_db, f)
        f.close()


def appendSplicingData(accessions):    
    splicing_db = splicingSequenceDatabaseBySplicingID()
    error_list = []
    for accession in set(accessions):
        return_value = splicing_db.add(accession, list(splicing_df[splicing_df['uniprot_accession'] == accession]['seqio_object']))
        if return_value:
            error_list.append(return_value)




# load splicing dataset            
splicing_db = splicingSequenceDatabaseBySplicingID()


# Load Niklas' databases
handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_3DID_database.pickle', 'rb')
pipeline_3DID_database = pickle.load(handle)
handle.close()

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_core_templates.pickle', 'rb')
pipeline_core_templates = pickle.load(handle)
handle.close()

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_interaction_database.pickle', 'rb')
pipeline_interaction_database = pickle.load(handle)
handle.close()

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_uniprot_sequences.pickle', 'rb')
pipeline_uniprot_sequences = pickle.load(handle)
handle.close()


# Parse splice sequence file
handle = open('/home/alexey/working/pipeline-splicing/databases/4alexey.fasta')

uniprot_accession = []
splicing_accession = []
splicing_accession_prefix = []
splice_confidence = []
seqio_object = []
seqio_length = []


for sequence in SeqIO.parse(handle, 'fasta'):
    uniprot_accession.append(sequence.id.split('//')[1][:-3])
    splicing_accession.append(sequence.id.split('//')[0])
    splicing_accession_prefix.append(sequence.id.split('//')[0][0])
    splice_confidence.append(float(sequence.description.split('\t')[-1]))
    seqio_object.append(sequence)
    seqio_length.append(len(sequence.seq))  

pipeline_core = []
pipeline_interaction = []

for key in uniprot_accession:
    if key in pipeline_core_templates:
        pipeline_core.append(True)
    else:
        pipeline_core.append(False)
    if key in pipeline_interaction_database:
        pipeline_interaction.append(True)
    else:
        pipeline_interaction.append(False)

splicing_df_all = pd.DataFrame({'uniprot_accession': uniprot_accession,
                              'splicing_accession': splicing_accession,
                              'inc_exc': splicing_accession_prefix,
                              'splice_confidence': splice_confidence,
                              'seqio_object': seqio_object,
                              'seqio_length': seqio_length,
                              'pipeline_core': pipeline_core,
                              'pipeline_interaction': pipeline_interaction})

# Get the number of splice variats with given properties
print len(splicing_df_all[(splicing_df_all['pipeline_interaction'] == True)])

print len(splicing_df_all[(splicing_df_all['splice_confidence'] > 0.60)])

print len(splicing_df_all[(splicing_df_all['seqio_length'] < 3000)])
    
print len(splicing_df_all[(splicing_df_all['pipeline_interaction'] == True) &
    (splicing_df_all['splice_confidence'] > 0.60)])

print len(splicing_df_all[(splicing_df_all['pipeline_interaction'] == True) &
    (splicing_df_all['splice_confidence'] > 0.60) & 
    (splicing_df_all['seqio_length'] < 3000)])
    
print len(splicing_df_all)

## Subset all splicing variants to thoes with desired properties   
splicing_df_filtered = splicing_df_all[(splicing_df_all['pipeline_interaction'] == True) &
    (splicing_df_all['splice_confidence'] > 0.60) & 
    (splicing_df_all['seqio_length'] < 3000)]
splicing_df_filtered.index = range(0,len(splicing_df_filtered))


## Add sequences to database (only have to do this once)
#for index in range(0,len(splicing_df_filtered)):
#    splicing_db.add((splicing_df_filtered.ix[index]['uniprot_accession'],
#                    splicing_df_filtered.ix[index]['splicing_accession'],),
#                    splicing_df_filtered.ix[index]['seqio_object'])


## Export sequences to files (only have to do this once)                    
#splicing_db.export_sequences()
#splicing_db.export_pickle()


## Make alignments and put them to a folder
#splicing_db.allign_in_tcoffee()
#splicing_db.export_pickle()

## Analyse alignments
#splicing_db.analyse_alignments()
#splicing_db.export_pickle()


class buildPD:
    def __init__(self):       
        self.splicing_db = {'uniprotID':[], 'splicingID':[], 'pFamID1':[], 'pFamID2':[], 'start1':[], 'end1':[], 'start2':[], 'end2':[], 'interactingAA':[],\
                            'domain_length':[], 'interface_length':[], 'pFam_occurance':[], 'mdfreq':[], 'mdfrac':[], 'mifreq':[], 'mifrac':[], 'in_3DID':[]}

    def add(self, uniprotID, splicingID, pFamID1, pFamID2, start1, end1, start2, end2, interactingAA, domain_length, interface_length, pFam_occurance, \
                    mdfreq, mdfrac, mifreq, mifrac, in_3DID):
        self.splicing_db['uniprotID'].append(uniprotID)
        self.splicing_db['splicingID'].append(splicingID)
        self.splicing_db['pFamID1'].append(pFamID1)
        self.splicing_db['pFamID2'].append(pFamID2)
        self.splicing_db['start1'].append(start1)
        self.splicing_db['end1'].append(end1)
        self.splicing_db['start2'].append(start2)
        self.splicing_db['end2'].append(end2)
        self.splicing_db['interactingAA'].append(interactingAA)
        self.splicing_db['domain_length'].append(domain_length)
        self.splicing_db['interface_length'].append(interface_length)
        self.splicing_db['pFam_occurance'].append(pFam_occurance)
        self.splicing_db['mdfreq'].append(mdfreq)
        self.splicing_db['mdfrac'].append(mdfrac)
        self.splicing_db['mifreq'].append(mifreq)
        self.splicing_db['mifrac'].append(mifrac)
        self.splicing_db['in_3DID'].append(in_3DID)

    def export_as_db(self, filename):
        print 'saving the splicing database'
        f = open(filename, 'wb')    
        pickle.dump(self.splicing_db, f)
        f.close()
    
    def export_as_df(self, filename):
        print 'saving the splicing dataframe'
        f = open(filename, 'wb')    
        pickle.dump(pd.DataFrame(self.splicing_db), f)
        f.close()

handle = open('splicing_sequence_alignmnets_by_splicingID.pickle', 'r')
splicing_db = pickle.load(handle)
handle.close()

#output_handle_domain = open('/home/alexey/subgroup meetings/splicing/131009/spliced-domain-list.txt', 'w')
#output_csv_domain = csv.writer(output_handle_domain, delimiter= '\t')
#output_handle_interface = open('/home/alexey/subgroup meetings/splicing/131009/spliced-interface-list.txt', 'w')
#output_csv_interface = csv.writer(output_handle_interface, delimiter= '\t')
#
#mutation_domain_frequency = []
#mutation_domain_frequency_usable = []
#mutation_domain_fraction = []
#mutation_domain_fraction_usable = []
#mutation_interface_frequency = []
#mutation_interface_frequency_usable = []
#mutation_interface_fraction = []
#mutation_interface_fraction_usable = []

##unique_proteins = set()
##unique_interactions = set()
#
#no_3DID = 0
#no_3DID_reverse = 0
#num_dimers = 0
#both_no_3DID_and_dimers = 0
                
splicing_df_nr_mutated = buildPD()

for dict_key in splicing_db:
    pFam_list = []
    uniprotID = dict_key[0]
    splicingID = dict_key[1]
    unique_interactions = set()
    for interaction in pipeline_interaction_database[uniprotID]:
        pFamID1 = interaction[2]
        pFamID2 = interaction[7]
        
        start1 = int(interaction[3].strip().split('-')[0])
        end1 = int(interaction[3].strip().split('-')[-1])
        start2 = int(interaction[8].strip().split('-')[0])
        end2 = int(interaction[8].strip().split('-')[-1])
        interactingAA = [int(mutAA) for mutAA in interaction[4].strip().split(',') if mutAA != '']
        domain_length = end1 - start1 - 1 # not counting the first and last AA in interaction domain
        interface_length = len(interactingAA)
        
        num_mutations_in_domain = 0
        num_mutations_in_interface = 0
        for mutation in splicing_db[dict_key][3]:
            if mutation > start1 and mutation < end1: # have to confirm if the range is inclusive or not
                num_mutations_in_domain += 1
            if mutation in interactingAA:
                num_mutations_in_interface += 1
                    
        mdfreq = num_mutations_in_domain
        mdfrac = float(num_mutations_in_domain)/domain_length if domain_length > 0 else 0
        mifreq = num_mutations_in_interface
        mifrac = float(num_mutations_in_interface)/interface_length if interface_length > 0 else 0
        
        if pipeline_3DID_database.has_key((pFamID1, pFamID2)):
            in_3DID = 'forward'
        elif pipeline_3DID_database.has_key((pFamID2, pFamID1)):
            in_3DID = 'reverse'
        else:
            in_3DID = 'absent'
        
        if tuple(interactingAA) not in unique_interactions and interface_length > 0:
            unique_interactions.add(tuple(interactingAA))
            pFam_list.append(pFamID1)
            pFam_occurance = pFam_list.count(pFamID1)
            splicing_df_nr_mutated.add(uniprotID, splicingID, pFamID1, pFamID2, start1, end1, start2, end2, interactingAA, domain_length, interface_length, pFam_occurance, \
                    mdfreq, mdfrac, mifreq, mifrac, in_3DID)
            if pFamID1 == 'dsrm':
                print pFam_list
                print pFam_occurance
                print interactingAA
                print


splicing_df_nr_mutated.export_as_df('splicing_interactions_nr_df.pickle')

handle = open('splicing_interactions_df.pickle', 'r')
splicing_df = pickle.load(handle)
handle.close()

handle = open('splicing_interactions_nr_df.pickle', 'r')
splicing_nr_df = pickle.load(handle)
handle.close()


#output_handle_domain.close()
#output_handle_interface.close()
                
#mdfreq = [mut_freq for mut_freq in mutation_domain_frequency_usable if mut_freq > 0]
#mdfrac = [mut_frac for mut_frac in mutation_domain_fraction_usable if mut_frac > 0]
#
#mifreq = [mut_freq for mut_freq in mutation_interface_frequency_usable if mut_freq > 0]
#mifrac = [mut_frac for mut_frac in mutation_interface_fraction_usable if mut_frac > 0]


# plotting
if False:
    f = pl.figure(figsize=(8,6), dpi=150)
    
    pl.subplot(2,2,1)
    n12, bins12, patches12 = pl.hist(splicing_df[splicing_df['mdfreq'] > 0 ]['mdfreq'], np.arange(0, 310, 15))
    pl.xlim(0, 305)
    pl.xlabel('No. of AA of pFam domain removed')
    pl.ylabel('No. of splicing events')
    
    pl.subplot(2,2,2)
    n22, bins22, patches22 = pl.hist(splicing_df[splicing_df['mdfreq'] > 0 ]['mdfrac'], np.arange(0, 1.1, 0.05))
    pl.xlim(0, 1.05)
    pl.xlabel('Fraction of pFam domain removed')
    
    pl.subplot(2,2,3)
    n12, bins12, patches12 = pl.hist(splicing_df[splicing_df['mifreq'] > 0 ]['mifreq'], np.arange(0, 310, 15))
    pl.xlim(0, 305)
    pl.xlabel('No. of interacting AA removed')
    pl.ylabel('No. of splicing events')
    
    pl.subplot(2,2,4)
    n22, bins22, patches22 = pl.hist(splicing_df[splicing_df['mifrac'] > 0 ]['mifrac'], np.arange(0, 1.1, 0.05))
    pl.xlim(0, 1.05)
    pl.xlabel('Fraction of interacting AA removed')
    
    pl.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.98, wspace=0.15, hspace=0.30)
    pl.show()













#


         
#        splicing_df[(splicing_df['pipeline_interaction'] == True)]['seqio_length'].hist()

#        return

#        splicing_db.allign_in_tcoffee([hit for hit in set(high_confidence_hits)])
#        splicing_db.allign_in_tcoffee(['all', ])

#        ## Part for multiprocessing ##
#        # see: http://doughellmann.com/2009/04/pymotw-multiprocessing-part-1.html
#        # Establish communication queues
#        tasks = multiprocessing.JoinableQueue()
#        results = multiprocessing.Queue()
#        
#        # create a logger instance
#        # I started to play with the logger a little bit but didn't have the time
#        # to fully implement it to be really usefull. It works, one just has
#        # to set the desired logging with the information where needed
#        log = MultiProcessingLog('results/multiprocessingLog.log', mode='a', maxsize=999, rotate=5)
#        
#        # create the pool to control the number of t_coffee instances
#        # that are allowed to run in parallel
#        pool = ActivePool()
#        s = multiprocessing.Semaphore(1) # int for number of parallel runs of tCoffee
#        
#        # Start consumers
#        print 'Creating %d consumers' % num_consumers
#        proc_name = [ 'Consumer-' + str(i) for i in range(1, num_consumers + 1) ]
#        consumers = [ Consumer(proc_name[i-1], tasks, results, runTime, pool, s, DEBUG, outputPath, log)
#                      for i in range(1, num_consumers + 1) ]
#
#        for w in consumers:
#            w.start()
#    
#        # check if a skip file is present
#        # I normally didn't use it but the idea is that if you run the pipeline
#        # and it reached the maximum runtime, you could easily re-run and skip
#        # the already calculated entries.
#        if os.path.isfile('processed.log'):
#            skipFile = open('processed.log', 'r')
#            skip = list()
#            for line in skipFile:
#                skip.append(line.strip())
#            skipFile.close()
#            SKIP = True
#        else:
#            SKIP = False
#        
#        DO_MULTIPROCESSING = True
#         
#
#mutation_uniprot = True # True if you are given a uniprot sequence and AA to be mutated, rather than a crystal structure
#outputPath = 'results/'
#savePDB = 'results/pdbFiles/'
#tmpPath = '/tmp/pipeline-splicing/'
#pdbPath = '/home/alexey/working/lib/pdb_database/structures/divided/pdb/'
#matrix = MatrixInfo.blosum80
#gap_start = -16
#gap_extend = -4
#modeller_runs = 1 # for modeller
#buildModel_runs = 1 # for FoldX
#foldX_WATER = '-IGNORE'
#interaction_database    = getInteractions('')
#threeDID_database       = get3DID('', '')
#get_uniprot_sequence    = getUniprotSequence('')
#pdb_resolution_database = pdb_resolution('/home/alexey/working/pipeline-splicing/databases/pdbtosp.txt')
#core_template_database  = core_Templates('')
#
#print 'got this far'
#multiTask = multi.Task(mutation_uniprot,
#                        'Q9Y297', 	# 'O14733', # uniprotKB, 
#                        'F514Y',    # 'P41C', # mutation,
#                        savePDB, 
#                        tmpPath, 
#                        outputPath,
#                        pdbPath,
#                        matrix, 
#                        gap_start, 
#                        gap_extend, 
#                        modeller_runs,
#                        buildModel_runs,
#                        foldX_WATER,
#                        interaction_database,
#                        threeDID_database,
#                        get_uniprot_sequence,
#                        pdb_resolution_database,
#                        core_template_database
#                        )
#
#print multiTask
#
#
#try:
#    print multiTask.findTemplate()
#except:
#    print 'something went wrong'
#    raise
#
#print 'done part one'
#
#try:
#    print multiTask()
#except:
#    print 'something went wrong'
#    raise




