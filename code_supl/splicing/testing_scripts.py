# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:42:41 2013

@author: alexey
"""

import pickle

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

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_uniprot_sequences.pickle', 'rb')
pipeline_uniprot_sequences = pickle.load(handle)
handle.close()

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_uniprot_sequences.pickle', 'rb')
pipeline_splicing_sequences = pickle.load(handle)
handle.close()


# Looking at Sabastians interaction database
#cd ~/working/pipeline-splicing/code
#run splicing_data_prep

counter = 0
while counter < 10:
    for uniprotKb in pipeline_interaction_database:
        unique_pFam = set()
        interacting_AA = set()
        for interaction in pipeline_interaction_database[uniprotKb]:
            if (interaction[2], interaction[3]) in unique_pFam:
                interfaceAA = [int(mutAA) for mutAA in interaction[4].strip().split(',') if mutAA != '']
                for mutation in interfaceAA:
                    if mutation not in interacting_AA:
                        counter += 1
                        print interacting_AA
                        print interaction
                        print
                        break
            if (interaction[2], interaction[3]) not in unique_pFam:
                interfaceAA = [int(mutAA) for mutAA in interaction[4].strip().split(',') if mutAA != '' and mutAA != 'Contact_surface']
                if len(interfaceAA) > 0:
                    unique_pFam.add((interaction[2], interaction[3]))
                    for mutation in interfaceAA:
                        interacting_AA.add(mutation)





import csv

df_keys = [key for key in splicing_df[(splicing_df['mifreq'] > 0) & (splicing_df['mifrac'] < 0.8)]]

with open('/home/alexey/subgroup meetings/splicing/131009/short_spliced.csv', 'a') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(['uniprotID', 'splicingID', 'pFamID1', 'pFamID2', 'mifrac'])
    used_pairs = set()
    for i, row in enumerate(splicing_df[(splicing_df['mifreq'] > 0) & (splicing_df['mifrac'] < 0.8)].values):
        if row[df_keys.index('pFamID1')] == row[df_keys.index('pFamID2')]:
            if tuple([row[df_keys.index('pFamID1')], row[df_keys.index('pFamID2')]]) in used_pairs:
                continue
            else:
                used_pairs.add(tuple((row[df_keys.index('pFamID1')], row[df_keys.index('pFamID2')])))
        csv_writer.writerow([row[df_keys.index('uniprotID')], row[df_keys.index('splicingID')],\
                                row[df_keys.index('pFamID1')], row[df_keys.index('pFamID2')], "%.2f" % row[df_keys.index('mifrac')]])

with open('/home/alexey/subgroup meetings/splicing/131009/long_spliced.csv', 'a') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(['uniprotID', 'splicingID','pFamID1', 'pFamID2', 'mifrac'])
    used_pairs = set()
    for i, row in enumerate(splicing_df[(splicing_df['mifreq'] > 0) & (splicing_df['mifrac'] >=0.8)].values):
        if row[df_keys.index('pFamID1')] == row[df_keys.index('pFamID2')]:
            if tuple([row[df_keys.index('pFamID1')], row[df_keys.index('pFamID2')]]) in used_pairs:
                continue
            else:
                used_pairs.add(tuple([row[df_keys.index('pFamID1')], row[df_keys.index('pFamID2')]]))
        csv_writer.writerow([row[df_keys.index('uniprotID')], row[df_keys.index('splicingID')],
                         row[df_keys.index('pFamID1')], row[df_keys.index('pFamID2')], "%.2f" % row[df_keys.index('mifrac')]])


 
#df_keys = [key for key in splicing_df[(splicing_df['mifreq'] > 0) & (splicing_df['mifrac'] < 0.8)]]
#df_keys.index('pFamID2')




from class_get_uniprot_template_core_and_interface import get_template_interface
self.tmpPath = '/tmp/pipeline-splicing/'
self.unique = 'Consumer-1'
pdbPath
savePDB
saveAlignments
pool
semaphore
get_uniprot_seq
get_interactions
get_3did_entries
get_resolution

new_template = get_template_interface():
    """
    Check if a mutation falls into an interface and tries to find the best
    structural template
    
    """
    
    def __init__(self, tmpPath, unique, pdbPath, savePDB, saveAlignments, pool, semaphore, get_uniprot_seq, get_interactions, get_3did_entries, get_resolution):
        """
        input
        
        tmpPath             type 'str'
        unique              type 'str'
        pdbPath             type 'str'
        savePDB             type 'str'
        saveAlignments      type 'str'
        pool                type class '__main__.ActivePool'
        semaphore           type class 'multiprocessing.synchronize.Semaphore'
        get_uniprot_seq     <__main__.getUniprotSequence instance>
        get_interactions    <__main__.getInteractions instance>
        get_3did_entries    <__main__.get3DID instance>
        get_resolution      <__main__.pdb_resolution instance>

 







