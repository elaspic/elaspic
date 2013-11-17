# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 13:24:52 2013

@author: alexey
"""

import modeller
import Bio
import pickle

f = open('pipeline_3DID_database.pickle','rb')
pipeline_3DID_database = pickle.load(f)
f.close()

f = open('pipeline_core_templates.pickle','rb')
pipeline_core_templates = pickle.load(f)
f.close()

f = open('pipeline_interaction_database.pickle','rb')
pipeline_interaction_database = pickle.load(f)
f.close()

f = open('pipeline_uniprot_sequences.pickle','rb')
pipeline_uniprot_sequences = pickle.load(f)
f.close()