# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 17:13:11 2015

@author: alexey
"""


#%% Imports ommon to all test files
from __future__ import print_function
import os
import sys

try:
    code_path = os.path.join(os.path.dirname(__file__), '..')
except:
    code_path = os.path.dirname(os.getcwd())

if code_path not in sys.path:
    sys.path.insert(0, code_path)

try:
    from elaspic import Pipeline
    import helper_functions as hf
except ImportError:
    del sys.modules['elaspic']
    from elaspic import Pipeline
    import helper_functions as hf

default_config_file = os.path.join(code_path, '../config/config_file.ini')
temp_path = hf.get_temp_path()
pipeline = Pipeline(default_config_file)



#%%
import os
import sys

code_path = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, code_path)

import logging
from common import constants
import pdb_template
import helper_functions as hf

default_config_file = os.path.join(code_path, '../config/config_file.ini')


#%%
pdb_id = '4f3l'
pdb_chains = ['B', 'A']
pdb_domain_defs = ['71:241', '42:186']
tmp_path = hf.get_temp_path()


def test_1():
    pdb_template.get_pdb()


def test_2():
    pdb_template_options = {
        'pdb_path': constants.pdb_database_path,
        'pdb_id': pdb_id,
        'chain_ids': pdb_chains,
        'domain_defs': pdb_domain_defs,
        'output_path': tmp_path,
        'tmp_path': tmp_path,
        'logger': logging}
    pdb = pdb_template.PDBTemplate(**pdb_template_options)
    pdb.extract() # chain becomes domain after running ``pdb.extract()``
    pdb.save_structure()

pdb_template.get_pdb(pdb_id, constants.pdb_database_path, tmp_path=tmp_path, pdb_type='cif', use_external=True)




#%%
# Problem when modelling `Q16512` because the template pdb contains a phosphorylated serine

from elaspic import conf, helper_functions, pdb_template

pdb_id = '3qc4'
chain_ids = ['A', 'B']
domain_defs = ['76:357', '-5:359']

tmp_path = conf.get_temp_path(temp_path_suffix='test_pdb_template') + '/'
logger = helper_functions.get_logger()
pdb_path = '/home/kimlab1/database_data/pdb/data/data/structures/divided/pdb/'

pdb = pdb_template.PDBTemplate(pdb_path, pdb_id, chain_ids, domain_defs, tmp_path, tmp_path, logger)
pdb.extract()
pdb.save_sequences()
pdb.save_structure()


#%%

