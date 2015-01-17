# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 17:13:11 2015

@author: alexey
"""

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