# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 17:13:11 2015

@author: alexey
"""

#%% Imports ommon to all test files
from __future__ import print_function

import sys
biopython_develop_paths = [
    '/home/kimlab1/strokach/working/packages/biopython',
    '/home/kimlab1/strokach/working/packages/biopython/build/lib.linux-x86_64-3.4',
]
for path in biopython_develop_paths:
    if path not in sys.path:
        sys.path.insert(0, path)    

import os
import os.path as op

from common import constants

from elaspic import conf, pipeline, helper_functions as hf, pdb_template
    
try:
    code_path = op.abspath(op.join(os.path.dirname(__file__), '..'))
except:
    code_path = op.abspath(op.dirname(os.getcwd()))


#%%
%load_ext autoreload
%autoreload 2



#%%
PDB_DATABASE_PATH = '/home/kimlab1/database_data/pdb/data/data/structures/divided/pdb/'
TEMP_PATH = conf.get_temp_path()
CONFIG_FILE = op.join(code_path, 'config', 'config_file.ini')
LOGGER = hf.get_logger()

pipeline = pipeline.Pipeline(CONFIG_FILE)



#%%
def test_1():
    pdb_template.get_pdb()


def test_2():
    pdb_template_options = {
        'pdb_id': '4f3l',
        'chain_ids': ['B', 'A'],
        'domain_defs': ['71:241', '42:186'],
        'pdb_path': PDB_DATABASE_PATH,
        'output_path': TEMP_PATH,
        'tmp_path': TEMP_PATH,
        'logger': LOGGER
    }
    pdb = pdb_template.PDBTemplate(**pdb_template_options)
    pdb.extract() # chain becomes domain after running ``pdb.extract()``
    pdb.save_structure()


def test_3():
    """
    """
    pdb_template_options = {
        'pdb_id': '3vkh',
        'chain_ids': ['B'],
        'domain_defs': ['3212:2662'],
        'pdb_path': PDB_DATABASE_PATH,
        'output_path': TEMP_PATH,
        'tmp_path': TEMP_PATH,
        'logger': LOGGER
    }
    pdb = pdb_template.PDBTemplate(**pdb_template_options)
    pdb.extract() # chain becomes domain after running ``pdb.extract()``
    pdb.save_structure()

    
#%%
test_3()


#%%

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
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
fixer = PDBFixer(filename='myfile.pdb')
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.removeHeterogens(True)
fixer.addMissingHydrogens(7.0)
fixer.addSolvent(fixer.topology.getUnitCellDimensions())
PDBFile.writeFile(fixer.topology, fixer.positions, open('output.pdb', 'w'))

