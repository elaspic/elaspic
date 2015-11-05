# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 14:28:57 2015

@author: strokach
"""

import os
import os.path as op
import shutil
import urllib.request
import pytest

def pytest_runtest_makereport(item, call):
    if "incremental" in item.keywords:
        if call.excinfo is not None:
            parent = item.parent
            parent._previousfailed = item

def pytest_runtest_setup(item):
    if "incremental" in item.keywords:
        previousfailed = getattr(item.parent, "_previousfailed", None)
        if previousfailed is not None:
            pytest.xfail("previous test failed (%s)" %previousfailed.name)


#%% Directories
try:
    TESTS_BASE_DIR = op.dirname(__file__)
except NameError:
    TESTS_BASE_DIR = os.getcwd()
assert op.split(TESTS_BASE_DIR)[-1] == 'tests'

TESTS_DATA_DIR = op.join(TESTS_BASE_DIR, 'data')


#%%
import logging
import logging.config
logger = logging.getLogger(__name__)

LOGGING_CONFIGS = {
    'version': 1,              
    'disable_existing_loggers': False,  # this fixes the problem

    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        },
        'clean': {
            'format': '%(message)s',
        },
    },
    'handlers': {
        'default': {
            'level':'DEBUG',    
            'class':'logging.StreamHandler',
            'formatter':'clean',
        },  
    },
    'loggers': {
        '': {                  
            'handlers': ['default'],        
            'level': 'DEBUG',  
            'propagate': True  
        }
    }
}
logging.config.dictConfig(LOGGING_CONFIGS)



#%% Functions
PDB_URL = 'http://www.rcsb.org/pdb/files/{}.pdb'

def get_pdb(pdb_id, input_folder, output_folder, use_remote=False):
    """Move PDB structure to the local working directory.
    """
    input_pdb_filename = op.join(input_folder, pdb_id + '.pdb')
    output_pdb_filename = op.join(output_folder, pdb_id + '.pdb')

    # If the PDB already exists, do nothing...
    if op.isfile(output_pdb_filename):
        message = (
            'PDB file already exists in the temp folder ({}). No need to move or download.'
            .format(output_pdb_filename)
        )
        logger.debug(message)
        return output_pdb_filename

    # Look for PDB file in the same folder
    if not op.isfile(input_pdb_filename):
        if use_remote:
            logger.info('Downloading PDB from the web to {}...'.format(input_pdb_filename))
            response = urllib.request.urlopen(PDB_URL.format(pdb_id.upper()))
            with open(input_pdb_filename, 'wb') as ofh:
                ofh.write(response.read())
        else:
            raise Exception('No PDB input file found!')

    logger.info('Copying {} to {}...'.format(input_pdb_filename, output_pdb_filename))
    shutil.copy(input_pdb_filename, output_pdb_filename)
    return output_pdb_filename



def download_pdb(pdb_id, output_folder):
    """Move PDB structure to the local working directory.
    """
    output_pdb_filename = op.join(output_folder, pdb_id + '.pdb')

    # If the PDB already exists, do nothing...
    if op.isfile(output_pdb_filename):
        logger.debug('PDB file already exists in the folder: {}.'.format(output_folder))
        return output_pdb_filename

    # Download the PDB file from the internet...
    logger.info('Downloading PDB {}...'.format(pdb_id + '.pdb'))
    response = urllib.request.urlopen(PDB_URL.format(pdb_id))
    with open(output_pdb_filename, 'wb') as ofh:
        ofh.write(response.read())

    return output_pdb_filename



#%%
print('Conftest loaded!')
