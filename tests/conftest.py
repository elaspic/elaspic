# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 14:28:57 2015

@author: strokach
"""

import os
import os.path as op
import shutil
import pytest
import logging
import logging.config
from elaspic import structure_tools, sequence
logger = logging.getLogger(__name__)


# %% Directories
try:
    TESTS_BASE_DIR = op.dirname(__file__)
except NameError:
    TESTS_BASE_DIR = os.getcwd()
assert op.split(TESTS_BASE_DIR)[-1] == 'tests'

TESTS_DATA_DIR = op.join(TESTS_BASE_DIR, 'data')


# %% Logger
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
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'clean',
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


# %% Functions
PDB_URL = 'http://www.rcsb.org/pdb/files/{}.pdb'


def get_structure(pdb_id, input_folder, output_folder, use_remote=True):
    """Move PDB structure to the local working directory.
    """
    input_file = op.join(input_folder, pdb_id + '.pdb')
    output_file = op.join(output_folder, pdb_id + '.pdb')

    # If the PDB already exists, do nothing...
    if op.isfile(output_file):
        logger.debug('Structure file {} already exists!'.format(output_file))
        return output_file

    # Look for PDB file in the same folder
    if not op.isfile(input_file):
        if use_remote:
            input_file = structure_tools.download_pdb_file(pdb_id, input_folder)
        else:
            raise Exception('No PDB input file found!')

    logger.info('Copying {} to {}...'.format(input_file, output_file))
    shutil.copy(input_file, output_file)
    return output_file


def get_sequence(uniprot_id, input_dir, output_dir, use_remote=True):
    """Move PDB structure to the local working directory.
    """
    input_file = op.join(input_dir, uniprot_id + '.fasta')
    output_file = op.join(output_dir, uniprot_id + '.fasta')

    # If the PDB already exists, do nothing...
    if op.isfile(output_file):
        logger.debug('Sequence file {} already exists!'.format(output_file))
        return output_file

    # Look for PDB file in the same folder
    if not op.isfile(input_file):
        if use_remote:
            input_file = sequence.download_uniport_sequence(uniprot_id, input_dir)
        else:
            raise Exception('No PDB input file found!')

    logger.info('Copying {} to {}...'.format(input_file, output_file))
    shutil.copy(input_file, output_file)
    return output_file


# %% Pytest
def pytest_runtest_makereport(item, call):
    if "incremental" in item.keywords:
        if call.excinfo is not None:
            parent = item.parent
            parent._previousfailed = item


def pytest_runtest_setup(item):
    if "incremental" in item.keywords:
        previousfailed = getattr(item.parent, "_previousfailed", None)
        if previousfailed is not None:
            pytest.xfail("previous test failed (%s)" % previousfailed.name)


# %%
print('Conftest loaded!')
