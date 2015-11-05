# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 13:39:23 2015

@author: strokach
"""
import os.path as op
import logging
import pytest

from elaspic import (
    conf, sequence, model, structure_tools, structure_analysis, pipeline, call_tcoffee, helper,
    predictor, local_pipeline
)

logger = logging.getLogger(__name__)
basedir = op.dirname(__file__)


# %%
# Source of good PDB stuctures: http://www.rcsb.org/pdb/101/motm_archive.do
pdb_mutatations = {
    # only one chain
    '3M7R': {
        'A': [
             'I271A',  # core
             'Y401A',  # core
        ]
    },
    # this one has three symetric chains
    '1THJ': {
        'B': [
            # interface with A & C
            # all of these should had one core and two interface predictions
            'P36A',
            'E37A',
            'P54A',
            'M55A',
            'G77A',
            'Q118A',
            'Q136A',
            'P152A',
            'R153A',
            'G170A',
        ],
    },
    # this one has two chains and DNA in it...
    '3OS0': {
        'A': [
            'K36A',  # surface, no interface
            'V58A',  # core
            'V89A',  # core
            'P247A',  # interface with B_S175
            'Q250A',  # interface with B_S175
            'L251A',  # interface with B_S175
            'N275A',  # interface with B_I178
        ],
        'B': [
            'S175A',  # interface with A_P247, A_Q250, A_L251
            'I178A',  # interface with A_N275
        ]
    },
}


working_dirs = {
    '3M7R': '/tmp/elaspic/x5szda5d',
    '1THJ': '/tmp/elaspic/7siaev7u',
    '3OS0': '/tmp/elaspic/sz9wlfc2',
}


# %% Fixtures
@pytest.fixture(scope='session', params=list(pdb_mutatations.keys()))
def pdb_id(request):
    return request.param


# %%
def test_local_pipeline(pdb_id):
    working_dir = working_dirs.get(pdb_id, None)
    config_file = op.join(basedir, 'test_local_pipeline.ini')
    conf.read_configuration_file(config_file, unique_temp_dir=working_dir)
    configs = conf.Configs()
    pdb_file = structure_tools.download_pdb_file(pdb_id, configs['unique_temp_dir'])
    for chain_id in pdb_mutatations[pdb_id]:
        for mutation in pdb_mutatations[pdb_id][chain_id]:
            mutation_pdb = '{}_{}'.format(chain_id, mutation)
            lp = local_pipeline.LocalPipeline(pdb_file, mutations=mutation_pdb)
            lp.run_all_models()
            lp.run_all_mutations()


# %%
if __name__ == '__main__':
    import pytest
    pytest.main(['test_local_pipeline.py', '-svx'])
