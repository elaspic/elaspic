# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 13:39:23 2015

@author: strokach
"""
import os.path as op
import logging
import pytest
from elaspic import (
    conf, sequence, structure_tools, local_pipeline
)
from conftest import TESTS_BASE_DIR

logger = logging.getLogger(__name__)


# %% Constants
QUICK = (hasattr(pytest, "config") and pytest.config.getoption('--quick'))
logger.debug("\nRunning quickly: {}".format(QUICK))


# %%
# Source of good PDB stuctures: http://www.rcsb.org/pdb/101/motm_archive.do
pdb_mutatations = {
    # test_1; only one chain
    '3M7R': {
        'A': [
             'I271A',  # core
             'Y401A',  # core
        ]
    },
    # test_2; this one has three symetric chains
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


sequence_mutations = {
    ('2Z5Y', 'Q5NU32'): {
        '1': [
            'H12A',
            'M13A',
            'F14A',
            'D15A',
            'V16A',
            'V17A',
            'V18A',
            'I19A',
            'G20A',
            'G21A',
            'G22A',
            'I23A',
            'S24A',
            'G25A',
            'L26A',
            'S27A',
            'A28G',
            'A29G',
            'K30A',
            'L31A',
        ]
    },
}


# %%
if QUICK:
    # pdb_mutatations
    for pdb_id in pdb_mutatations:
        for chain_id in pdb_mutatations[pdb_id]:
            mutations = pdb_mutatations[pdb_id][chain_id]
            break
        break
    pdb_mutatations = {pdb_id: {chain_id: [mutations[0]]}}
    # sequence_mutations
    for key in sequence_mutations:
        for chain_id in sequence_mutations[key]:
            mutations = sequence_mutations[key][chain_id]
            break
        break
    sequence_mutations = {key: {chain_id: [mutations[0]]}}


# %% Fixtures
@pytest.fixture(scope='session', params=list(pdb_mutatations.keys()))
def pdb_id(request):
    return request.param


@pytest.fixture(scope='session', params=list(sequence_mutations.keys()))
def pdb_id_sequence(request):
    return request.param


# %%
def test_pdb_mutation_pipeline(pdb_id):
    working_dir = None  # can set to something if don't want to rerun entire pipeline
    config_file = op.join(TESTS_BASE_DIR, 'test_local_pipeline.ini')
    conf.read_configuration_file(config_file, unique_temp_dir=working_dir)
    configs = conf.Configs()
    pdb_file = structure_tools.download_pdb_file(pdb_id, configs['unique_temp_dir'])
    for chain_id in pdb_mutatations[pdb_id]:
        for mutation in pdb_mutatations[pdb_id][chain_id]:
            mutation_pdb = '{}_{}'.format(chain_id, mutation)
            lp = local_pipeline.LocalPipeline(
                pdb_file, mutations=mutation_pdb)
            lp.run_all_sequences()
            lp.run_all_models()
            lp.run_all_mutations()


def test_sequence_mutation_pipeline(pdb_id_sequence):
    pdb_id, sequence_id = pdb_id_sequence
    working_dir = None  # can set to something if don't want to rerun entire pipeline
    config_file = op.join(TESTS_BASE_DIR, 'test_local_pipeline.ini')
    conf.read_configuration_file(config_file, unique_temp_dir=working_dir)
    configs = conf.Configs()
    pdb_file = structure_tools.download_pdb_file(pdb_id, configs['unique_temp_dir'])
    sequence_file = sequence.download_uniport_sequence(sequence_id, configs['unique_temp_dir'])
    for chain_pos in sequence_mutations[pdb_id_sequence]:
        for mutation in sequence_mutations[pdb_id_sequence][chain_pos]:
            mutation_sequence = '{}_{}'.format(chain_pos, mutation)
            lp = local_pipeline.LocalPipeline(
                pdb_file, sequence_file, mutations=mutation_sequence)
            lp.run_all_sequences()
            lp.run_all_models()
            lp.run_all_mutations()


# %%
if __name__ == '__main__':
    pytest.main([__file__, '-vsx', '--quick'])
