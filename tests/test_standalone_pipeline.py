""".

.. todo::

    - Test that the standalone pipeline folder is relocatable.
"""
import os
import os.path as op
import logging
import shutil
import pytest
from elaspic import conf
import helper_fns

logger = logging.getLogger(__name__)

# Constants
CONFIG_FILE = op.join(op.dirname(__file__), 'standalone_pipeline_local.ini')

if hasattr(pytest, "config"):
    QUICK = pytest.config.getoption('--quick')
    CONFIG_FILE = pytest.config.getoption('--config-file') or CONFIG_FILE
else:
    QUICK = False

logger.debug('Running quick: {}'.format(QUICK))
logger.debug('Config file: {}'.format(CONFIG_FILE))


# Source of good PDB stuctures: http://www.rcsb.org/pdb/101/motm_archive.do
pdb_mutatations = {
    # This one was giving me trouble
    '1S1Q': {
        'A': [
            'V43A',
            'F44A',
            'N45A',
            'D46A',
            'F88A',
        ],
    },
    # test_1; only one chain
    '3M7R': {
        'A': [
            'I271A',  # core
            'Y401A',  # core
        ]
    },
    # test_2; this one has three symetric chains
    '1THJ': {
        # Chain B interacts with A & C
        # All of these mutations should had one core and two interface predictions
        # Mutation type 1
        'B': [
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
        # Mutation type 2
        'B': [
            'P37A',
            'E38A',
            'P55A',
            'M56A',
            'G78A',
            'Q119A',
            'Q137A',
            'P153A',
            'R154A',
            'G171A',
        ],
        # Mutation type 3
        '2': [
            'P37A',
            'E38A',
            'P55A',
            'M56A',
            'G78A',
            'Q119A',
            'Q137A',
            'P153A',
            'R154A',
            'G171A',
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
    ('2FOY', 'P23280'): {
        #  Mutation type 1
        'A': [
            'Q15A',
            'N24A',
        ],
        #  Mutation type 2
        'A': [
            'H34A',
            'G43A',
        ],
        #  Mutation type 3
        '1': [
            'H34A',
            'G43A',
        ],
    },
    ('2Z5Y', 'Q5NU32'): {
        # Mutation type 1 / 2
        'A': [
            'H12A',
            'M13A',
        ],
        #  Mutation type 3
        '1': [
            'H12A',
            'M13A',
            'K30A',
            'L31A',
        ]
    },
}


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


# Fixtures
@pytest.fixture(scope='session', params=list(pdb_mutatations.keys()))
def pdb_id(request):
    return request.param


@pytest.fixture(scope='session', params=list(sequence_mutations.keys()))
def pdb_id_sequence(request):
    return request.param


def test_pdb_mutation_pipeline(pdb_id):
    # Canonical folder
    unique_temp_dir = op.join(op.splitext(__file__)[0], pdb_id, '.elaspic')
    os.makedirs(unique_temp_dir, exist_ok=True)
    conf.read_configuration_file(CONFIG_FILE, DEFAULT={'unique_temp_dir': unique_temp_dir})
    os.chdir(unique_temp_dir)
    helper_fns.run_pdb_mutation_pipeline(pdb_id, pdb_mutatations)

    # Make sure it works even if we copy folder
    unique_temp_dir_old = unique_temp_dir

    # Precalculated sequence
    logger.debug('\n\n' + '*' * 80 + 'Precalculated sequence\n')
    working_dir = op.join(op.splitext(__file__)[0], pdb_id + '_precalc_sequence')
    unique_temp_dir = op.join(working_dir, '.elaspic')
    os.makedirs(unique_temp_dir, exist_ok=True)
    try:
        shutil.copy2(
            op.join(unique_temp_dir_old, 'sequence.json'),
            op.join(unique_temp_dir, 'sequence.json'))
        shutil.copytree(
            op.join(unique_temp_dir_old, 'sequence'),
            op.join(unique_temp_dir, 'sequence'))
        conf.read_configuration_file(CONFIG_FILE, DEFAULT={'unique_temp_dir': unique_temp_dir})
        os.chdir(unique_temp_dir)
        helper_fns.run_pdb_mutation_pipeline(pdb_id, pdb_mutatations)
    except:
        raise
    finally:
        shutil.rmtree(unique_temp_dir)

    # Precalculated model
    logger.debug('\n\n' + '*' * 80 + 'Precalculated model\n')
    working_dir = op.join(op.splitext(__file__)[0], pdb_id + '_precalc_model')
    unique_temp_dir = op.join(working_dir, '.elaspic')
    os.makedirs(unique_temp_dir, exist_ok=True)
    try:
        shutil.copy2(
            op.join(unique_temp_dir_old, 'sequence.json'),
            op.join(unique_temp_dir, 'sequence.json'))
        shutil.copytree(
            op.join(unique_temp_dir_old, 'sequence'),
            op.join(unique_temp_dir, 'sequence'))
        shutil.copy2(
            op.join(unique_temp_dir_old, 'model.json'),
            op.join(unique_temp_dir, 'model.json'))
        shutil.copytree(
            op.join(unique_temp_dir_old, 'model'),
            op.join(unique_temp_dir, 'model'))
        conf.read_configuration_file(CONFIG_FILE, DEFAULT={'unique_temp_dir': unique_temp_dir})
        os.chdir(unique_temp_dir)
        helper_fns.run_pdb_mutation_pipeline(pdb_id, pdb_mutatations)
    except:
        raise
    finally:
        shutil.rmtree(unique_temp_dir)


def test_sequence_mutation_pipeline(pdb_id_sequence):
    unique_temp_dir = op.join(op.splitext(__file__)[0], '.'.join(pdb_id_sequence), '.elaspic')
    os.makedirs(unique_temp_dir, exist_ok=True)
    conf.read_configuration_file(CONFIG_FILE, DEFAULT={'unique_temp_dir': unique_temp_dir})
    os.chdir(unique_temp_dir)
    return helper_fns.run_sequence_mutation_pipeline(pdb_id_sequence, sequence_mutations,)
