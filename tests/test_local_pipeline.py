import logging
import pytest
from elaspic import conf
from conftest import DEFAULT_LOCAL_CONFIG as CONFIG_FILE
import helper_fns

logger = logging.getLogger(__name__)


# %% Constants
if hasattr(pytest, "config"):
    QUICK = pytest.config.getoption('--quick')
    CONFIG_FILE = (
        pytest.config.getoption('--config-file') or CONFIG_FILE
    )
else:
    QUICK = False

conf.read_configuration_file(CONFIG_FILE, unique_temp_dir=None)

logger.debug('Running quick: {}'.format(QUICK))
logger.debug('Config file: {}'.format(CONFIG_FILE))


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


def test_pdb_mutation_pipeline(pdb_id):
    return helper_fns.run_pdb_mutation_pipeline(pdb_id, pdb_mutatations)


def test_sequence_mutation_pipeline(pdb_id_sequence):
    return helper_fns.run_sequence_mutation_pipeline(pdb_id_sequence, sequence_mutations,)


# %%
if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-svx', '--quick'])
