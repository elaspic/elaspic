import os.path as op
import shutil
import tempfile

import pytest

import elaspic
import elaspic.tools
from kmtools import py_tools, sequence_tools

logger = py_tools.get_logger(__name__)

ANALYZERS = [elaspic.tools.Conservation, elaspic.tools.Provean, ]
SEQUENCE_MUTATIONS = [
    ('UPI0000047181', 'M0A:Q9A:A10G:T19A:F20A:I99A:R439A:P440A'),
    ('Q8IVG9', 'M0A:A1G:P2A:R3A:G4A:F5A:T12A:K20A:R21A:R22A:A23G'),
]

elaspic.CONFIGS['sequence_dir'] = op.join(tempfile.mkdtemp(prefix='sequence-'))
elaspic.CONFIGS['n_cores'] = 1
elaspic.CONFIGS['blast_db_dir'] = '/home/kimlab1/database_data/blast/db'


@pytest.mark.parametrize("sequence_id, mutations", SEQUENCE_MUTATIONS)
@pytest.mark.parametrize("Analyzer", [elaspic.tools.Conservation])
def test_1(Analyzer, sequence_id, mutations):
    # Prepare structure
    s = sequence_tools.fetch_sequence(sequence_id)
    # Prepare analyzer
    analyzer = Analyzer(s)
    analyzer.build()
    assert analyzer.done
    for mutation in mutations.split(':'):
        logger.debug('{}-{}', sequence_id, mutation)
        result = analyzer.analyze(mutation)
        logger.debug("{}\n", result)


@pytest.mark.parametrize("sequence_id, mutations", SEQUENCE_MUTATIONS)
@pytest.mark.parametrize("Analyzer", [elaspic.tools.Provean])
def test_2(Analyzer, sequence_id, mutations):
    # Prepare structure
    s = sequence_tools.fetch_sequence(sequence_id)
    # Prepare analyzer
    analyzer = Analyzer(s)
    # Set precalculated data
    supporting_set_file = op.abspath(op.join(
        op.dirname(__file__),
        'test_{}'.format(analyzer.__class__.__name__.lower()),
        '{}_provean_supset'.format(sequence_id)))
    shutil.copy(supporting_set_file, analyzer.supporting_set_file)
    analyzer.result = {
        'supporting_set_file': analyzer.supporting_set_file,
        'supporting_set_length': analyzer.supporting_set_length,
    }
    assert analyzer.done
    # Build
    analyzer.build()
    assert analyzer.done
    for mutation in mutations.split(':'):
        logger.debug('{}-{}', sequence_id, mutation)
        result = analyzer.analyze(mutation)
        logger.debug("{}\n", result)
