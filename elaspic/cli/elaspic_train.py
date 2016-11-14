"""ELASPIC TRAIN
"""
import os.path as op
import json
import pickle
import logging
import pandas as pd
from elaspic import DATA_DIR, CACHE_DIR, elaspic_predictor


logger = logging.getLogger(__name__)


def elaspic_train(args):
    """Train ELASPIC core and interface predictors."""
    _train_predictor('core')
    _train_predictor('interface')


def _train_predictor(prefix):
    """Train ELASPIC core (`prefix` == 'core') OR interface (`prefix` == interface) predictor."""
    logger.debug("_train_predictor(%s)", prefix)
    assert prefix in ['core', 'interface']
    # Load data
    data = pd.read_csv(
        op.join(DATA_DIR, '{}_training_set.tsv.gz'.format(prefix)),
        sep='\t', low_memory=False)
    with open(op.join(DATA_DIR, '{}_features.json'.format(prefix)), 'rt') as fh:
        features = json.load(fh)
    with open(op.join(DATA_DIR, '{}_options.json'.format(prefix)), 'rt') as fh:
        options = json.load(fh)
    # Train
    core_predictor = elaspic_predictor.get_final_predictor(data, features, options)
    # Save
    with open(op.join(CACHE_DIR, '{}_clf.pickle'.format(prefix)), 'wb') as fh:
        pickle.dump(core_predictor, fh, pickle.HIGHEST_PROTOCOL)


def configure_train_parser(sub_parsers):
    help = "Train the ELASPIC classifiers"
    description = help + """
"""
    example = """
Examples:

    elaspic train

"""
    parser = sub_parsers.add_parser(
        'train',
        help=help,
        description=description,
        epilog=example,
    )
    parser.set_defaults(func=elaspic_train)
