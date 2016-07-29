import os
import os.path as op
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

BASE_DIR = op.abspath(op.dirname(__file__))
DATA_DIR = op.join(BASE_DIR, 'data')
CACHE_DIR = op.join(BASE_DIR, 'data')
os.makedirs(CACHE_DIR, exist_ok=True)

# Imports have to happen after constants have been defined
from .elaspic_predictor import CorePredictor, InterfacePredictor

# Don't autoload submodules requiring database configuration
blacklist = [
    'elaspic_database',
    'elaspic_dtabase_tables',
    'database_pipeline',
]

__all__ = [
    op.splitext(f)[0]  # remove .py extension
    for f in os.listdir(BASE_DIR)  # list contents of current dir
    if not f.startswith('_') and
    ((op.isfile(op.join(BASE_DIR, f)) and f.endswith('.py')) or
     (op.isdir(op.join(BASE_DIR, f)) and op.isfile(op.join(BASE_DIR, f, '__init__.py')))) and
    not any(f.startswith(pre) for pre in blacklist)
]
from . import *
