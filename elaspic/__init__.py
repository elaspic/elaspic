"""ELASPIC.
"""
# flake8: noqa
from .sequence import Sequence
from .model import Model

__all__ = [
    'conf',
    'exc',
    'pipeline',
    'tools',
]
from . import *

__version__ = "0.2.0"

import os
import os.path as op

BASE_DIR = op.abspath(op.dirname(__file__))
DATA_DIR = op.join(BASE_DIR, 'predictor', 'data')
CACHE_DIR = op.join(BASE_DIR, 'predictor', 'data')
os.makedirs(CACHE_DIR, exist_ok=True)

from .conf import CONFIGS

del os, op
