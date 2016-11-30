"""ELASPIC.
"""
__all__ = [
    'conf',
    'exc',
    'tools',
]
__version__ = "0.1.40"

import os
import os.path as op

from . import *
from ._pipeline import _Pipeline
from .model import Model
from .predictor import Predictor
from .sequence import Sequence
from .standalone_pipeline import StandalonePipeline

BASE_DIR = op.abspath(op.dirname(__file__))
DATA_DIR = op.join(BASE_DIR, 'predictor')
CACHE_DIR = op.join(BASE_DIR, 'predictor')
os.makedirs(CACHE_DIR, exist_ok=True)

CONFIGS = {}
