__version__ = "0.1.46.dev0"

import logging
import os
import os.path as op

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

BASE_DIR = op.abspath(op.dirname(__file__))
DATA_DIR = op.join(BASE_DIR, "data")
CACHE_DIR = op.join(BASE_DIR, "data")
os.makedirs(CACHE_DIR, exist_ok=True)
