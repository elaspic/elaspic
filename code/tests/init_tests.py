# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 22:35:03 2015

@author: Alexey Strokach
"""
#%% Imports common to all test files
import os
import os.path as op
import sys

try:
    code_path = op.join(op.dirname(op.abspath(__file__)), '..')
except:
    code_path = op.join(op.dirname(os.getcwd()), '..')

print('code_path: {}'.format(code_path))
if code_path not in sys.path:
    sys.path.insert(0, code_path)

try:
    from elaspic import Pipeline
    import helper_functions as hf
except ImportError:
    del sys.modules['elaspic']
    from elaspic import Pipeline
    import helper_functions as hf

default_config_file = os.path.join(code_path, '../config/config_file.ini')
temp_path = hf.get_temp_path()


