# -*- coding: utf-8 -*-
"""
"""
import os

#: Location where the ELASPIC package is installed
base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')

#: Location of the ELASPIC package code
code_path = os.path.join(base_path, 'elaspic/')

#: Location of the ELASPIC package tests
tests_path = os.path.join(base_path, 'tests/')

#: Location of the ELASPIC package executables
#: To be used only if `bin_path` is not set in the configuration file
bin_path = os.path.join(base_path, 'tests/')

