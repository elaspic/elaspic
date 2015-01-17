# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 13:57:36 2014

@author: alexey
"""
#%%
import os
import sys

import __main__
if hasattr(__main__, '__file__'):
    code_path = os.path.dirname(__file__)
else:
    code_path = os.path.dirname(os.getcwd())

if code_path not in sys.path:
    sys.path.insert(0, code_path)


import elaspic
import helper_functions as hf

default_config_file = os.path.join(code_path, '../config/config_file.ini')
temp_path = hf.get_temp_path()
pipeline = elaspic.Pipeline(default_config_file)




def test_1():
    """Test cases where there exists a model but the chains are not interacting
    """
    uniprot_id = 'P15056'
    mutations = 'R575F'
    run_type = 3
    n_cores = 1

    uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutations, run_type, n_cores)

    print('Done!')



def test_2():
    """Test for error: "Exception: Expected and actual FoldX amino acids do not match!"
    """
    uniprot_id = 'O14745'
    mutations = 'D57G,V60G,Y38H,D183V,I179V,P152S,P187R'
    run_type = 3
    n_cores = 1

    uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutations, run_type, n_cores)

    print('Done!')



def test_3():
    """Test for error: "Exception: Expected and actual FoldX amino acids do not match!"
    """
    uniprot_id = 'O14594'
    mutations = 'E1167K,F1177L,A69T,D72Y,D201G,D1155N,A339T,D1183E'
    run_type = 3
    n_cores = 1

    uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutations, run_type, n_cores)

    print('Done!')















#%%



