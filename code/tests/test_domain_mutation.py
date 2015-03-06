# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 13:57:36 2014

@author: alexey
"""
from __future__ import print_function
#%% Common to all test files
from elaspic_tests import Pipeline, default_config_file
pipeline = Pipeline(default_config_file)


#%%
def test_1():
    """Test cases where there exists a model but the chains are not interacting
    """
    uniprot_id = 'P15056'
    mutations = 'R575F'
    results = pipeline(uniprot_id, mutations, 3)
    print('Done!')


def test_2():
    """Test for error: "Exception: Expected and actual FoldX amino acids do not match!"
    """
    uniprot_id = 'O14745'
    mutations = 'D57G,V60G,Y38H,D183V,I179V,P152S,P187R'
    results = pipeline(uniprot_id, mutations, 3)
    print('Done!')


def test_3():
    """Test for error: "Exception: Expected and actual FoldX amino acids do not match!"
    """
    uniprot_id = 'O14594'
    mutations = 'E1167K,F1177L,A69T,D72Y,D201G,D1155N,A339T,D1183E'
    results = pipeline(uniprot_id, mutations, 3)
    print('Done!')



#%%



