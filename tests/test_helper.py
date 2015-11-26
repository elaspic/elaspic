# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 14:20:21 2015

@author: strokach
"""
import os

import pytest

from elaspic import helper


# %%
@pytest.fixture(scope='function')
def cleanup(request):
    def fin():
        try:
            os.remove('deleteme.txt')
        except FileNotFoundError:
            pass
    request.addfinalizer(fin)


def _worker(x):
    with helper.open_exclusively('deleteme.txt') as ofh:
        ofh.write('\t'.join(str(x) for i in range(10)) + '\n')


def test_open_exclusively(cleanup):
    # Set up
    NUMBER_OF_WRITES = 5000

    #
    from multiprocessing import Pool
    with Pool(processes=36) as pool:
        pool.map(_worker, range(NUMBER_OF_WRITES))

    with open('deleteme.txt', 'r') as ofh:
        data = ofh.readlines()
    assert len(data) == NUMBER_OF_WRITES


# %%
if __name__ == '__main__':
    import pytest
    pytest.cmdlined.main(['test_helper.py', '-vsx', '--quick'])
