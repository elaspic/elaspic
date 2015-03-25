#!/usr/bin/env python
"""
Created on Fri Mar  6 18:18:24 2015

@author: Alexey Strokach

In order to run this script from a Spyder console, you first need to import modeller using the
``import modeller`` command. If your paths are set up to use modeller from python 3.4, but
you're actually using python 2.7, you also need to update your ``sys.path`` variable (which is
equivalent to the ``$PYTHONPATH`` variable in bash).

.. code-block: python

import sys
sys.path = [
    c if "lib/x86_64-intel8/python3" not in c
    else '/'.join(c.split('/')[:-1]) + "/python2.5"
    for c in sys.path
]
import modeller


"""
from __future__ import unicode_literals

import os
import argparse

from elaspic.pipeline import Pipeline

# read which configFile to use
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', required=True)
parser.add_argument('-i', '--input_file')
parser.add_argument('-u', '--uniprot_id')
parser.add_argument('-m', '--mutations', nargs='?', default=['',])
parser.add_argument('-t', '--run_type', nargs='?', type=int, default=5)
parser.add_argument('-n', '--n_cores', nargs='?', type=int, default=1)
args = parser.parse_args()

assert os.path.isfile(args.config_file)
pipeline = Pipeline(args.config_file)

if args.input_file \
and os.path.isfile(args.input_file):
    uniprot_ids = []
    mutations = []
    with open(args.input_file, 'r') as fh:
        for line in fh:
            # Can skip lines by adding spaces or tabs before them
            if line[0][0] == ' ' or line[0][0] == '\t':
                continue

            row = [ l.strip() for l in line.split('\t') ]

            # AS: Mutation does not necessarily have to be specified
            if len(row) > 1:
                uniprot_id, mutation = row[0], row[1]
            elif len(row) == 1:
                uniprot_id = row[0]
                mutation = ''
            #
            uniprot_ids.append(uniprot_id)
            mutations.append(mutation)

elif args.uniprot_id:
    uniprot_ids = [args.uniprot_id,]
    mutations = args.mutations if args.mutations is not None else ''
else:
    error_message = (
        'Need to supply either a list of uniprot_mutation combos '
        'or a flatfile with the same!')
    raise Exception(error_message)

run_type = args.run_type
n_cores = args.n_cores

#~ # Run jobs
#~ for uniprot_id, mutation in zip(uniprot_ids, mutations):
    #~ print uniprot_id
    #~ print mutation
    #~ print run_type
    #~ uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutation, run_type, n_cores)

# Run jobs
# uniprot_id, mutations as a comma-separated string,
uniprot_domains_and_domain_pairs = pipeline(uniprot_ids[0], mutations, run_type, n_cores)

