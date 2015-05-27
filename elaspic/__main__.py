#!/usr/bin/env python
"""
Created on Fri Mar  6 18:18:24 2015
@author: Alexey Strokach


In order to run this script from a Spyder console, you first need to import modeller using the
``import modeller`` command. If your paths are set up to use modeller from python 3.4, but
you're actually using python 2.7, you also need to update your ``sys.path`` variable (which is
equivalent to the ``$PYTHONPATH`` variable in bash):

.. code-block: python

    import sys
    sys.path = [
        c if "lib/x86_64-intel8/python3" not in c
        else '/'.join(c.split('/')[:-1]) + "/python2.5"
        for c in sys.path
    ]
    import modeller


"""
#%%
from __future__ import unicode_literals

import os
import sys
import argparse

from elaspic.pipeline import Pipeline


#%% Parse arguments

def get_parser():

    description = """
    Run the ELASPIC pipeline
    """
    
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter)
        
    parser.add_argument(
        '-c', '--config_file', required=True, 
        help='ELASPIC configuration file')
    parser.add_argument(
        '-u', '--uniprot_id', 
        help="The uniprot_id of the protein that you want to mutate (e.g. 'P28223')")
    parser.add_argument(
        '-m', '--mutations', nargs='?', default=['',], 
        help="Mutation(s) that you wish to evaluate (e.g. 'D172E,R173H')")
    parser.add_argument(
        '-p', '--uniprot_domain_pair_ids', default='',
        help="List of uniprot_domain_pair_ids to analyse "
            "(useful if you want to restrict your analysis to only a handful of domains) " )
    parser.add_argument(
        '-f', '--input_file', 
        help=("A tab separated file of uniprot_ids, mutations, and optionally, uniprot_domain_pair_ids \n"
            "(optional; to be used instead of `--uniprot_id` and `--mutations`)"))
    parser.add_argument(
        '-t', '--run_type', nargs='?', type=int, default=5, choices=[1,2,3,4,5],
        help=('Type of analysis to perform: \n'
            '  1: Calculate Provean only \n'
            '  2: Create homololgy models only \n'
            '  3: Evaluate mutations only \n'
            '  4: Create homology models and evaluate mutations \n'
            '  5: Calculate Provean, create homology models, and evaluate mutations \n'))
    return parser



#%%
def validate_args(args):
    if not os.path.isfile(args.config_file):
        raise Exception('The configuration file {} does not exist!'.format(args.config_file))
        
    if args.input_file and not os.path.isfile(args.input_file):
        raise Exception('The input file {} does not exist!'.format(args.input_file))
        
    if ((args.uniprot_id is None and args.input_file is None) or 
        (args.uniprot_id is not None and args.input_file is not None)):
            raise Exception('One of `--uniprot_id` or `--input_file` must be specified!')


def parse_input_file(input_file):
    """
    Does not work! Do not use!
    """
    uniprot_ids = []
    mutations = []
    uniprot_domain_pair_ids = []
    with open(input_file, 'r') as fh:
        for line in fh:
            # Can skip lines by adding spaces or tabs before them
            if line[0][0] == ' ' or line[0][0] == '\t':
                print('Skipping line: {}'.format(line))
                continue
            row = [ l.strip() for l in line.split('\t') ]
            # Specifying the mutation is optional
            if len(row) == 2:
                uniprot_id, mutation, uniprot_domain_pair_id = row[0], row[1], row[2]
            elif len(row) == 1:
                uniprot_id, mutation, uniprot_domain_pair_id = row[0], row[1], ''
            elif len(row) == 1:
                uniprot_id, mutation, uniprot_domain_pair_id = row[0], '', ''
            uniprot_ids.append(uniprot_id)
            mutations.append(mutation)
            uniprot_domain_pair_ids.append(uniprot_domain_pair_id)
    return uniprot_ids, mutations, uniprot_domain_pair_ids


def main():
    parser = get_parser()
    args = parser.parse_args()
    validate_args(args)
    pipeline = Pipeline(args.config_file)
    if args.input_file:
        uniprot_ids, mutations, uniprot_domain_pair_ids = parse_input_file(args.input_file)
    else:
        uniprot_ids = [args.uniprot_id,]
        mutations = args.mutations if args.mutations else ''
        uniprot_domain_pair_ids = args.uniprot_domain_pair_ids
        
    uniprot_domain_pair_ids_asint = (
        [int(x) for x in uniprot_domain_pair_ids.split(',') if x]
        if uniprot_domain_pair_ids 
        else []
    )
    pipeline(
        uniprot_id=uniprot_ids[0],         
        mutations=mutations, 
        run_type=args.run_type, 
        uniprot_domain_pair_ids=uniprot_domain_pair_ids_asint,
    )

    #~ # Run jobs
    #~ for uniprot_id, mutation in zip(uniprot_ids, mutations):
        #~ print uniprot_id
        #~ print mutation
        #~ print run_type
        #~ uniprot_domains_and_domain_pairs = pipeline(uniprot_id, mutation, run_type, n_cores)
        

#%%
if __name__ == '__main__':
    main()

