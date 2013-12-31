# -*- coding: utf-8 -*-

import logging
import gzip
import numpy as np
import os
import argparse

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionException
mmcif_parser = MMCIFParser()


# File paths
file_path = '/home/kimlab1/database_data/pdb/ftp/data/structures/divided/mmCIF'
output_folder = '/home/kimlab1/strokach/working/databases/mysql-query-outputs/RealPDBPfam27Table'
temp_filename = output_folder + '/' + 'tmp_structure.cif'

if False:
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('input_file', help='Full path to input file.')
    arguments = argument_parser.parse_args()
    input_file = arguments.input_file
else:
    input_file = output_folder + '/' + 'RealPDBFam27Table.tsv'


# Set up the logger
FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(filename = output_folder + '/' + "errors.log", 
                    filemode = 'w',
                    delay = True,
                    level = logging.DEBUG,
                    format=FORMAT)


# Some constants
aa_translations = {
'ALA': 'A',
'ARG': 'R',
'ASN': 'N',
'ASP': 'D',
'CYS': 'C',
'GLU': 'E',
'GLN': 'Q',
'GLY': 'G',
'HIS': 'H',
'ILE': 'I',
'LEU': 'L',
'LYS': 'K',
'MET': 'M',
'PHE': 'F',
'PRO': 'P',
'SER': 'S',
'THR': 'T',
'TRP': 'W',
'TYR': 'Y',
'VAL': 'V',
}


# Inputs
inputs = np.recfromtxt(input_file, delimiter='\t')

# Outputs
ids_full = []
errors_occured_full = []
sequences_full = []

for line in inputs:
    
    pdb_id = line[0]
    pdb_chain = line[3]
    domain_defs = line[4].split(',')
    domain_defs = [(int(domain.split(':')[0]), int(domain.split(':')[1])) for domain in domain_defs]
    
    #temp = np.recfromcsv('/home/alexey/working/databases/mysql-query-outputs/RealPDBFam27Table.tsv', delimiter='\t')
    
    error_list = []
    
    # Read the gunzipped file
    f_in = gzip.open(file_path + '/' + pdb_id[1:3] + '/' + pdb_id + '.cif.gz', 'rb')
    f_out = open(temp_filename, 'w')
    f_out.writelines(f_in.readlines())
    f_out.close()
    f_in.close()
    try:
        structure = mmcif_parser.get_structure(pdb_id, temp_filename)
    except PDBConstructionException:
        errors_occured = True
        print 'xxx'
        pass
    os.remove(temp_filename)
    
    # Get the correct pdb chain
    model = structure[0]
    try:
        chain = model[pdb_chain]
    except KeyError:
        print 'yyy'
        pass
    
    
    # Get the domain sequence and whether or not it looks suspicious
    domainAA_list = []
    passed_start = False
    reached_end = False
    errors_occured = False
    id_for_logger = pdb_id + '_' + pdb_chain + '_' + '_'.join(['-'.join([str(domain) for domain in domain_def]) for domain_def in domain_defs])
    for domain_def in domain_defs:
    
        for residue in chain.get_residues():  
            if residue.id[1] == domain_def[0]+1:
                passed_start = True
            if residue.id[1] == domain_def[1]+1:
                reached_end = True
                if not passed_start:
                    logging.debug(id_for_logger + ' reached the last AA of a domain without going through the first!')
                    errors_occured = True
                    
            if passed_start and residue.id[0] == ' ':
#                print residue.id
                if residue.resname in aa_translations.keys():
                    domainAA_list.append(aa_translations[residue.resname])
                else:
                    domainAA_list.append('-')
                    logging.debug(id_for_logger + ' an unusual AA was found in the domain!')
                    errors_occured = True
            
            if reached_end:
                break
                    
        if not reached_end:
            logging.debug(id_for_logger + ' reached the end of chain before reaching the last AA of domain!')
            errors_occured = True
        
        # Move on to the next disjoint segment of the domain
        passed_start = False
        reached_end = False
        domainAA_list.append('/')
    
    domainAA_string = ''.join(domainAA_list[:-1])
    
    if len(domainAA_string) - domainAA_string.count('/') != sum([y+1-x for (x,y) in domain_defs]):
        logging.debug(id_for_logger + ' length of AA sequence does not match the length of domain!')
        errors_occured = True
        
    if errors_occured:
        error_list.append((pdb_id, pdb_chain, domain_defs))

    ids_full.append(id_for_logger)
    errors_occured_full.append(errors_occured)
    sequences_full.append(domainAA_string)



    

#
#console_handler = StreamHandler()
#console_handler.setFormatter(default_formatter)
#
#error_handler = FileHandler(temp_folder + '/' + "errors.log", "w")
#error_handler.setLevel(logging.ERROR)
#error_handler.setFormatter(default_formatter)
#
#root = logging.getLogger()
#root.addHandler(console_handler)
#root.addHandler(error_handler)
#root.setLevel(logging.DEBUG)


#residues = chain.child_dict
#for domain_def in domain_defs:
#    for res_number in range(domain_def[0]+1, domain_def[1]+1): # pdbs start numbering from 1, domain_defs are from 0
#        if residues.has_key((' ', res_number, ' ')) and residues[(' ', res_number, ' ')].resname in aa_translations.keys():
#            domainAA_list.append(aa_translations[residues[(' ', res_number, ' ')].resname])
#        else:
#            domainAA_list.append('-')
#    domainAA_list.append('/')            