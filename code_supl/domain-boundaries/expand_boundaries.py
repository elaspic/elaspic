# -*- coding: utf-8 -*-

from collections import namedtuple
import cPickle as pickle


domain_domain_interaction_info = namedtuple('domain_domain_interaction_info',
"pdb_id pdb_chain_1 pdb_definition_1 pdb_contact_residues_1 pdb_chain_2 "
"pdb_definition_2 pdb_contact_residues_2 pdb_domain_sequence_1 pdb_domain_sequence_2 ")


uniprot_ids = ['P21802', 'Q02750', 'Q93074', 'P48735', 'P42771', 'P52333', 
'Q9UM73', 'P04198', 'Q16236', 'P46531', 'P04626', 'Q9Y4A5', 'Q02548', 'Q969H0', 
'P60484', 'P31749', 'P22607', 'P58012', 'P07949', 'P00533', 'Q9UKF5', 'Q5JWF2', 
'P46439', 'Q01081', 'Q13315', 'P16234', 'P42336', 'P63092', 'Q9UMS4', 'P08581', 
'P10721', 'O96017', 'Q8IXJ9', 'P40337', 'P30153', 'Q6N021', 'Q01628', 'Q15831', 
'P40763', 'P22681', 'P84243', 'P50148', 'Q99835', 'P35222', 'P04637', 'Q13485', 
'P36888', 'P40238', 'Q01130', 'P00519', 'Q9Y6K1', 'O60674', 'P29992', 'P63000', 
'Q06124', 'P01112', 'P01111', 'P01116', 'O75874', 'O95302', 'P15056', 'O75533', 
'Q92793']


class foo(object):
    
    def __init__(self):
        self.temp = domain_domain_interaction_info(*uniprot_ids[:9])
        
    def __call__(self):
        return self.temp
        
    def display(self):
        print self.temp
        
    def export_pickle(self):
        with open('deleteme.pickle', 'wb') as fh:
            pickle.dump(self.temp, fh)
            
            
bar = foo()


##

reload(pipeline)
temp = pipeline.domain_domain_interactions('/home/kimlab1/strokach/working/databases/mysql-query-outputs/domain_domain_interactions.tsv')

mismatch = []
temp_keys = temp.db.keys()[:100]
for key in temp_keys:
    temp2 = temp(key)
    domain_1_len = list(temp2)[0].pdb_definition_1[0][1] - list(temp2)[0].pdb_definition_1[0][0] + 1
    sequence_1_len = len(list(temp2)[0].pdb_domain_sequence_1)
    domain_2_len = list(temp2)[0].pdb_definition_2[0][1] - list(temp2)[0].pdb_definition_2[0][0] + 1
    sequence_2_len = len(list(temp2)[0].pdb_domain_sequence_2)
    print (sequence_1_len - domain_1_len), (sequence_2_len - domain_2_len), '\n'
    mismatch.append(sequence_1_len - domain_1_len)
    mismatch.append(sequence_2_len - domain_2_len)
    
