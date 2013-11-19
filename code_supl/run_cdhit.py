# -*- coding: utf-8 -*-

import pandas as pd
import cPickle as pickle


db_path = '/home/alexey/working/databases/uniprot-yanqi/'
db_filename1 = 'yanqi_all_seqrecords1.pickle'
db_filename2 = 'yanqi_all_seqrecords2:delete_overlaps.pickle'

fh = open(db_path + db_filename2, 'rb')
sequence_db = pickle.load(fh)
fh.close()



db_path = '/home/alexey/working/databases/sebastian-daniel-alexey/'


input_filename = 'pfamA_sebastian_select_columns.tsv'
output_filename = 'pfamA_sebastian_list_of_repeat_domains.tsv'

fh = open(db_path + input_filename,'r')
my_df = pd.read_csv(fh, sep='\t', header=None, 
                    names=['autoPfamA', 'pfamID', 'pfam_name', 
                    'pfam_description', 'pfam_source', 'pfam_type']);
fh.close()

my_df = my_df[my_df['pfam_type'] == 'Repeat']


#my_df.to_cvs()
#my_db2.to_csv()
#db_filename2 = 'yanqi_all_seqrecords2_delete_overlaps.pickle'
#db_filename3 = 'yanqi_all_seqrecords3_link_repeats.pickle'
#db_filename4 = 'yanqi_all_seqrecords3_df.pickle'
#db = make_uniprot_pfam_database()