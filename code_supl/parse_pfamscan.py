# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from os.path import isfile
from collections import defaultdict, Counter

import cPickle as pickle
import pandas as pd


class parse_header_uniprot():
    pass

    def get_uniprotID(self, header_seq):
        return header_seq.split('|')[1].split('-')[0]
    
    def get_dict_key(self, header_seq):
        uniprotID = self.get_uniprotID(header_seq)
        splicingID = header_seq.split('|')[1]
        return (uniprotID, splicingID)


class parse_header_splicing():
    pass

    def get_uniprotID(self, header_seq):
        return header_seq.split('//')[1].split('-')[0]
    
    def get_dict_key(self, header_seq):
        uniprotID = self.get_uniprotID(header_seq)
        splicingID = header_seq.split('//')[0]
        return (uniprotID, splicingID)


class make_uniprot_pfam_database():
    
    def __init__(self):
        self.sequence_dict = defaultdict(list)
        self.sequence_df = pd.DataFrame()
        
        
    def read_uniprot(self, filename, id_type):
        """ Convert sequences from fasta files into a dictionary of SeqRecord objects
        """
        # The header file is formatted differently depending on whether the file
        # comes from uniprot or from our colaborators
        if id_type == 'uniprot':
            parse_header = parse_header_uniprot()
        elif id_type == 'splicing':
            parse_header = parse_header_splicing()
        else:
            print('Unrecognised id_type!')
            return
        
        uniprot_fasta_filehandle = open(filename, 'r')
        
        for record in SeqIO.parse(uniprot_fasta_filehandle, "fasta"):
            # Extract uniprotID from the header
            dict_key = parse_header.get_dict_key(record.id)
            # Add record to the dictionary
            self.sequence_dict[dict_key].append(record)
            
        uniprot_fasta_filehandle.close()
    
    
    
    def read_pfam(self, filename, id_type):
        """ Convert sequences from pfam_scan.pl output into a dictionary of SeqRecord objects
        """
        # The header file is formatted differently depending on whether the file
        # comes from uniprot or from our colaborators
        if id_type == 'uniprot':
            parse_header = parse_header_uniprot()
        elif id_type == 'splicing':
            parse_header = parse_header_splicing()
        else:
            print('Unrecognised id_type!')
            return   
    
        # Information given by pfam_scan.pl in the header line
        annotation_keys = ['alignment_start', 'alignment_end', 'envelope_start', 'envelope_end', 
                           'hmm_acc', 'hmm_name', 'type', 'hmm_start', 'hmm_end', 'hmm_length', 
                           'bit_score', 'E_value', 'significance', 'clan']
        
        pfamscan_results_filehandle = open(filename, 'r')
        
        keep_sequence = False
        for line in pfamscan_results_filehandle:
            
            # Remove carriage return from the end of the line (while keeping spaces)
            line = line.rstrip('\r\n')
            
            # Skip comment lines that begin with '#' followed by a space, or empty lines
            if line.split() == [] or line.split()[0] == '#': 
                continue
            
            # The only lines that don't begin with '#' are the header lines
            if line[0] != '#':
                
                # Add the previously-processed sequence to the dictionary
                if keep_sequence:
                    self.sequence_dict[dict_key].append(SeqRecord(domain_sequence))
                    self.sequence_dict[dict_key][-1].id = seq_id
                    self.sequence_dict[dict_key][-1].name = seq_id
                    self.sequence_dict[dict_key][-1].annotations = annotations
                    self.sequence_dict[dict_key][-1].letter_annotations = letter_annotations
                
                # Initiate values for the new sequence
                values = line.split()
                seq_id = values[0]
                dict_key = parse_header.get_dict_key(values[0])
                annotations = {}
                for key, value in zip(annotation_keys, values[1:]):
                    annotations[key] = value
                letter_annotations = {}
                
                # Decide whether the sequence meets preliminary inclusion cutoffs
                uniprot_length = float(len(self.sequence_dict[dict_key][0].seq))
                domain_length = float(annotations['alignment_end']) - int(annotations['alignment_start'])
                overlap = domain_length / uniprot_length
                if overlap > 0.10 and int(annotations['significance']) == 1:
                    keep_sequence = True
                else:
                    keep_sequence = False
            
            # Add letter_annotations from subsequent lines
            if line[0:11] == '#HMM       ':
                letter_annotations['hmm'] = line[11:]
            if line[0:11] == '#MATCH     ':
                letter_annotations['match'] = line[11:]
            if line[0:11] == '#PP        ':
                letter_annotations['pp'] = line[11:]
            if line[0:11] == '#SEQ       ':
                domain_sequence = Seq(line[11:])
            if line[0:11] == '#CS        ':
                letter_annotations['cs'] = line[11:]
                
        pfamscan_results_filehandle.close()
    
    
    def remove_overlapping_domains(self):
        """ If two or more domains overlap by more than 5 AA, keep only the longest domain
        """
        for key in self.sequence_dict.keys():
            
            if len(self.sequence_dict[key]) == 1:
                # No Pfam domains were found
                continue
            
            # Keep track of indices of Pfam domain SeqRecords to be deleted
            delete_overlaps_index = set()
            
            for i in range(1, len(self.sequence_dict[key])): # -1 includes the last one
                domainA_counter = Counter(range(int(self.sequence_dict[key][i].annotations['alignment_start'])-1, 
                                                int(self.sequence_dict[key][i].annotations['alignment_end'])))
                
                for j in range(i+1, len(self.sequence_dict[key])): # -1 includes the last one
                    domainB_counter = Counter(range(int(self.sequence_dict[key][j].annotations['alignment_start'])-1, 
                                                    int(self.sequence_dict[key][j].annotations['alignment_end'])))
                    
                    overlap             = list((domainA_counter & domainB_counter).elements())
                    domainA_remainder   = list((domainA_counter - domainB_counter).elements())
                    domainB_remainder   = list((domainB_counter - domainA_counter).elements())
                    
                    if overlap >= 5:
                        # Pfam domains overlap, delete the second domain
                        delete_overlaps_index.add(j)
                        
                        if len(domainB_remainder) > len(domainA_remainder):
                            # DomainB is longer, replace domainA with domainB
                            self.sequence_dict[key][i] = self.sequence_dict[key][j]
                            
                        elif len(domainB_remainder) == len(domainA_remainder):
                            # Domains are of equal length, use E_value to break ties
                            if self.sequence_dict[key][j].annotations['E_value'] < self.sequence_dict[key][i].annotations['E_value']:
                                # DomainB has a lower (better) E_value, replace domainA with domainB
                                self.sequence_dict[key][i] = self.sequence_dict[key][j]
                
            # Delete the seqrecords that were marked for deletion (overlapping sequences that are shorter)
            # you have to delete from highest index to lowest, so that the indices of subsequent items don't change
            for del_index in sorted(delete_overlaps_index, reverse=True):
                del self.sequence_dict[key][del_index]

    
    def link_repeating_domains(self):
        pass
    
    
    def make_data_frame(self):
        uniprotID = []
        splicingID = []
        alignment_start = []
        alignment_end = []
        hmm_acc = []
        hmm_name = []
        hmm_type = [] # 'type'
        hmm_clan = [] # 'clan'
        for key in self.sequence_dict.keys():
            if len(self.sequence_dict[key]) == 1:
                # No Pfam domains were found
                uniprotID.append(key[0])
                splicingID.append(key[1])
                alignment_start.append('')
                alignment_end.append('')
                hmm_acc.append('')
                hmm_name.append('')
                hmm_type.append('')
                hmm_clan.append('')
            else:
                for seqrecord in self.sequence_dict[key][1:]:
                    # Save every Pfam domain that was found for this (uniprotID, spliceID) tuple
                    uniprotID.append(key[0])
                    splicingID.append(key[1])
                    alignment_start.append(seqrecord.annotations['alignment_start'])
                    alignment_end.append(seqrecord.annotations['alignment_end'])
                    hmm_acc.append(seqrecord.annotations['hmm_acc'])
                    hmm_name.append(seqrecord.annotations['hmm_name'])
                    hmm_type.append(seqrecord.annotations['type'])
                    hmm_clan.append(seqrecord.annotations['clan'])
                    
        # Compile results into a pandas DataFrame
        self.sequence_df = pd.DataFrame({'uniprotID': uniprotID,
                                    'splicingID': splicingID,
                                    'alignment_start': alignment_start,
                                    'alignment_end': alignment_end,
                                    'hmm_acc': hmm_acc,
                                    'hmm_name': hmm_name,
                                    'hmm_type': hmm_type,
                                    'hmm_clan': hmm_clan})
    
    
    def save_dict(self, filename):
        fh = open(filename, 'wb')
        pickle.dump(self.sequence_dict)
        fh.close()

        
    def save_df(self, filename):
        fh = open(filename, 'wb')
        pickle.dump(self.sequence_df)
        fh.close()

        
    def export_df_to_tsv(self, filename):
        fh = open(filename, 'w')
        self.sequence_df.to_csv(fh, sep='\t')
        fh.close()



db_path = '/home/alexey/working/databases/'
db_filename = 'yanqi_all_seqrecords1.pickle'
db_filename2 = 'yanqi_all_seqrecords2_delete_overlaps.pickle'
db_filename3 = 'yanqi_all_seqrecords3_link_repeats.pickle'
db_filename4 = 'yanqi_all_seqrecords3_df.pickle'

db = make_uniprot_pfam_database()

db.read_uniprot(db_path + 'yanqi-human-uniprot-with-varsplic.fasta', 'uniprot')
print len(db.sequence_dict)

db.read_uniprot(db_path + 'yanqi-splicing.fasta', 'splicing')
print len(db.sequence_dict)

db.read_pfam(db_path + 'yanqi-human-uniprot-with-varsplic.pfamscan', 'uniprot')
print len(db.sequence_dict)

db.read_pfam(db_path + 'yanqi-splicing.pfamscan', 'splicing')
print len(db.sequence_dict)

db.save_dict(db_path + db_filename)

db.remove_overlapping_domains()
db.save_dict(db_path + db_filename2)

db.link_repeating_domains()
db.save_dict(db_path + db_filename3)

db.make_data_frame()
db.save_df(db_path + db_filename4)
