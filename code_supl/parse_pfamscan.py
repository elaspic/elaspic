# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import defaultdict, Counter

import cPickle as pickle
import pandas as pd
import logging


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
    
    def __init__(self, pfamA_clans_filename, pfamA_repeating_domains_filename, pfamA_supra_domains_filename):
        
        # Initiate data types for storing information
        self.sequence_dict = defaultdict(list)
        self.sequence_df = pd.DataFrame()
        
        # Read the "Pfam_A.clans.tsv" file from pfam.janelia.org
        self.pfamA_clans = {}
        fh = open(pfamA_clans_filename, 'r')
        for line in fh:
            row = line.strip().split('\t')
            self.pfamA_clans[row[1]] = row[2]
        fh.close()
        
        # Read in a list of repeating domains, separated by "\n"
        fh = open(pfamA_repeating_domains_filename)
        self.pfamA_repeating_domains = [line.strip() for line in fh.readlines()]
        fh.close()
        # Delete the header line
        del self.pfamA_repeating_domains[0]

        # Read in a list of supra domains (domains linked by +), separated by "\n"
        fh = open(pfamA_supra_domains_filename)
        pfamA_supra_domains = [line.strip() for line in fh.readlines()]
        fh.close()
        # Delete the header line
        del pfamA_supra_domains[0]
        pfamA_supra_domains = set(pfamA_supra_domains)
        self.pfamA_supra_domains = pfamA_supra_domains.copy()
        for domain in pfamA_supra_domains:
            subdomains = domain.split('+')
            for i in range(2, len(subdomains)):
                for j in range(0, len(subdomains) + 1 - i):
                    self.pfamA_supra_domains.add('+'.join(subdomains[j:i+j]))
                
    
        
    def read_uniprot(self, filename, id_type):
        """ Convert sequences from fasta files into a dictionary of SeqRecord objects
        """
        logging.info('read_pfam(self, ' + filename + ', ' + id_type + ')')
        
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
        logging.info('read_pfam(self, ' + filename + ', ' + id_type + ')')
        
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
        pfamscan_results_filehandle = open(filename, 'r')
        
        keep_sequence = False
        for line_number, line in enumerate(pfamscan_results_filehandle):
            
            # Remove carriage return from the end of the line (while keeping spaces)
            line = line.rstrip('\r\n')
            
            # Skip comment lines that begin with '#' followed by a space, or empty lines
            if line.split() == [] or line.split()[0] == '#': 
                continue
            
            # The only lines that don't begin with '#' are the header lines
            if line[0] != '#':
                
                # Add the previously-processed sequence to the dictionary
                if keep_sequence:
                    while str(domain_sequence).find('-') != -1:
                        dash_index = str(domain_sequence).find('-')
                        domain_sequence = domain_sequence[:dash_index] + domain_sequence[dash_index+1:]
                        for key in letter_annotations:
                            letter_annotations[key] = letter_annotations[key][:dash_index] + letter_annotations[key][dash_index+1:]
                    assert (len(domain_sequence) == annotations['alignment_defs'][0][1]+1 - annotations['alignment_defs'][0][0])
                    self.sequence_dict[dict_key].append(SeqRecord(domain_sequence))
                    self.sequence_dict[dict_key][-1].id = seq_id
                    self.sequence_dict[dict_key][-1].name = seq_id
                    self.sequence_dict[dict_key][-1].annotations = annotations
                    self.sequence_dict[dict_key][-1].letter_annotations = letter_annotations
                
                # Initiate values for the new sequence
                values = line.split()
                seq_id = values[6]
                dict_key = parse_header.get_dict_key(values[0])
                   
                annotations = {}               
                annotations['alignment_defs'] = [(int(values[1]), int(values[2])), ]
                annotations['envelope_defs'] = [(int(values[3]), int(values[4])), ]
                annotations['hmm_acc'] = [values[5], ]
                annotations['hmm_name'] = [values[6], ]
                annotations['type'] = [values[7], ]
                annotations['hmm_defs'] = [(int(values[8]), int(values[9])), ]
                annotations['hmm_length'] = [int(values[10]), ]
                annotations['bit_score'] = [float(values[11]), ]
                annotations['E_value'] = [float(values[12]), ]
                annotations['significance'] = [int(values[13]), ]
                if values[14] != 'No_clan':
                    annotations['clan_acc'] = [values[14], ]
                    annotations['clan_name'] = [self.pfamA_clans[values[14]], ]
                else:
                    annotations['clan_acc'] = [None, ]
                    annotations['clan_name'] = [None, ]
                
                letter_annotations = {}
                
                if int(annotations['significance'][0]) == 1:
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
        logging.info('remove_overlapping_domains(self)')
        
        for key in self.sequence_dict.keys():
            
            if len(self.sequence_dict[key]) < 3:
                # One pfam domain or less, nothing to remove
                continue
            
            # Sort seqrecords in ascending order based on the evalue
            pfam_seqrecords_copy = self.sequence_dict[key][1:]
            pfam_seqrecords_copy = sorted(pfam_seqrecords_copy, key=lambda k: k.annotations['E_value'][0]) # Sort in ascending order
            self.sequence_dict[key][1:] = pfam_seqrecords_copy
            
            ## Delete overlapping domains
            i = 1 # index of domain A
            while i < len(self.sequence_dict[key])-1: # less than, so different than range
                domainA_counter = Counter(range(self.sequence_dict[key][i].annotations['alignment_defs'][0][0]-1, 
                                                self.sequence_dict[key][i].annotations['alignment_defs'][0][1]))
                j = i + 1 # index of domain B
                while j < len(self.sequence_dict[key]): # less than, so different than range
                    domainB_counter = Counter(range(self.sequence_dict[key][j].annotations['alignment_defs'][0][0]-1, 
                                                    self.sequence_dict[key][j].annotations['alignment_defs'][0][1]))
                    
                    overlap             = list((domainA_counter & domainB_counter).elements())
                    domainA_remainder   = list((domainA_counter - domainB_counter).elements())
                    domainB_remainder   = list((domainB_counter - domainA_counter).elements())
                    
                    if len(overlap) > 0.1 * len(domainA_counter):
                        # Domain A overlaps with domain B by more than 10%
                    
                        # Debugging message to track which domains are removed
                        if len(domainA_remainder) > 10 and len(domainB_remainder) > 10:
                            logging.debug('Domains with large overhangs were removed!')
                            logging.debug(key)
                            logging.debug('i = ' + str(i) + '   j = ' + str(j) + ';')
                            logging.debug('overlap = ' + str(len(overlap)))
                            logging.debug('domainA overhang = ' + str(len(domainA_remainder)))
                            logging.debug('domainB overhang = ' + str(len(domainB_remainder)) + '\n')
                                    
                        del self.sequence_dict[key][j]
                        continue
                    j += 1
                # Finished looking at overlaps for the ith domain, move on to i+1
                i += 1

    
    def link_repeating_domains(self):
        """ If consecutive domains are known to be repeats: join them and 
        replace hmm_name with clan_name, hmm_acc with clan_acc
        """
        logging.info('link_repeating_domains(self)')
        
        for key in self.sequence_dict.keys():
            
            if len(self.sequence_dict[key]) < 3:
                # One pfam domain or less, nothing to link
                continue
            
            # Sort seqrecords in ascending order based on the start of their domain in the uniprot sequence
            pfam_seqrecords_copy = self.sequence_dict[key][1:]
            pfam_seqrecords_copy = sorted(pfam_seqrecords_copy, key=lambda k: k.annotations['alignment_defs'][0][0]) # Sort in ascending order
            self.sequence_dict[key][1:] = pfam_seqrecords_copy
            
            ## Join repeating domains
            i = 1 # index of seq_record
            while i < len(self.sequence_dict[key])-1: # less than, so different than range
                if ((self.sequence_dict[key][i].id in self.pfamA_repeating_domains) and
                    (self.sequence_dict[key][i+1].id in self.pfamA_repeating_domains)):
                    # If both domain A and domain B are repeating domains...
                    if (self.sequence_dict[key][i].annotations['clan_acc'][0] ==
                        self.sequence_dict[key][i+1].annotations['clan_acc'][0]):   
                        # ...and domain A and domain B are in the same clan...                                
                        if self.sequence_dict[key][i].annotations['clan_acc'][0]:
                            # If that clan actually exists (is not None)...                                
                            # ...then join the domains
                            self.sequence_dict[key][i] = self.__join_domains(
                                                            self.sequence_dict[key][0],
                                                            self.sequence_dict[key][i],
                                                            self.sequence_dict[key][i+1],
                                                            'repeating')
                            del self.sequence_dict[key][i+1]
                            continue
                        else:
                            # The clan does not exist (is None)...
                            if (self.sequence_dict[key][i].id == 
                                self.sequence_dict[key][i+1].id):
                                # If domain A and domain B are the same domain...
                                # ...then join the domains
                                self.sequence_dict[key][i] = self.__join_domains(
                                                                self.sequence_dict[key][0],
                                                                self.sequence_dict[key][i],
                                                                self.sequence_dict[key][i+1],
                                                                'repeating')
                                del self.sequence_dict[key][i+1]
                                continue
                # If no domains were deleted, move on to the next pair
                i += 1
            
            ## Join supra domains
            i = 1 # index of seq_record
            while i < len(self.sequence_dict[key])-2:
                if (self.sequence_dict[key][i].id + '+' +
                    self.sequence_dict[key][i+1].id in self.pfamA_supra_domains):
                    # If domain A and domain B are found in PDBFam linked with +...
                    # ...then join them
                    self.sequence_dict[key][i] = self.__join_domains(
                                                    self.sequence_dict[key][0],
                                                    self.sequence_dict[key][i],
                                                    self.sequence_dict[key][i+1],
                                                    'supra')
                    del self.sequence_dict[key][i+1]
                    continue
                # If no domains were deleted, move on to the next pair
                i += 1          
                       

    def __join_domains(self, seq_record_canonical, seq_record1, seq_record2, join_type):
        
        if join_type == 'supra':
            new_seq_id = seq_record1.id + '+' + seq_record2.id
        elif join_type == 'repeating':
            assert(seq_record1.annotations['clan_name'][0] == seq_record2.annotations['clan_name'][0])
            if seq_record1.annotations['clan_name'][0]:
                new_seq_id = seq_record1.annotations['clan_name'][0]
            else:
                assert(seq_record1.id == seq_record2.id)
                new_seq_id = seq_record1.id

        
        new_annotations = {}
        new_annotations['alignment_defs'] = seq_record1.annotations['alignment_defs'] + seq_record2.annotations['alignment_defs']
        new_annotations['envelope_defs'] = seq_record1.annotations['envelope_defs'] + seq_record2.annotations['envelope_defs']
        new_annotations['hmm_acc'] = seq_record1.annotations['hmm_acc'] + seq_record2.annotations['hmm_acc']
        new_annotations['hmm_name'] = seq_record1.annotations['hmm_name'] + seq_record2.annotations['hmm_name']
        new_annotations['type'] = seq_record1.annotations['type'] + seq_record2.annotations['type']
        new_annotations['hmm_defs'] = seq_record1.annotations['hmm_defs'] + seq_record2.annotations['hmm_defs']
        new_annotations['hmm_length'] = seq_record1.annotations['hmm_length'] + seq_record2.annotations['hmm_length']
        new_annotations['bit_score'] = seq_record1.annotations['bit_score'] + seq_record2.annotations['bit_score']
        new_annotations['E_value'] = seq_record1.annotations['E_value'] + seq_record2.annotations['E_value']
        new_annotations['significance'] = seq_record1.annotations['significance'] + seq_record2.annotations['significance']
        new_annotations['clan_acc'] = seq_record1.annotations['clan_acc'] + seq_record2.annotations['clan_acc']
        new_annotations['clan_name'] = seq_record1.annotations['clan_name'] + seq_record2.annotations['clan_name']


        new_domain_sequence = seq_record_canonical.seq[new_annotations['alignment_defs'][0][0]-1:new_annotations['alignment_defs'][-1][1]]

        
        seq_record1_alignment_length = seq_record1.annotations['alignment_defs'][-1][1] + 1 - seq_record1.annotations['alignment_defs'][0][0]
        seq_record2_alignment_length = seq_record2.annotations['alignment_defs'][-1][1] + 1 - seq_record2.annotations['alignment_defs'][0][0]
        
        logging.debug('new_domain_sequence: ' + str(new_domain_sequence))
        logging.debug('seq_record1.seq: ' + str(seq_record1.seq))
        logging.debug('new_domain_sequence[:seq_record1_alignment_length]: ' + str(new_domain_sequence[:seq_record1_alignment_length]))
        logging.debug('seq_record2.seq: ' + str(seq_record2.seq))
        logging.debug('new_domain_sequence[-seq_record2_alignment_length:]: ' + str(new_domain_sequence[-seq_record2_alignment_length:]))

        assert(str(seq_record1.seq).upper() == str(new_domain_sequence[:seq_record1_alignment_length]).upper())
        assert(str(seq_record2.seq).upper() == str(new_domain_sequence[-seq_record2_alignment_length:]).upper())
        
        domainA_counter = Counter(range(seq_record1.annotations['alignment_defs'][0][0]-1, 
                                        seq_record1.annotations['alignment_defs'][-1][1]))
                                        
        domainB_counter = Counter(range(seq_record2.annotations['alignment_defs'][0][0]-1, 
                                        seq_record2.annotations['alignment_defs'][-1][1]))
                    
        overlap = list((domainA_counter & domainB_counter).elements())
        
        
        middle_pad = ''
        numAA_to_truncate = 0        
        if len(overlap) == 0:
            middle_pad = ':'*(seq_record2.annotations['alignment_defs'][0][0] - (seq_record1.annotations['alignment_defs'][-1][1] + 1))
        else:
            numAA_to_truncate = len(overlap)
        
        logging.debug('middle_pad: ' + str(middle_pad))
        logging.debug('numAA_to_truncate: ' + str(numAA_to_truncate))
        
        
        new_letter_annotations = {}
        new_letter_annotations['hmm'] = seq_record1.letter_annotations['hmm'] + middle_pad + seq_record2.letter_annotations['hmm'][numAA_to_truncate:]
        new_letter_annotations['match'] = seq_record1.letter_annotations['match'] + middle_pad + seq_record2.letter_annotations['match'][numAA_to_truncate:]
        new_letter_annotations['pp'] = seq_record1.letter_annotations['pp'] + middle_pad + seq_record2.letter_annotations['pp'][numAA_to_truncate:]
        if seq_record1.letter_annotations.has_key('cs') and seq_record2.letter_annotations.has_key('cs'):
            new_letter_annotations['cs'] = seq_record1.letter_annotations['cs'] + middle_pad + seq_record2.letter_annotations['cs'][numAA_to_truncate:]

                 
        new_seq_record = SeqRecord(new_domain_sequence)
        new_seq_record.id = new_seq_id
        new_seq_record.name = new_seq_id
        new_seq_record.annotations = new_annotations
        new_seq_record.letter_annotations = new_letter_annotations
        
        return new_seq_record
    
    
    def make_data_frame(self):
        uniprotID = []
        splicingID = []
        alignment_start = []
        alignment_end = []
        alignment_defs = []
        hmm_acc = []
        hmm_name = []
        hmm_type = [] # 'type'
        hmm_clan_acc = [] # 'clan'
        hmm_clan_name = []
        seq_id = []
        for key in self.sequence_dict.keys():
            if len(self.sequence_dict[key]) == 1:
                # No Pfam domains were found
                uniprotID.append(key[0])
                splicingID.append(key[1])
                alignment_start.append('')
                alignment_end.append('')
                alignment_defs.append('')
                hmm_acc.append('')
                hmm_name.append('')
                hmm_type.append('')
                hmm_clan_acc.append('')
                hmm_clan_name.append('')
                seq_id.append('')
            else:
                for seqrecord in self.sequence_dict[key][1:]:
                    # Save every Pfam domain that was found for this (uniprotID, spliceID) tuple
                    uniprotID.append(key[0])
                    splicingID.append(key[1])
                    alignment_start.append(seqrecord.annotations['alignment_defs'][0][0])
                    alignment_end.append(seqrecord.annotations['alignment_defs'][-1][1])
                    alignment_defs.append(seqrecord.annotations['alignment_defs'])
                    hmm_acc.append(seqrecord.annotations['hmm_acc'])
                    hmm_name.append(seqrecord.annotations['hmm_name'])
                    hmm_type.append(seqrecord.annotations['type'])
                    hmm_clan_acc.append(seqrecord.annotations['clan_acc'])
                    hmm_clan_name.append(seqrecord.annotations['clan_name'])
                    seq_id.append(seqrecord.id)
                    
        # Compile results into a pandas DataFrame
        self.sequence_df = pd.DataFrame({'uniprotID': uniprotID,
                                    'splicingID': splicingID,
                                    'alignment_start': alignment_start,
                                    'alignment_end': alignment_end,
                                    'alignment_defs': alignment_defs,
                                    'hmm_acc': hmm_acc,
                                    'hmm_name': hmm_name,
                                    'hmm_type': hmm_type,
                                    'hmm_clan_acc': hmm_clan_acc,
                                    'hmm_clan_name': hmm_clan_name,
                                    'seq_id': seq_id})
    
    
    def save_dict(self, filename):
        fh = open(filename, 'wb')
        pickle.dump(self.sequence_dict, fh)
        fh.close()

        
    def save_df(self, filename):
        fh = open(filename, 'wb')
        pickle.dump(self.sequence_df, fh)
        fh.close()

        
    def export_df_to_tsv(self, filename):
        fh = open(filename, 'w')
        self.sequence_df.to_csv(fh, sep='\t')
        fh.close()




# Set up logging file for debugging
db_path = '/home/alexey/working/databases/uniprot-yanqi/'

log_filename    = 'parse_pfamscan.log'

pfamA_clans_filename = '/home/alexey/working/databases/pfam.janelia.org/Pfam-A.clans.tsv'
pfamA_repeating_domains_filename = '/home/alexey/working/databases/mysql-query-outputs/pfamA_repeating_domains.txt'
pfamA_supra_domains_filename = '/home/alexey/working/databases/mysql-query-outputs/pfamA_supra_domains.txt'

db_filename     = 'yanqi_all_seqrecords1.pickle'
db_filename2    = 'yanqi_all_seqrecords2.delete_overlaps.pickle'
db_filename3    = 'yanqi_all_seqrecords3.link_repeats.pickle'
db_filename4    = 'yanqi_all_seqrecords3.link_repeats_df.pickle'


FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(filename = db_path + log_filename, 
                    filemode = 'w',
                    delay = True,
                    level = logging.INFO,
                    format=FORMAT)

db = make_uniprot_pfam_database(pfamA_clans_filename, pfamA_repeating_domains_filename, pfamA_supra_domains_filename)

#db.read_uniprot(db_path + 'yanqi-human-uniprot-with-varsplic.fasta', 'uniprot')
#print len(db.sequence_dict)
#
#db.read_uniprot(db_path + 'yanqi-splicing.fasta', 'splicing')
#print len(db.sequence_dict)
#
#db.read_pfam(db_path + 'yanqi-human-uniprot-with-varsplic.pfamscan', 'uniprot')
#print len(db.sequence_dict)
#
#db.read_pfam(db_path + 'yanqi-splicing.pfamscan', 'splicing')
#print len(db.sequence_dict)
#
#db.save_dict(db_path + db_filename)
#
#db.remove_overlapping_domains()
#db.save_dict(db_path + db_filename2)

fh = open(db_path + db_filename2, 'r')
db.sequence_dict = pickle.load(fh)
fh.close()

db.link_repeating_domains()
db.save_dict(db_path + db_filename3)

db.make_data_frame()
db.save_df(db_path + db_filename4)

#my_df[my_df.seq_id.str.contains('.+\+.+\+.+\+.+\+.+')]
