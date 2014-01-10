# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import defaultdict, OrderedDict, Counter

import cPickle as pickle
import pandas as pd
import numpy as np


# Set up the logger (not the best place to do this, will fix later)
import logging
logging.basicConfig(filename = '/home/kimlab1/strokach/working/databases/uniprot-yanqi/parse_pfamscan.log',
                    filemode = 'w',
                    delay = True,
                    level = logging.CRITICAL, # highest level... DEBUG / INFO / WARNING / ERROR / CRITICAL
                    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')


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
    
    def __init__(self, pfamA_clans_filename='/home/kimlab1/strokach/working/databases/pfam.janelia.org/Pfam-A.clans.tsv',
                 pfamA_repeating_domains_filename='/home/kimlab1/strokach/working/databases/mysql-query-outputs/pfamA_repeating_domains.txt',
                 pfamA_superdomains_filename='/home/kimlab1/strokach/working/databases/mysql-query-outputs/pfamA_superdomains.txt'):
        """ Initialise files required for merging repeating domains and superdomains
        """
        
        # Initiate data types for storing information
        self.sequence_dict = defaultdict(list)
        self.sequence_df = pd.DataFrame()
        
        # Read the "Pfam_A.clans.tsv" file from pfam.janelia.org
        # key: pfam_id, value: clan_id ('\N' for None)
        self.pfamA_clans = {}
        with open(pfamA_clans_filename, 'r') as fh:
            for line in fh:
                row = line.strip().split('\t')
                self.pfamA_clans[row[1]] = row[2]
        
        # Read in a list of repeating domains, separated by "\n"
        with open(pfamA_repeating_domains_filename) as fh:
            self.pfamA_repeating_domains = [line.strip() for line in fh.readlines()]
        # Delete the header line
        del self.pfamA_repeating_domains[0]

        # Read in a list of superdomains (domains linked by +), separated by "\n"
        with open(pfamA_superdomains_filename) as fh:
            pfamA_superdomains = [line.strip() for line in fh.readlines()]
        # Delete the header line
        del pfamA_superdomains[0]
        # Add all partial combinations of superdomains to the set of superdomains
        pfamA_superdomains = set(pfamA_superdomains)
        self.pfamA_superdomains = pfamA_superdomains.copy()
        for domain in pfamA_superdomains:
            subdomains = domain.split('+')
            for i in range(2, len(subdomains)):
                for j in range(0, len(subdomains) + 1 - i):
                    self.pfamA_superdomains.add('+'.join(subdomains[j:i+j]))



    def run(self):
        """
        Read in human uniprot and splicing fasta and pfam_scan.pl output files.
        Link repeating and overlapping domains.
        """
        db_path = '/home/kimlab1/strokach/working/databases/uniprot-yanqi/'
        
        self.read_uniprot(db_path + 'yanqi-human-uniprot-with-varsplic.fasta', 'uniprot')
        print len(self.sequence_dict)        
        self.read_uniprot(db_path + 'yanqi-splicing.fasta', 'splicing')
        print len(self.sequence_dict)
        self.read_pfamscan(db_path + 'yanqi-human-uniprot-with-varsplic.pfamscan', 'uniprot')
        print len(self.sequence_dict)
        self.read_pfamscan(db_path + 'yanqi-splicing.pfamscan', 'splicing')
        print len(self.sequence_dict)
        self.remove_overlapping_domains()
        print len(self.sequence_dict)
        self.link_repeating_domains()
        print len(self.sequence_dict)
        print "Finished reading and linking domains..."

        
        
    def _choose_header_parser(self, id_type):
        # The header file is formatted differently depending on whether the file
        # comes from uniprot or from our colaborators
        if id_type == 'uniprot':
            parse_header = parse_header_uniprot()
        elif id_type == 'splicing':
            parse_header = parse_header_splicing()
        else:
            print('Unrecognised id_type!')
            return
        return parse_header


        
    def read_uniprot(self, fasta_filename, id_type):
        """ Convert sequences from fasta files into a dictionary of SeqRecord objects
        """
        logging.info('read_uniprot(self, ' + fasta_filename + ', ' + id_type + ')')
        
        parse_header = self._choose_header_parser(id_type)
        
        # Add sequences to internal data dictionary
        # key: tuple(uniprot_id, uniprot_splicing_id or collaborator_id)
        with open(fasta_filename, 'r') as fh:
            for record in SeqIO.parse(fh, "fasta"):
                # Extract uniprotID from the header
                dict_key = parse_header.get_dict_key(record.id)
                # Add record to the dictionary
                self.sequence_dict[dict_key].append(record)
      

    
    def read_pfamscan(self, pfamscan_results_filename, id_type):
        """
        Convert sequences from pfam_scan.pl output into a dictionary of SeqRecord objects
        
        pfamscan should be ran with the following settings:
        global evalue: < 1e-4
        c-evalue: < 1e-4
        i-evalue: < 1e-3
        This function removes domains that fall below pfamscan quality criteria (i.e. '?' instead of '!' in hmmer)
        """
        logging.info('read_pfamscan(self, ' + pfamscan_results_filename + ', ' + id_type + ')')
        
        parse_header = self._choose_header_parser(id_type)
    
        # Information given by pfam_scan.pl in the header line
        keep_sequence = False
        pfamscan_results_filehandle = open(pfamscan_results_filename, 'r')
        for line_number, line in enumerate(pfamscan_results_filehandle):
            
            # Remove carriage return from the end of the line (while keeping spaces)
            line = line.rstrip('\r\n')
            
            # Skip comment lines that begin with '#' followed by a space, or empty lines
            if line.split() == [] or line.split()[0] == '#': 
                continue
            
            # The only lines that don't begin with '#' are the header lines
            if line[0] != '#':
                
                # If the previously-processed sequence met the criteria for inclusion...
                if keep_sequence:
                    # Remove gaps in the pfam_scan alignment and corresponding annotations
                    while str(domain_sequence).find('-') != -1:
                        dash_index = str(domain_sequence).find('-')
                        domain_sequence = domain_sequence[:dash_index] + domain_sequence[dash_index+1:]
                        for key in letter_annotations:
                            letter_annotations[key] = letter_annotations[key][:dash_index] + letter_annotations[key][dash_index+1:]
                    assert (len(domain_sequence) == annotations['alignment_defs'][0][1]+1 - annotations['alignment_defs'][0][0])
                    
                    # Add the previously-processed sequence to the dictionary
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
        """
        1) Order the domains based on evalues (lower to higher)
        2) Remove a domain if it's overlap with the previous domain is > 10%
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
                    # Compare with the next domain (with a higher evalue)
                    j += 1
                # Finished looking at overlaps for the ith domain, move on to i+1
                i += 1


    
    def link_repeating_domains(self):
        """
        If consecutive domains are known to be repeats: 
        join them and replace hmm_name with clan_name, hmm_acc with clan_acc
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
                if (self.sequence_dict[key][i].id + '+' + self.sequence_dict[key][i+1].id 
                    in self.pfamA_superdomains):
                    # If domain A and domain B are found in PDBFam linked with +...
                    # ...then join them
                    self.sequence_dict[key][i] = self.__join_domains(
                                                    self.sequence_dict[key][0],
                                                    self.sequence_dict[key][i],
                                                    self.sequence_dict[key][i+1],
                                                    'super')
                    del self.sequence_dict[key][i+1]
                    continue
                # If no domains were deleted, move on to the next pair
                i += 1          
                       


    def __join_domains(self, seq_record_canonical, seq_record1, seq_record2, join_type):
        
        # The new domain name will be different depending on whether you are joining repeating domains or superdomains
        if join_type == 'super':
            # Link pfam names with a '+'
            new_seq_id = seq_record1.id + '+' + seq_record2.id
        elif join_type == 'repeating':
            assert(seq_record1.annotations['clan_name'][0] == seq_record2.annotations['clan_name'][0])
            if seq_record1.id == seq_record2.id:
                # If both pfam names are the same, keep that pfam name
                new_seq_id = seq_record1.id
            else:
                # Else use the clan name
                assert seq_record1.annotations['clan_name'][0]
                new_seq_id = seq_record1.annotations['clan_name'][0]


        # Since all annotations are lists, to get annotations for the merged domains you just join the lists
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
        
        
        # Cut the canonical uniprot sequence to the boundaries of the new domain
        new_domain_sequence = seq_record_canonical.seq[new_annotations['alignment_defs'][0][0]-1:new_annotations['alignment_defs'][-1][1]]
        
        seq_record1_alignment_length = seq_record1.annotations['alignment_defs'][-1][1] + 1 - seq_record1.annotations['alignment_defs'][0][0]
        seq_record2_alignment_length = seq_record2.annotations['alignment_defs'][-1][1] + 1 - seq_record2.annotations['alignment_defs'][0][0]
        
        domainA_counter = Counter(range(seq_record1.annotations['alignment_defs'][0][0]-1, seq_record1.annotations['alignment_defs'][-1][1]))
        domainB_counter = Counter(range(seq_record2.annotations['alignment_defs'][0][0]-1, seq_record2.annotations['alignment_defs'][-1][1]))
        
        overlap = list((domainA_counter & domainB_counter).elements())
        
        logging.debug('new_domain_sequence: ' + str(new_domain_sequence))
        logging.debug('seq_record1.seq: ' + str(seq_record1.seq))
        logging.debug('new_domain_sequence[:seq_record1_alignment_length]: ' + str(new_domain_sequence[:seq_record1_alignment_length]))
        logging.debug('seq_record2.seq: ' + str(seq_record2.seq))
        logging.debug('new_domain_sequence[-seq_record2_alignment_length:]: ' + str(new_domain_sequence[-seq_record2_alignment_length:]))

        assert(str(seq_record1.seq).upper() == str(new_domain_sequence[:seq_record1_alignment_length]).upper())
        assert(str(seq_record2.seq).upper() == str(new_domain_sequence[-seq_record2_alignment_length:]).upper())
        

        # Join the letter annotations so that they are of the same length
        # (maybe use them later on for better alignments?)
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

        
        # Make a seqrecord object for the merged domain
        new_seq_record = SeqRecord(new_domain_sequence)
        new_seq_record.id = new_seq_id
        new_seq_record.name = new_seq_id
        new_seq_record.annotations = new_annotations
        new_seq_record.letter_annotations = new_letter_annotations
        
        return new_seq_record
    
   
    def get_dataframe(self, canonical=True):
        """ Extract the most important information from the dictionary and 
        convert it into a dataframe
        This can be reversed using: 
        {k: list(v) for k,v in df.groupby("Address")["ID"]}
        (http://stackoverflow.com/a/20112846)
        """
        temp_dict = OrderedDict()
        temp_dict['uniprot_id'] = []
#        temp_dict['splicing_id'] = []
        temp_dict['seq_id'] = []
        temp_dict['alignment_def'] = []
        temp_dict['alignment_defs'] = []
        temp_dict['hmm_acc'] = []
        temp_dict['hmm_name'] = []
        temp_dict['hmm_type'] = [] # 'type'
        temp_dict['hmm_clan_acc'] = [] # 'clan'
        temp_dict['hmm_clan_name'] = []
        
        for key in self.sequence_dict.keys():
            # Output data only for canonical sequences
            if canonical and key[0] != key[1]:
                continue
                
            if len(self.sequence_dict[key]) == 1:
                # If no pfam domains were found, still keep an empty record 
                # for future reference 
                temp_dict['uniprot_id'].append(key[0])
#                temp_dict['splicing_id'].append(key[1])
                temp_dict['seq_id'].append('')
                temp_dict['alignment_def'].append('')
                temp_dict['alignment_defs'].append('')
                temp_dict['hmm_acc'].append('')
                temp_dict['hmm_name'].append('')
                temp_dict['hmm_type'].append('')
                temp_dict['hmm_clan_acc'].append('')
                temp_dict['hmm_clan_name'].append('')
                continue
            
            for counter, seqrecord in enumerate(self.sequence_dict[key]):
                # The first record in the list contains the canonical uniprot sequence
                if counter == 0:
                    continue
                
                # Save every Pfam domain that was found for this (uniprotID, spliceID) tuple
                temp_dict['uniprot_id'].append(key[0])
#                temp_dict['splicing_id'].append(key[1])
                temp_dict['seq_id'].append(seqrecord.id)
                temp_dict['alignment_def'].append((seqrecord.annotations['alignment_defs'][0][0],seqrecord.annotations['alignment_defs'][-1][1],))
                temp_dict['alignment_defs'].append(seqrecord.annotations['alignment_defs'])
                temp_dict['hmm_acc'].append(seqrecord.annotations['hmm_acc'])
                temp_dict['hmm_name'].append(seqrecord.annotations['hmm_name'])
                temp_dict['hmm_type'].append(seqrecord.annotations['type'])
                temp_dict['hmm_clan_acc'].append(seqrecord.annotations['clan_acc'])
                temp_dict['hmm_clan_name'].append(seqrecord.annotations['clan_name'])
        
        # Compile results into a pandas DataFrame
        pfam_df = pd.DataFrame(temp_dict)
        return pfam_df
    
        
    def pickle_dict(self, filename):
        with open(filename, 'wb') as fh:
            pickle.dump(self.sequence_dict, fh)



if __name__ == "__main__":

    # Instantiate the database object
    pfamA_clans_filename = '/home/kimlab1/strokach/working/databases/pfam.janelia.org/Pfam-A.clans.tsv'
    pfamA_repeating_domains_filename = '/home/kimlab1/strokach/working/databases/mysql-query-outputs/pfamA_repeating_domains.txt'
    pfamA_superdomains_filename = '/home/kimlab1/strokach/working/databases/mysql-query-outputs/pfamA_superdomains.txt'
    
    db = make_uniprot_pfam_database(pfamA_clans_filename, pfamA_repeating_domains_filename, pfamA_superdomains_filename)
    db.run()
    
    
    # Save data incrementally for debugging
    
#    db_path = '/home/kimlab1/strokach/working/databases/uniprot-yanqi/'
#    db_filename = 'yanqi_all_seqrecords1.pickle'
#    db_filename2 = 'yanqi_all_seqrecords2.delete_overlaps.pickle'
#    db_filename3 = 'yanqi_all_seqrecords3.link_repeats.pickle'
#    db_filename4 = 'yanqi_all_seqrecords3.link_repeats_df.pickle'
#    df_filename = 'yanqi_all_seqrecords3.link_repeats_df.tsv'
#    df_oject_filename = 'yanqi_all_seqrecords3.link_repeats_object.pickle'
#    
#    db.read_uniprot(db_path + 'yanqi-human-uniprot-with-varsplic.fasta', 'uniprot')
#    print len(db.sequence_dict)
#    
#    db.read_uniprot(db_path + 'yanqi-splicing.fasta', 'splicing')
#    print len(db.sequence_dict)
#    
#    db.read_pfamscan(db_path + 'yanqi-human-uniprot-with-varsplic.pfamscan', 'uniprot')
#    print len(db.sequence_dict)
#    
#    db.read_pfamscan(db_path + 'yanqi-splicing.pfamscan', 'splicing')
#    print len(db.sequence_dict)
#    db.pickle_dict(db_path + db_filename)
#    
#    db.remove_overlapping_domains()
#    print len(db.sequence_dict)
#    db.pickle_dict(db_path + db_filename2)
#    
#    db.link_repeating_domains()
#    print len(db.sequence_dict)
#    db.pickle_dict(db_path + db_filename3)
#    
#    db.make_data_frame()
#    db.pickle_dict(db_path + db_filename4)
#    
#    db.export_df_to_tsv(db_path + df_filename)
#    
#    with open(db_path + df_oject_filename, 'wb') as fh:
#        pickle.dump(db, fh)
