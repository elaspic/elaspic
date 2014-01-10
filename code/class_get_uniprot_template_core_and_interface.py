# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:25:45 2013

@author: niklas
"""

from class_pdbTemplate import pdbTemplate
import class_callTcoffee as tc
import class_error as error

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.PDB.PDBParser import PDBParser

#from Bio.Align import MultipleSeqAlignment

import gzip
from math import sqrt
from math import fabs

import datetime
import subprocess


class GetTemplate():
    """
    Parent class holding functions for finding the correct template given a
    uniprot sequence
    """
    def __init__(self, tmpPath, unique, pdbPath, savePDB, saveAlignments,
                  pool, semaphore, uniprot_sequence_database, pdb_resolution_database, 
                  domain_definition_database, log, path_to_archive):
        """
        input:
        tmpPath             type 'str'
        unique              type 'str'
        pdbPath             type 'str'
        savePDB             type 'str'
        saveAlignments      type 'str'
        pool                type class '__main__.ActivePool'
        semaphore           type class 'multiprocessing.synchronize.Semaphore'
        get_uniprot_seq     type __main__.getUniprotSequence instance
        """
        self.tmpPath = tmpPath
        self.unique = unique + '/'
        self.pdbPath = pdbPath
        self.savePDB = savePDB
        self.saveAlignments = saveAlignments
        self.pool = pool
        self.semaphore = semaphore
        
        self.uniprot_sequence_database = uniprot_sequence_database
        
        self.pdb_resolution_database = pdb_resolution_database
        self.domain_definition_database = domain_definition_database
        
        # get the logger from the parent and add a handler
        self.log = log
        self.path_to_archive = path_to_archive
        self.bad_pdbs = ['3C4D', '3LB7', '3NSV', '2NP8', '2WN0']
        
    
    def make_SeqRec_object(self, sequence, domain, ID):
        """
        return a Biopython SeqRec object
        cuts the sequence (given as string) to the domain boundaries and sets the ID
        
        input:
        sequence    type class 'Bio.SeqRecord.SeqRecord' or str
        domain      type 'list' of int
        ID          type 'str'

        return      type class 'Bio.SeqRecord.SeqRecord'
        """
        
        if isinstance(sequence, SeqRecord):
            sequence = sequence[domain[0]-1:domain[1]]
        elif isinstance(sequence, Seq):
            sequence = SeqRecord(sequence[domain[0]-1:domain[1]])
        elif isinstance(sequence, str):
            sequence = SeqRecord(Seq(sequence[domain[0]-1:domain[1]]))
        sequence.id = ID

        return sequence
    
#    def close(self):
#        self.get_uniprot_seq.close()
    
    ##############
    # for map_to_uniprot
    #
    def map_to_uniprot(self, uniprotKB, uniprot_sequence, domain_uniprot, pdbCode, chain, domain_pdb, saveAlignments, refine):
        """
        Align the uniprot and the pdb sequence 
        
        
        input:
        uniprotKB           type 'str'
        uniprot_sequence    type class 'Bio.SeqRecord.SeqRecord'
        domain_uniprot      type 'list' of 'int'
        pdbCode             type 'str'
        chain               type 'str'
        domain_pdb          type 'list' of 'int'
        saveAlignments      type 'str'
        refine              type 'bool'
        
        
        return:
        alignment           class 'Bio.Align.MultipleSeqAlignment'
        score               type 'float'
        pdb_sequence_id     type 'str'
        domain_uniprot      type 'list' of 'int'

        """
        # 1st: cut the uniprot sequence to the domain boundaries
        # 2nd: check if the uniprot sequences "overlaps", if yes, extend domain
        #      and repeat alignment

        
        # get the sequence from the pdb
        pdb_sequence, domain_pdb, chainNumberingDomain = self.get_pdb_sequence(pdbCode, chain, domain_pdb)
        
        
        
        # get the uniprot sequence
        uniprot_sequence_domain = self.make_SeqRec_object(uniprot_sequence, domain_uniprot, uniprotKB)
        
        # get the alignment
        alignment, score, pdb_sequence_id = self.map_to_uniprot_helper(uniprot_sequence_domain, pdb_sequence, saveAlignments)


        
        # AS Start ############################################################
        # It was assumed that the uniprot domain boundaries are expanded and are correct.
        # However when using pfamscan output, this assumption is not correct.
        
#        domain_uniprot_expanded = self.expand_boundaries(alignment, domain_uniprot, )
        
        
        # AS End ##############################################################


        
        # if it is the last round optimise the alignment
        if refine:
            alignment, score, cut_uniprot, cut_pdb = self.align_shorten(alignment, score, pdb_sequence_id, saveAlignments)
            domain_uniprot = [ domain_uniprot[0] + cut_uniprot[0], domain_uniprot[1] - cut_uniprot[1] ]


        return alignment, score, pdb_sequence_id, domain_uniprot



    def map_to_uniprot_helper(self, uniprot_sequence, pdb_sequence, saveAlignments):
        """
        input:
        uniprot_sequence    type class 'Bio.SeqRecord.SeqRecord'
        pdb_sequence        type class 'Bio.SeqRecord.SeqRecord'
        saveAlignments      type 'str'
        
        return:
        alignment           type class 'Bio.Align.MultipleSeqAlignment'
        score               type 'float'
        pdb_sequence.id     type 'str'
        """
        # write both sequences to one file
        seqFiles = open(self.tmpPath + self.unique + 'seqfiles.fasta', 'w')
        SeqIO.write([uniprot_sequence, pdb_sequence], seqFiles, 'fasta')
        seqFiles.close()
        
        seqIDs = [uniprot_sequence.id, pdb_sequence.id]
        
        # do the alignment and get the score
        alignment, score = self.do_align(seqIDs, saveAlignments)

        return alignment, score, pdb_sequence.id
        
        
        
    def expand_boundaries(self, alignment, domain_uniprot, ):
        """
        AS:
        
        input:
        
        output:
        
        """
        start = domain_uniprot[0]
        end = domain_uniprot[1]
        
        
        # subtract number of dashes from start, add number of dashes to end
        
        
        return
        
    
    
    def do_align(self, seqIDs, saveAlignments):
        """
        Align the sequences in the seqFile.fasta file and return the alignment
        and the percentage identity
        
        input
        seqIDs              type 'list'     ;look like ['P01112', '1FOEB']
        saveAlignments      type 'str'
        
        
        alignments[0]       type class 'Bio.Align.MultipleSeqAlignment'
        score               type 'float'
        """
        # use the try statement to ensure that the sempaphore is released
        # in the end
        self.log.debug("Aligning: " + seqIDs[0] + ':' + seqIDs[1])
        
        if self.semaphore:
            try:
                with self.semaphore:
                    # register to keep track of the running t_coffee instances
                    # the maximal number is set in the configfile.
                    self.pool.makeActive(self.unique)
                    tcoffee = tc.tcoffee_alignment(self.tmpPath, self.unique,
                                                   saveAlignments, 
                                                   [self.tmpPath + self.unique + 'seqfiles.fasta', ], 
                                                   seqIDs)
                    alignments = tcoffee.align()
            except:
                raise
            finally:
                    self.pool.makeInactive(self.unique)
        else:
            tcoffee = tc.tcoffee_alignment(self.tmpPath, self.unique,
                               saveAlignments, 
                               [self.tmpPath + self.unique + 'seqfiles.fasta', ], 
                               seqIDs)
            alignments = tcoffee.align()
        
        
        self.log.debug("Done aligning: " + seqIDs[0] + ':' + seqIDs[1])
        
        # see http://www.nature.com/nmeth/journal/v10/n1/full/nmeth.2289.html#methods
        # getting the score like Aloy did, based on sequence identity and coverage
        score_id = self.get_identity(alignments[0])
        score_cov = self.get_coverage(alignments[0])
        a = 0.95
        score = a*score_id/100.0 * score_cov/100 + (1.0-a)*score_cov/100
        
        return alignments[0], score
        
        

    def get_sequence_from_alignment(self, aligned_sequence):
        """
        Removes the gaps (i.e. '-') from the sequence and creates a
        Biopython SeqRecord object.
        
        input:
        aligned_sequence    type class 'Bio.SeqRecord.SeqRecord'

        return:
                            type class 'Bio.SeqRecord.SeqRecord'
        """
        seq = ''
        for i in range(len(aligned_sequence)):
            if aligned_sequence[i] != '-':
                seq = seq + aligned_sequence[i]
        result      = SeqRecord(Seq(seq))
        result.id   = aligned_sequence.id
        result.name = aligned_sequence.name

        return result
    
    
    def check_for_loners(self, uniprot_alignment, pdb_alignment, left_or_right):
        """
        For a quality control of the alignment 
        
        input:
        uniprot_alignment   type class 'Bio.SeqRecord.SeqRecord'
        pdb_alignment       type class 'Bio.SeqRecord.SeqRecord'
        left_or_right       type 'str'
        
        return:
                            type class 'Bio.SeqRecord.SeqRecord', 'str'
        """
        # first check for which sequence to check
        if left_or_right == 'left':
            for i in range(len(uniprot_alignment)):
                if uniprot_alignment[i] == '-':
                    return uniprot_alignment, 'pdb'
                elif pdb_alignment[i] == '-':
                    return pdb_alignment, 'uniprot'
        elif left_or_right == 'right':
            for i in range(len(uniprot_alignment)-1,-1,-1):
                if uniprot_alignment[i] == '-':
                    return uniprot_alignment, 'pdb'
                elif pdb_alignment[i] == '-':
                    return pdb_alignment, 'uniprot'
        else:
            print 'You must specify left or right correctly!'
            return


            
        
    def align_shorten(self, alignment, score, pdb_sequence_id, saveAlignments):
        """
        T-Coffee has an issue where it is placing some amino acids strangely
        far away from the main central block. That looks something like this:
            
            AA-------TTTT-XXXXXXXX---       pdb sequence
            XXXXXXXXXXXXXXXXXXXX-XXXX       uniprot sequence
        
        In such a case it is advisable to shorten the uniprot sequence to get
        a better alignment.
        
        input:
        alignment               type class 'Bio.Align.MultipleSeqAlignment'
        score                   type 'float'
        pdb_sequence_id         type 'str'
        saveAlignments          type 'str'
        
        
        return:
        alignment               type class 'Bio.Align.MultipleSeqAlignment'
        score                   type 'float'
        cut_left_all_uniprot    type 'int'
        cut_right_all_uniprot   type 'int'
        cut_left_all_pdb        type 'int'
        cut_right_all_pdb       type 'int'
        """

        uniprot_alignment, pdb_alignment = self.pick_sequence(alignment, pdb_sequence_id)

        cut_left_all_uniprot = 0
        cut_right_all_uniprot = 0
        cut_left_all_pdb = 0
        cut_right_all_pdb = 0
        
        # shorten as long as the alignment appears to have to large gaps
        cut = True
        while cut:
            cut, cut_left, cut_right, uniprot_or_pdb_left, uniprot_or_pdb_right, uniprot_sequence, pdb_sequence = self.align_cut(alignment, pdb_sequence_id)

            if uniprot_or_pdb_left == 'uniprot':
                cut_left_all_uniprot += cut_left
            elif uniprot_or_pdb_left == 'pdb':
                cut_left_all_pdb += cut_left
            if uniprot_or_pdb_right == 'uniprot':
                cut_right_all_uniprot += cut_right
            elif uniprot_or_pdb_right == 'pdb':
                cut_right_all_pdb += cut_right

            if cut:
                alignment, score, pdb_sequence.id = self.map_to_uniprot_helper(uniprot_sequence, pdb_sequence, saveAlignments)

        return alignment, score, [cut_left_all_uniprot, cut_right_all_uniprot], [cut_left_all_pdb, cut_right_all_pdb]
        
    
    
    def align_cut(self, alignment, pdb_sequence_id):
        """
        Checks for loners and shortens the sequence
        
        input:
        alignment               type class 'Bio.Align.MultipleSeqAlignment'
        pdb_sequence_id         type 'str'
        
        
        return:
        cut                     type 'bool'
        cut_left                type 'int'
        cut_right               type 'int'
        uniprot_or_pdb_left     type 'str'
        uniprot_or_pdb_right    type 'str'
        uniprot_sequence        type class 'Bio.SeqRecord.SeqRecord'
        pdb_sequence            type class 'Bio.SeqRecord.SeqRecord'
        """
        uniprot_alignment, pdb_alignment = self.pick_sequence(alignment, pdb_sequence_id)
        uniprot_sequence                 = self.get_sequence_from_alignment(uniprot_alignment)
        uniprot_sequence_old             = uniprot_sequence
        pdb_sequence                     = self.get_sequence_from_alignment(pdb_alignment)
        
        ## check for loners
        # left
        seq_left, uniprot_or_pdb_left   = self.check_for_loners(uniprot_alignment, pdb_alignment, 'left')
        do_shorten_left, loners_left    = self.check_loners(seq_left, 'left')
        # right
        seq_right, uniprot_or_pdb_right = self.check_for_loners(uniprot_alignment, pdb_alignment, 'right')
        do_shorten_right, loners_right  = self.check_loners(seq_right, 'right')

        
        cut = False
        cut_left = 0
        cut_right = 0
        # get the amount to cut left
        if do_shorten_left == True:
            cut = True
            seq_begin, seq_end, gap_end = loners_left
            # nr. of amino acids gone astray
            aa = seq_end - seq_begin
            # shorten by length of the gap minus 2 (or three to be genourous)
            # gap_end is the position where the first gap ends, thus the number
            # of amino acids that occured until that gap have to be taken into account
            cut_left = gap_end - aa
            assert( cut_left >= 0 )
        # get the amount to cut right
        if do_shorten_right == True:
            cut = True
            seq_begin, seq_end, gap_end = loners_right
            # nr. of amino acids gone astray
            aa = seq_begin - seq_end
            cut_right = (len(uniprot_alignment) - gap_end) - aa
            assert( cut_right >= 0 )
        

        # only shorten the uniprot sequence
        # Sebastian determined the domain boundaries. It is more reliable for
        # proteins with structure and the domain boundaries of the pdbs are
        # thus not modified. By uncommenting the bit below this could be changed.
        if cut:
            if cut_right == 0:
                if uniprot_or_pdb_right == 'uniprot' and uniprot_or_pdb_left == 'uniprot':
                    uniprot_sequence = uniprot_sequence[cut_left:]
                #elif uniprot_or_pdb_right == 'uniprot' and uniprot_or_pdb_left == 'pdb':
                #    pdb_sequence = pdb_sequence[cut_left:]
                elif uniprot_or_pdb_right == 'pdb' and uniprot_or_pdb_left == 'uniprot':
                    uniprot_sequence = uniprot_sequence[cut_left:]
                elif uniprot_or_pdb_right == 'pdb' and uniprot_or_pdb_left == 'pdb':
                    pass
                #    pdb_sequence = pdb_sequence[cut_left:]
                else:
                    print 'You must specify uniprot or pdb!'
                    print 'uniprot_or_pdb_right', uniprot_or_pdb_right
                    print 'uniprot_or_pdb_left', uniprot_or_pdb_left
            else:
                if uniprot_or_pdb_right == 'uniprot' and uniprot_or_pdb_left == 'uniprot':
                    uniprot_sequence = uniprot_sequence[cut_left:-cut_right]
                elif uniprot_or_pdb_right == 'uniprot' and uniprot_or_pdb_left == 'pdb':
                    uniprot_sequence = uniprot_sequence[:-cut_right]
                #    pdb_sequence = pdb_sequence[cut_left:]
                elif uniprot_or_pdb_right == 'pdb' and uniprot_or_pdb_left == 'uniprot':
                    uniprot_sequence = uniprot_sequence[cut_left:]
                #    pdb_sequence = pdb_sequence[:-cut_right]
                elif uniprot_or_pdb_right == 'pdb' and uniprot_or_pdb_left == 'pdb':
                #    pdb_sequence = pdb_sequence[cut_left:-cut_right]
                    pass
                else:
                    print 'You must specify uniprot or pdb!'
                    print 'uniprot_or_pdb_right', uniprot_or_pdb_right
                    print 'uniprot_or_pdb_left', uniprot_or_pdb_left
        
        # since only the uniprot sequence is shortend, but above it works in
        # principle for both sequences simultaniously, a endless while loop can occur
        # to avoid this, check if the uniprot sequence is cut or not.
        if len(uniprot_sequence_old.seq) == len(uniprot_sequence.seq):
            cut = False

        return cut, cut_left, cut_right, uniprot_or_pdb_left, uniprot_or_pdb_right, uniprot_sequence, pdb_sequence


    
    def pick_sequence(self, alignment, pdb_sequence_id):
        """
        Pick the uniprot and pdb sequences from the alignment
        
        input:
        alignment           type class 'Bio.Align.MultipleSeqAlignment'
        pdb_sequence_id     type 'str'
        
        return:
        uniprot_alignment   type class 'Bio.SeqRecord.SeqRecord'
        pdb_alignemnt       type class 'Bio.SeqRecord.SeqRecord'
        """
        for align in alignment:
            if align.id != pdb_sequence_id:
                uniprot_alignment = align
            else:
                pdb_alignemnt = align
        
        return uniprot_alignment, pdb_alignemnt

    
    def get_pdb_sequence(self, pdbCode, chain, domain_pdb):
        """
        Return the pdb file sequence (not SEQRES)
        
        input:
        pdbCode                 type 'str'
        chain                   type 'str'
        domain_pdb              type 'list' of 'int'
        
        result                  type class 'Bio.SeqRecord.SeqRecord'
        domain_pdb              type 'tuple' of 'int'
        chainNumberingDomain    type 'list' of 'int'
        """
        domains = [domain_pdb, ]
        pdb = pdbTemplate(self.pdbPath, pdbCode, chain, domains, self.tmpPath + self.unique)
        
        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()

        chainNumberingDomain = pdb.getChainNumberingNOHETATMS(chain)
        if chainNumberingDomain == []:
            raise error.pdbError('Could not get the pdb numbering for ' + pdbCode + '_' + chain)

        # it happend that the given domain boundaries are larger than the
        # chain found in the pdb file. In this case, set the boundaries to the
        # maximum of the pdb chain
        if domain_pdb[0] < chainNumberingDomain[0]:
            domain_pdb[0] = chainNumberingDomain[0]
        if domain_pdb[1] > chainNumberingDomain[-1]:
            domain_pdb[1] = chainNumberingDomain[-1]

        # translate the pdb_domain from pdb numbering to sequence numbering
        try:
            domain_pdb = chainNumberingDomain.index(domain_pdb[0])+1, chainNumberingDomain.index(domain_pdb[1])+1
        except ValueError:
            raise error.pdbError('ValueError when mapping domain boundaries to sequence numbering: ' + pdbCode + '_' + chain)


        pdb_sequence = next(SeqIO.parse(self.tmpPath + self.unique + pdbCode + chain + '.seq.txt', 'fasta'))
        return pdb_sequence, domain_pdb, chainNumberingDomain

    
#    def get_pdb_sequence(self, pdbCode, chain, domain_pdb):
#        """
#        Return the pdb file sequence (not SEQRES)
#        
#        input:
#        pdbCode                 type 'str'
#        chain                   type 'str'
#        domain_pdb              type 'list' of 'int'
#        
#        result                  type class 'Bio.SeqRecord.SeqRecord'
#        domain_pdb              type 'tuple' of 'int'
#        chainNumberingDomain    type 'list' of 'int'
#        """
#        domains = [domain_pdb, ]
#        pdb = pdbTemplate(self.pdbPath, pdbCode, chain, domains, self.tmpPath + self.unique)
#        
#        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()
#
#        chainNumberingDomain = pdb.getChainNumberingNOHETATMS(chain)
#        if chainNumberingDomain == []:
#            raise error.pdbError('Could not get the pdb numbering for ' + pdbCode + '_' + chain)
#
#        # it happend that the given domain boundaries are larger than the
#        # chain found in the pdb file. In this case, set the boundaries to the
#        # maximum of the pdb chain
#        if domain_pdb[0] < chainNumberingDomain[0]:
#            domain_pdb[0] = chainNumberingDomain[0]
#        if domain_pdb[1] > chainNumberingDomain[-1]:
#            domain_pdb[1] = chainNumberingDomain[-1]
#
#        try:
#            # Raise the first domain boundary if that AA does not exist
#            if chainNumberingDomain.count(domain_pdb[0]) != 1:
#                for chainNumber in chainNumberingDomain:
#                    if chainNumber > domain_pdb[0]:
#                        domain_pdb[0] = chainNumber
#                        break
#            # Lower the second domain boundary if that AA does not exist
#            if chainNumberingDomain.count(domain_pdb[1]) != 1:
#                for chainNumber in reversed(chainNumberingDomain):
#                    if chainNumber < domain_pdb[1]:
#                        domain_pdb[1] = chainNumber
#                        break
#            # Translate the pdb_domain from pdb numbering to sequence numbering                
#            domain_pdb = chainNumberingDomain.index(domain_pdb[0])+1, chainNumberingDomain.index(domain_pdb[1])+1
#        except ValueError:
#            print "PDB code, chain:", pdbCode, chain
#            print "Domains:", domains
#            print "Dobain_pdb 1/2:", domain_pdb[0], domain_pdb[1]
#            print "ChainNumberingDomain:"
#            print chainNumberingDomain            
#            raise error.pdbError('ValueError when mapping domain boundaries to sequence numbering: ' + pdbCode + '_' + chain)
#
#        pdb_sequence = next(SeqIO.parse(self.tmpPath + self.unique + pdbCode + chain + '.seq.txt', 'fasta'))
#        return pdb_sequence, domain_pdb, chainNumberingDomain    
    
    
    def get_identity(self, alignment):
        """
        Return the sequence identity of two aligned sequences
        It takes the longer sequence as the reference
        takes a biopython alignment object as input. make sure that the alignment
        has only two sequences
        
        input
        alignment       type class 'Bio.Align.MultipleSeqAlignment'
        
        return:
                        type float
        """
        assert( len(alignment) == 2 )
        
        length_seq1 = len(alignment[0].seq.tostring().replace('-', ''))
        length_seq2 = len(alignment[1].seq.tostring().replace('-', ''))
        
        j = 0 # counts positions in first sequence
        i = 0 # counts identity hits
    
        if length_seq1 <= length_seq2:
            for amino_acid in alignment[0].seq:
                if amino_acid == '-':
                    pass
                else:
                    if amino_acid == alignment[1].seq[j]:
                        i += 1
                j += 1
        else:
            for amino_acid in alignment[1].seq:
                if amino_acid == '-':
                    pass
                else:
                    if amino_acid == alignment[0].seq[j]:
                        i += 1
                j += 1
        
        return float(100*i/min(length_seq1, length_seq2))
    
    
    def get_coverage(self, alignment):
        """
        Returns the coverage of the alginment in % 
        
        input
        alignment       type class 'Bio.Align.MultipleSeqAlignment'>
        """
        assert( len(alignment) == 2 )
        length_seq1 = len(alignment[0].seq.tostring().replace('-', ''))
        length_seq2 = len(alignment[1].seq.tostring().replace('-', ''))
        
        if length_seq1 >= length_seq2:
            return 100.0 * float(length_seq2) / float(length_seq1)
        else:
            return 100.0 * float(length_seq1) / float(length_seq2)
    
    

    #
    # ends here
    ###########

    ###########
    # for map_to_pdb_sequence
    #
    def map_to_pdb_sequence(self, alignment, alignmentID, position):
        """
        Given an alignment and a position of the uniprot sequence, this function
        maps the position of the uniprot sequence to the pdb sequence
        !! Note that the position numbering start with 1 !!
        
        input
        alignment       type class 'Bio.Align.MultipleSeqAlignment'
        alignmentID     type 'str'
        position        type 'int'
        
        
        pdb_position    type 'int'

        """
        if alignment[0].id == alignmentID:
            alignment_protein = alignment[0]
            alignment_uniprot = alignment[1]
        elif alignment[1].id == alignmentID:
            alignment_protein = alignment[1]
            alignment_uniprot = alignment[0]
        else:
            print 'Could not assign the alignment to pdb and uniprot correctly!'
            return 1

        # now get the position
        pdb_position = 0
        uniprot_position = 0
        for index in range(len(alignment_uniprot)):
            if uniprot_position >= int(position)-1 and not alignment_uniprot[index] == '-':
                check = index
                break
            elif alignment_uniprot[index] == '-':
                pdb_position += 1
                continue
            elif alignment_protein[index] == '-':
                uniprot_position += 1
                continue
            else:
                pdb_position += 1
                uniprot_position += 1
        
        # check if the uniprot position falls into a gap
        if alignment_protein[check] == '-':
            return 'in gap'
        else:
            return pdb_position + 1


    #
    # ends here
    ###########
    
    #
    # for checking the alignment
    #########
    

    def check_loners(self, aligned_sequence, left_or_right):
        """
        Check if there are amino acids moved strangely far away from the central block
        
        check_loners
        aligned_sequence    type class 'Bio.SeqRecord.SeqRecord'
        left_or_right       type 'str'

        """
        amino_acid = 'RHKDESTNQCGPAVILMFYW'
        FIRST = True
        AMINO_ACID_START = False
        GAP_STARTED = False
        NO_GAP = True
        seq_begin = 0
        seq_end = 0
        gap_end = 0
        if left_or_right == 'left':
            indexes = range(len(aligned_sequence))
        elif left_or_right == 'right':
            indexes = range(len(aligned_sequence)-1,-1,-1)
        else:
            print 'You must specify left or right!'
        for i in indexes:
            if aligned_sequence[i] == '-' and FIRST:
                continue
            elif aligned_sequence[i] in amino_acid and FIRST:
                FIRST = False
                AMINO_ACID_START = True
                seq_begin = i
            elif aligned_sequence[i] == '-' and AMINO_ACID_START:
                loners_distance = 1
                AMINO_ACID_START = False
                GAP_STARTED = True
                seq_end = i
            elif aligned_sequence[i] in amino_acid and GAP_STARTED:
                gap_end = i
                NO_GAP = False
                break
            else:
                pass
        
        if NO_GAP:
            return False, []
            
        loners          = seq_end - seq_begin
        loners_distance = gap_end - seq_end
        

        if fabs(loners_distance) > fabs(loners):
            return True, [seq_begin, seq_end, gap_end]
        elif fabs(loners) / fabs(loners_distance) <= 0.2:
            return True, [seq_begin, seq_end, gap_end]
        else:
            return False, []



    def chose_best_template(self, domain_template):
        """
        Selects the best template based on sequence identity and resolution
        of the pdb structure
        
        input:
        templates       type list
        
        return:
        compare __call__() method
        """

        # First sort by identity score:
        domain_template = sorted(domain_template, 
                                 key=lambda k: k['alignment_scores'][0],
                                 reverse=True)
        max_score = domain_template[0]['alignment_scores']
        
        # Next, sort by pdb resolution
        best_domain_interactions = list()
        for interaction in domain_template:
            if interaction['alignment_scores'][0] == max_score:
                best_domain_interactions.append(interaction)
        best_domain_interactions = sorted(best_domain_interactions, 
                                          key=lambda k: (k['pdb_type'], k['pdb_resolution'],),
                                          reverse=False)
        
        return best_domain_interactions[0]

    
    
class GetTemplateCore(GetTemplate):
    """ retrieve the info from Sebastians file and do the alignment

    """
    
    def __init__(self, tmpPath, unique, pdbPath, savePDB, saveAlignments, pool, 
             semaphore, uniprot_sequence_database, pdb_resolution_database, 
             domain_definition_database, log, path_to_archive):
        """
        input
        tmpPath                     type 'str'
        unique                      type 'str'
        pdbPath                     type 'str'
        savePDB                     type 'str'
        saveAlignments              type 'str'  
        pool                        type class '__main__.ActivePool'
        semaphore                   type class 'multiprocessing.synchronize.Semaphore'
        get_uniprot_seq             <__main__.getUniprotSequence instance>
        core_template_database      <__main__.core_Templates instance>
        """
        # call the __init__ from the parent class
        GetTemplate.__init__(self, tmpPath, unique, pdbPath, savePDB, saveAlignments,
                  pool, semaphore, uniprot_sequence_database, pdb_resolution_database, 
                  domain_definition_database, log, path_to_archive)



    def __call__(self, protein_domain):
        """
        input:
        uniprotKB       type 'str'
        mutation        type 'str'
        
        return:
        pdbCode                             type 'str'
        chain                               type 'str'
        domain_pdb                          type 'list' of 'int'
        score                               type 'float'
        alignment                           type class 'Bio.Align.MultipleSeqAlignment'
        str(mutation_pdb)                   type 'str'
        uniprot_sequence_domain             type class 'Bio.SeqRecord.SeqRecord'
        mutation_position_domain_uniprot    type 'int'
        """
        # get all possible templates
        protein_domain_templates = self.run(protein_domain)
        
        self.log.info("Done getting domain templates...")
#        self.log.debug("Templates:")
#        self.log.debug(protein_interaction_templates)
        
        # select the best template
        best_template = self.chose_best_template(protein_domain_templates)
#        self.log.debug("Best templates:")
#        self.log.debug(best_template)
        
        best_template = self.run(protein_domain, best_template)
        self.log.info('Done refining interface templates...')

        # Export the alignment to the output folder
        path_to_data = (
            'domain' + '/' + 
            best_template['uniprot_id'] + '/' + 
            best_template['pfam_name'] + '*' + 
            '_'.join(['-'.join([str(i) for i in x]) for x in best_template['domain_def']]))
        
        # Folder for storing files for export to output
        tmp_save_path = self.tmpPath + self.unique + '/' + path_to_data + '/'
        subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
        subprocess.check_call('cp ' + self.alignments_path + best_template['alignment_name'] +
                                ' ' + tmp_save_path + best_template['alignment_name'], shell=True)
        if self.path_to_archive:
            archive_save_path = self.path_to_archive + path_to_data
            subprocess.check_call('mkdir -p ' + archive_save_path, shell=True)
            subprocess.check_call('cp ' + self.alignments_path + best_template['alignment_name'] +
                                    ' ' + archive_save_path + best_template['alignment_name'], shell=True)
        
        return best_template        
        
        

    def run(self, protein_domain, refine):
        """
        """
        
        # get the templates from Sebastians file
        domain_definitions = self.domain_definition_database.get_dicts(protein_domain['pfam_name'])
        
        if not domain_definitions:
            raise error.NoStructuralTemplates(str(datetime.datetime.now().date()) + ': no templates found')
        
        protein_definitions = []
        for definition in domain_definitions:
            if refine and definition['idx'] != refine['idx']:
                continue
            
             # there are some obsolete pdbs, ignore them... or 2NP8 has only one chain
            if definition['pdb_id'] in self.bad_pdbs:
                continue    
            
            definition['protein_definition_idx'] = protein_domain['idx']
            definition['domain_definition_idx'] = definition['idx']
            
            uniprot_sequence = self.uniprot_sequence_database(definition['uniprot_id'])
            if not uniprot_sequence:
                raise error.NoSequenceFound
           
            alignment, score, alignment_id, domain_def = \
                self.map_to_uniprot(definition['uniprot_id'], uniprot_sequence, 
                                    definition['alignment_def'], definition['pdb_id'],
                                    definition['pdb_chain'], definition['pdb_domain_def'],
                                    self.saveAlignments, refine)
            
            uniprot_sequence_domain = self.make_SeqRec_object(uniprot_sequence, domain_def, definition['uniprot_id'])
    
            pdb_type, pdb_resolution = self.pdb_resolution_database(definition['pdb_id'])
                            
            new_data = {
                'pdb_type': pdb_type,
                'pdb_resolution': pdb_resolution,
                'uniprot_domain_sequence': uniprot_sequence_domain,
                'alignment' : alignment,
                'alignment_id': pdb_sequence_id,
                'alignment_name': alignment[0].id + '_' + alignment[1].id + '.aln',
                'alignment_scores' : (score, score, 0,)}
            
            definition.update(new_data)
            protein_definitions.append(definition)

#       return (pdbCode, chain, domain_pdb, score, alignment, str(mutation_pdb), uniprot_sequence_domain, mutation_position_domain_uniprot, pfamID1, domain_uniprot), [[uniprot_id_1, uniprot_sequence], ]
        return protein_definitions



class GetTemplateInterface(GetTemplate):
    """
    Check if a mutation falls into an interface and tries to find the best
    structural template
    
    """
    def __init__(self, tmpPath, unique, pdbPath, savePDB, saveAlignments, pool, 
                 semaphore, uniprot_sequence_database, pdb_resolution_database, 
                 domain_definition_database, domain_interaction_database, log,
                 path_to_archive):
        """
        input
        
        tmpPath             type 'str'
        unique              type 'str'
        pdbPath             type 'str'
        savePDB             type 'str'
        saveAlignments      type 'str'
        pool                type class '__main__.ActivePool'
        semaphore           type class 'multiprocessing.synchronize.Semaphore'
        get_uniprot_seq     <__main__.getUniprotSequence instance>
        get_interactions    <__main__.getInteractions instance>
        get_3did_entries    <__main__.get3DID instance>
        get_resolution      <__main__.pdb_resolution instance>
        """
        # call the __init__ from the parent class
        GetTemplate.__init__(self, tmpPath, unique, pdbPath, savePDB, saveAlignments,
                  pool, semaphore, uniprot_sequence_database, pdb_resolution_database, 
                  domain_definition_database, log, path_to_archive)
        
        # set the databases
        self.domain_interaction_database = domain_interaction_database
    
    
    def __call__(self, protein_pair):
        """
        input:
        uniprot_id_1   type 'str'
        mutation    type 'str'
        
        return
        type 'dict'
        
        entries of one element of the list (same for the run() method):
        type        description
        'str'       pdbCode: of the template found
        'str'       chain1: chains to be used as template 
        'str'       chain2: 
        'float'     score: based on sequence identity and coverage
        'float'     score1: score of chain 1
        'float'     score2: score of chain 2
        'str'       pfam family of interaction partner 1 (i.e. the query uniprot), 
        'str'       pfam family of the interaction partner 2 (from Sebastians file)
        'str'       uniprot_id_1_2: of the second interaction partner
        'Bio.SeqRecord.SeqRecord'           uniprot_id_1_sequence1_domain: uniprot sequence for the domain only
        'Bio.SeqRecord.SeqRecord'           uniprot_id_2_sequence2_domain:
        'Bio.Align.MultipleSeqAlignment'    alignment1: biopython alignment obejct with the alignment of uniprot seq1 with chain1
        'Bio.Align.MultipleSeqAlignment'    alignment2:
        'int'       mutation_position_domain: mutation position in the domain (shortened pdb sequence)
        'str'       pdb_domain1, looks like '74-390'
        'str'       pdb_domain2, looks like '74-390'
        'int'       exp. measurement method; 0 means X-ray, 2 NMR, 3 other
        'float'     resolution of the structure
        """
        
        # get all possible templates
        protein_interaction_templates = self.run(protein_pair)
        self.log.info("Done getting interface templates...")
#        self.log.debug("Templates:")
#        self.log.debug(protein_interaction_templates)
        
        # select the best template
        best_template = self.chose_best_template(protein_interaction_templates)
#        self.log.debug("Best templates:")
#        self.log.debug(best_template)
        
        best_template = self.run(protein_pair, best_template)
        self.log.info('Done refining interface templates...')

        # Export the alignment to the output folder
        path_to_data = (
            'domain-complex' + '/' + 
            best_template['uniprot_id_1'] + '_' + best_template['uniprot_id_2'] + '/' + 
            best_template['pfam_name_1'] + '*' + 
            '_'.join(['-'.join([str(i) for i in x]) for x in best_template['domain_def_1']]) + '*' + 
            best_template['pfam_name_2'] + '*' + 
            '_'.join(['-'.join([str(i) for i in x]) for x in best_template['domain_def_2']]))
        
        # Folder for storing files for export to output
        tmp_save_path = self.tmpPath + self.unique + '/' + path_to_data + '/'
        subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
        subprocess.check_call('cp ' + self.alignments_path + best_template['alignment_name_1'] +
                                ' ' + tmp_save_path + best_template['alignment_name_1'], shell=True)
        subprocess.check_call('cp ' + self.alignments_path + best_template['alignment_name_2'] +
                                ' ' + tmp_save_path + best_template['alignment_name_2'], shell=True)
        if self.path_to_archive:
            archive_save_path = self.path_to_archive + path_to_data
            subprocess.check_call('mkdir -p ' + archive_save_path, shell=True)
            subprocess.check_call('cp ' + self.alignments_path + best_template['alignment_name_1'] +
                                    ' ' + archive_save_path + best_template['alignment_name_1'], shell=True)
            subprocess.check_call('cp ' + self.alignments_path + best_template['alignment_name_2'] +
                                    ' ' + archive_save_path + best_template['alignment_name_2'], shell=True)
        
        # Remember where you stored the data
        best_template['path_to_data'] = path_to_data
        
        return best_template
        
    
    def run(self, protein_pair, refine=None):
        """
        """
            
        domain_interactions = self.domain_interaction_database.get_dicts([protein_pair['pfam_name_1'], protein_pair['pfam_name_2']])
        
        if not domain_interactions:
            raise error.NoStructuralTemplates(str(datetime.datetime.now().date()) + ': no templates found')

        protein_interactions = []
        for interaction in domain_interactions:
            
            if refine and interaction['idx'] != refine['idx']:
                continue
            
            # there are some obsolete pdbs, ignore them... or 2NP8 has only one chain
            if interaction['pdb_id'] in self.bad_pdbs:
                continue
            
            interaction['protein_interaction_idx'] = protein_pair['idx']
            interaction['domain_interaction_idx'] = interaction['idx']
            
            if interaction['reversed'] is True:
                # The pair of uniprots appears in reverse order in the domain interaction database
                # In order to merge these two tables we have to flip protein-pair data
                interaction['uniprot_id_1'], interaction['uniprot_id_2'] = protein_pair['uniprot_id_2'], protein_pair['uniprot_id_1']
                interaction['pfam_name_1'], interaction['pram_name_2'] = protein_pair['pram_name_2'], protein_pair['pfam_name_1']
                interaction['aligment_def_1'], interaction['aligment_def_2'] = protein_pair['aligment_def_2'], protein_pair['aligment_def_1']
            else:
                # The pair of uniprots appears in reverse order in the domain interaction database
                # In order to merge these two tables we have to flip protein-pair data
                interaction['uniprot_id_1'], interaction['uniprot_id_2'] = protein_pair['uniprot_id_2'], protein_pair['uniprot_id_1']
                interaction['pfam_name_1'], interaction['pram_name_2'] = protein_pair['pram_name_2'], protein_pair['pfam_name_1']
                interaction['aligment_def_1'], interaction['aligment_def_2'] = protein_pair['aligment_def_2'], protein_pair['aligment_def_1']

            
            # Get sequences of the first and second protein          
            uniprot_sequence_1 = self.uniprot_sequence_database(interaction['uniprot_id_1'])
            if uniprot_sequence_1 == 'no sequence':
                raise error.NoSequenceFound('No sequence found for ' + interaction['uniprot_id_1'])
            
            uniprot_sequence_2 = self.uniprot_sequence_database(interaction['uniprot_id_2'])
            if uniprot_sequence_2 == 'no sequence':
                raise error.NoSequenceFound('No sequence found for ' + interaction['uniprot_id_2'])


            # Align the first and second uniprots to the corresponding pdb domains
            alignment_1, score_1, alignment_id_1, domain_def_1 = \
                self.map_to_uniprot(interaction['uniprot_id_1'], uniprot_sequence_1, 
                                    interaction['alignment_def_1'], interaction['pdb_id'], 
                                    interaction['pdb_chain_1'], interaction['pdb_domain_def_1'],
                                    self.saveAlignments, refine)

            alignment_2, score_2, alignment_id_2, domain_def_2 = \
                self.map_to_uniprot(interaction['uniprot_id_2'], uniprot_sequence_2, 
                                    interaction['alignment_def_2'], interaction['pdb_id'], 
                                    interaction['pdb_chain_2'], interaction['pdb_domain_def_1'], 
                                    self.saveAlignments, refine)
            
            
            # Cut sequence to boundaries and set sequence ID
            uniprot_domain_sequence_1 = self.make_SeqRec_object(uniprot_sequence_1, domain_def_1, interaction['uniprot_id_1'])
            uniprot_domain_sequence_2 = self.make_SeqRec_object(uniprot_sequence_2, domain_def_2, interaction['uniprot_id_1'])
            
            pdb_type, pdb_resolution = self.pdb_resolution_database(interaction['pdb_id'])
            
            # Save calculation results
            new_data = {'pdb_type': pdb_type,
                        'pdb_resolution': pdb_resolution,
                        'alignment_1': alignment_1,
                        'alignment_id_1': alignment_id_1,
                        'alignment_name_1': alignment_1[0].id + '_' + alignment_1[1].id + '.aln',
                        'domain_def_1': domain_def_1,
                        'uniprot_domain_sequence_1': uniprot_domain_sequence_1,
                        'alignment_2': alignment_2,
                        'alignment_id_2': alignment_id_2,
                        'alignment_name_2': alignment_2[0].id + '_' + alignment_2[1].id + '.aln',
                        'domain_def_2': domain_def_2,
                        'uniprot_domain_sequence_2': uniprot_domain_sequence_2,
                        'alignment_scores': (float(score_1) + float(score_2), float(score_1), float(score_2),)
                        }
            interaction.update(new_data)
            
            protein_interactions.append(interaction)
            
        # protein_interactions contains information about interactions between
        # two specific domains in two specific uniprots
        return protein_interactions


    ###########
    # for check_structure
    #
    def check_structure(self, pdbCode, chainID, mutation):
        """ checks if the mutation falls into the interface, i.e. is in contact with
        another chain
        'mutation' has to be of the form A_T70H, mutation in chain A, from Tyr at
        position 70 to His
        NOTE: takes the mutation as numbered ins sequence! The conversion is done
        within this function!
        
        input
        pdbCode     type 'str'
        chainID     type 'str'
        mutation    type 'str'      ; B_Q61L
        
        return:
        contacts    type 'dict'     ; {'C': False, 'B': True}
                                      key:   chainID                type 'str'
                                      value: contact to chainID     type boolean
        """
        structure = self.getPDB(pdbCode, self.pdbPath)
        model = structure[0]
    
        chains   = [ chain for chain in model]
        chainIDs = [ chain.id for chain in model]
    
        fromAA = mutation[2]
        position = self.convert_mutation_position(model, mutation) # convert the position numbering
        contacts = { chainID: False for chainID in chainIDs if not chainID == mutation[0] }
        
        for i in range(len(chains)):
            if chains[i].id == mutation[0]:
                chain = chains[i]
                # use list expansion to select only the 'opposing chains'
                oppositeChains = [ x for x in chains if x != chains[i] ]
        
        # If the residues do not match, issue a warning.
        # To obtain a better model one could restrict to templates that have
        # the same amino acid as the uniprot sequence at the position of the mutation.
#        if chain[position].resname != self.convert_aa(fromAA):
#            print 'Residue missmatch while checking the structure!'
#            print 'pdbCode', pdbCode
#            print 'mutation', mutation
#            print chain[position].resname, self.convert_aa(fromAA)
#            print 'position', position

       
        for oppositeChain in oppositeChains:
            # check each residue
            for residue in oppositeChain:
                # for each residue each atom of the mutated residue has to be checked
                for atom1 in chain[position]: # chain[position] is the residue that should be mutated
                    # and each atom
                    for atom2 in residue:
                        r = self.distance(atom1, atom2)
                        if r <= 5.0:
                            contacts[oppositeChain.id] = True

        return contacts


    def convert_mutation_position(self, model, mutation):
        """ maps the mutation sequence position of the pdb to pdb numbering
        
        input
        model       class 'Bio.PDB.Model.Model'
        mutation    type 'str'                      ; B_Q61L
        
        return:
        chainNumbering[position-1]      type 'int'
        """
        chain = model[mutation[0]]
        position = int(mutation[3:-1])
        
        chainNumbering = self.getChainNumberingNOHETATMS(chain)

        return chainNumbering[position-1]
    
    
    def getChainNumberingNOHETATMS(self, chain):
        """
        returns a list with the numbering of the chains
        
        input:
        chain               class 'Bio.PDB.Chain.Chain'
        
        return:
        chainNumbering      type 'list' of 'int'
        """
        chainNumbering = list()
        amino_acids = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
        for residue in chain:
            if residue.resname in amino_acids and residue.id[0] == ' ':
                chainNumbering.append(residue.id[1])

#        chainNumbering = [residue.id[1] for residue in chain if is_aa(residue, standard=True)]
        return chainNumbering
    
    
    def getPDB(self, pdbCode, pdbPath):
        """
        parse a pdb file with biopythons PDBParser() and return the structure
        
        input: pdbCode  type String     four letter code of the PDB file
        
        return: Biopython pdb structure
        
        input:
        pdbCode     type 'str'
        pdbPath     type 'str'
        
        return:
        result      type class 'Bio.PDB.Structure.Structure'
        """
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        pdbFile = pdbPath + pdbCode[1:3].lower() + '/pdb' + pdbCode.lower() + '.ent.gz'
        pdbFileUncompressed = gzip.open(pdbFile, 'r')
        result = parser.get_structure('ID', pdbFileUncompressed)

        return result
    
     
    def distance(self, atom1, atom2):
        """
        returns the distance of two points in three dimensional space
        
        input: atom instance of biopython: class 'Bio.PDB.Atom.Atom
        
        return:
                type 'float'
        """
        a = atom1.coord
        b = atom2.coord
        assert(len(a) == 3 and len(b) == 3)
        return sqrt(sum( (a - b)**2 for a, b in zip(a, b)))
    
    
    def convert_aa(self, aa):
        """
        input:
        aa      type 'str'
        
        return:
                type 'str'
        """
        A_DICT = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', \
                  'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', \
                  'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', \
                  'Y':'TYR', 'V':'VAL', 'U':'SEC', 'O':'PYL', \
                  'B':'ASX', 'Z':'GLX', 'J':'XLE', 'X':'XAA', '*':'TER'}
        
        AAA_DICT = dict([(value,key) for key,value in A_DICT.items()])
        
        if len(aa) == 3:
            try:
                return AAA_DICT[aa.upper()]
            except KeyError:
                print  'Not a valid amino acid'
                return
        if len(aa) == 1:
            try:
                return A_DICT[aa.upper()]
            except KeyError:
                print  'Not a valid amino acid'
                return
        print 'Not a valid amino acid'
    #
    # ends here
    ###########

if __name__ == '__main__':
    
    from pipeline import ActivePool
    from pipeline import getUniprotSequence
    import multiprocessing
    from Bio import AlignIO
    
    get_uniprot_seq = getUniprotSequence('doesntmatter')
    pool = ActivePool()
    semaphore = multiprocessing.Semaphore(2)
    
    tmpPath = '/tmp/test/'
    unique = 'unique'
    pdbPath = '/home/niklas/pdb_database/structures/divided/pdb/'
    savePDB = '/tmp/test/trash/'
    saveAlignments = '/tmp/test/alignments/'
    
    a = GetTemplate(tmpPath, unique, pdbPath, savePDB, saveAlignments, pool, semaphore, get_uniprot_seq)

    alignment = AlignIO.read('/home/niklas/tmp/alignment_cluttered.fasta', 'fasta')
#    print alignment
    score = '100'
    pdb_sequence_id = '1LFD'
    
    alignment, score = a.align_shorten(alignment, score, pdb_sequence_id, saveAlignments)
    
#    print alignment