# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:08:57 2013

@author: niklas
"""

from os.path import isfile
from subprocess import check_call


def prepareModeller(outFile,
                    alignments,
                    targetIDs, 
                    templateIDs, 
                    HETATM_SEQs,
                    chains
                    ):
    
    # since the output file is open in append mode prior files are deleted
    if isfile(outFile):
        check_call('rm ' + outFile, shell=True)
    
    # add the HETATM dots intstead of a gap and an X
    # thats why position and position as boundary when opening the sequences
    # since the seqences can be 
    HETATMs = [ 0 for x in range(0,len(alignments)) ]

    # To check whether the HETATMs should be added as seperate chains, the number
    # of residues in the alignment has to be known.
    HETATMboundarys = [ 0 for x in range(0,len(alignments)) ]
    for index in range(len(HETATMboundarys)):
        for residue in alignments[index][1].seq:
            if residue != '-':
                HETATMboundarys[index] += 1
    
    # item in the following is a dictionary with the chainID as key
    # and the positions of the HETATMs in a list as values
    for key in HETATM_SEQs:
        i = chains.index(key.keys()[0])
        if key[chains[i]] == []:
            # i.e. no HETATMs are in the chain
            continue

        for position in key[chains[i]]:
            if position < HETATMboundarys[i]:
                # 'position' is given for the unaligned sequence.
                # The alignment might introduce gaps that have to be added to
                # the 'position' in order to add the HETATMs at the correct
                # position
                index = 0
                while True:
                    if index == position or index == 10000: # prevent a dead loop
                        break
                    try:
                        if alignments[i][1].seq[index] == '-':
                            position = position + 1
                    except IndexError:
                        print '------'
                        print 'alignments[i][1].seq', alignments[i][1].seq
                        print 'i', i
                        print 'length', len(alignments[i][1].seq)
                        print 'Error: index', index
                        raise
                    index += 1
                
                # if the target sequence is longer than the template (structure)
                # the HETATMS that should be attached at the end of the chain
                # have to be shifted
                # it could happen though that A.--A, and if one would just
                # shift the '.' one would get A--.A
                # only the following cases should be shifted:
                # AAAAA..AAAA   target
                # AAAAA..----   template
                # to obtain
                # AAAAAAAAA..   target
                # AAAAA----..   template
                position_backup = position
                position_add = 0
                while True:
                    try:
                        if alignments[i][1].seq[position + position_add] == '-':
                            position_add += 1
                        else:
                            break
                    except IndexError:
                        break
                position = position + position_add
                # if an index error occurs, the HETATMs have to added to the end
                try:
                    if alignments[i][1].seq[position] in 'RHKDESTNQCUGPAVILMFYW':
                        position = position_backup
                    alignments[i][0].seq = alignments[i][0].seq[:position] + '.' + alignments[i][0].seq[position:]
                    alignments[i][1].seq = alignments[i][1].seq[:position] + '.' + alignments[i][1].seq[position:]
                except IndexError:
                    alignments[i][0].seq = alignments[i][0].seq + '.'
                    alignments[i][1].seq = alignments[i][1].seq + '.'
                
                
            else:
                HETATMs[i] += 1

    # cut the alignment overhangs
    cutAlignments = list()
    for alignment in alignments:
        cutAlignments.append((cut_alignments(alignment, templateIDs)))
    
    alignments = cutAlignments
    
    
    if alignments[0][0].id == targetIDs[0] or alignments[0][0].id == targetIDs[1]:
        seq_records = [ alignments[x][0] for x in range(0,len(alignments)) ]
        write_fasta_to_modeller(outFile, seq_records, 'sequence', '_'.join(targetIDs), '.', '.', '.', '.', HETATMs)
        
        seq_records = [ alignments[x][1] for x in range(0,len(alignments)) ]
        templateID = templateIDs[0] + ''.join( [ item[-1] for item in templateIDs[1:] ] )
#        write_fasta_to_modeller(outFile, seq_records, 'structure', templateID, 'FIRST', '@', 'END', '@', HETATMs)
        write_fasta_to_modeller(outFile, seq_records, 'structure', templateID, '.', '.', '.', '.', HETATMs)
    if len(alignments) > 1:
        if alignments[0][1].id == targetIDs[0] or alignments[0][1].id == targetIDs[1]:
            seq_records = [ alignments[x][1] for x in range(0,len(alignments)) ]
            write_fasta_to_modeller(outFile, seq_records, 'sequence', '_'.join(targetIDs), '.', '.', '.', '.', HETATMs)
            
            seq_records = [ alignments[x][0] for x in range(0,len(alignments)) ]
            templateID = templateIDs[0] + ''.join( [ item[-1] for item in templateIDs[1:] ] )
#            write_fasta_to_modeller(outFile, seq_records, 'structure', templateID, 'FIRST', '@', 'END', '@', HETATMs)
            write_fasta_to_modeller(outFile, seq_records, 'structure', templateID, '.', '.', '.', '.', HETATMs)
                
    return



def write_fasta_to_modeller(outfile, 
                            inSeqObject, 
                            seqType, 
                            seqName, 
                            res1, 
                            chain1, 
                            res2, 
                            chain2, 
                            HETATM, 
                            ):
    """
    writes a sequence object to a file ready for use with modeller
    
    input:  outfile     type: string        the output file
            inSeqObject type: Biopython Sequence Record Object containing the sequence
            seqType     type: string        should be 'sequence' or 'structure'
            seqName     type: string        Name used by modeller for adressing the sequence
            res1        type: string        amino acid to start modelling (see modeller manual)
            chain1      type: string        chain to use for modelling (see modeller manual)
            res2        type: string
            chain2      type: string
            HETATM      type: int           Number of HETATM residues
           
    return: None
    """
    # set the chain break symbol ('/') if HETATM are added to the alignment
    # HETFlag is used to indicate if the HETATM should be written, and if the
    # chain break symbol '/' should be inserted
    HETFlag = [ 0 for x in range(0,len(HETATM)) ]
    sequence = ''
    FIRST = 0 # in case its the first sequence don't add a chain break (/)
    for i in range(0,len(inSeqObject)):
        if HETATM[i] > 0:
            HETFlag[i] = 1
        sequence = sequence + FIRST*'/' + str(inSeqObject[i].seq) + HETFlag[i]*'/' + HETATM[i]*'.'
        FIRST = 1
        
    # need to append since the function is called twice    
    with open(outfile, 'a') as f:
        f.write('>P1;' + seqName + '\n')
        f.write(seqType + ':' + seqName + ':' + res1 + ':' + chain1 + ':' + res2 + ':' + chain2 + '::::\n')
        f.write(sequence + '*')
        f.write('\n\n')
    return


        
def cut_alignments(alignment, templateIDs, OVERHANG=2):
    """ removes overhanging alignments
    """
    # cut only if the uniprot sequence is longer, i.e. overhanging
    if alignment[0].id in templateIDs:
        pdb_alignment = alignment[0]
    else:
        pdb_alignment = alignment[1]
        
    cut_left  = 0
    cut_right = 0
        # cut the beginning
    for i in range(len(pdb_alignment)):
        if pdb_alignment == '-':
            pass
        else:
            cut_left = i
            break
    
    # cut the end
    for i in reversed(range(len(pdb_alignment))):
        if pdb_alignment == '-':
            pass
        else:
            cut_right = i
            break
    
    # make sure that the cutting indices do not exceed the limits 0 and length
    # of the sequences (with '-'). String slicing with negative indices is wrong
    # but does not raise an error!!    
    if cut_left <= OVERHANG:
        cut_left = 0
    else:
        cut_left = cut_left-OVERHANG
    if len(alignment[0]) < cut_right + OVERHANG:
        cut_right = len(alignment[0])
    else:
        cut_right = cut_right+OVERHANG+1

    # make the new, shortened alignment object
    return alignment[:,cut_left:cut_right]




