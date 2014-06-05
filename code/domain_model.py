# -*- coding: utf-8 -*-

import os
import subprocess
import numpy as np
import pandas as pd

import pdb_template
from pdb_template import PDBTemplate
import errors
import call_modeller

import sql_db
import analyze_structure
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO


def get_pdb_structure(path_to_pdb_file):
    parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
    structure = parser.get_structure('ID', path_to_pdb_file)
    return structure


def get_unique_resnum_and_sequence(interactions):
    interacting_resnum = []
    interacting_aa = ''
    unique_resnums = set()
    for interaction in interactions:
        if interaction[0] not in unique_resnums:
            unique_resnums.add(interaction[0])
            interacting_resnum.append(interaction[0])
            interacting_aa += interaction[1]
    return interacting_resnum, interacting_aa

from os.path import isfile
from subprocess import check_call


def prepare_modeller(outFile, alignments, targetIDs, templateIDs, HETATM_SEQs, chains):
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



def write_fasta_to_modeller(outfile, inSeqObject, seqType, seqName, res1, chain1, res2, chain2, HETATM):
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


class GetModel(object):
    """
    """

    def __init__(
            self, global_temp_path, temp_path, unique, pdb_path,
            db, log, n_cores, subprocess_ids,
            modeller_runs):
        """
        Produces a model of a single uniprot domain or a domain pair

        Parameters
        ----------
        temp_path : 'str'
        unique : 'str'
        pdb_path : 'str'
        savePDB : 'str'

        Returns
        -------
        model : sql.object
            Object having all the information about the produced model. The model
            itself is stored in a file on disk.
        """
        self.global_temp_path = global_temp_path
        self.temp_path = temp_path
        self.unique = unique
        self.unique_temp_folder = temp_path + unique + '/'
        self.pdb_path = pdb_path
        self.db = db
        self.subprocess_ids =subprocess_ids
        # get the logger from the parent and add a handler
        self.log = log
        self.modeller_runs = modeller_runs
        self.n_cores = n_cores


    def __call__(self, uniprot_domain, uniprot_template):
        # Folder for storing files for export to output
        model_errors = []
        save_path = self.temp_path + str(uniprot_domain.path_to_data)
        self.log.debug('save path:')
        self.log.debug(save_path)

        if type(uniprot_domain) == sql_db.UniprotDomain:
            self.log.debug('template pdb:')
            self.log.debug(uniprot_template.domain.pdb_id)
#            # Get the canonical uniprot sequence
#            uniprot_sequence = self.db.get_uniprot_sequence(uniprot_domain.uniprot_id)
#            # Cut it to domain boundaries
#            uniprot_sequence_domain = uniprot_sequence[
#                sql_db.decode_domain(uniprot_template.domain_def)[0]-1:
#                sql_db.decode_domain(uniprot_template.domain_def)[1]]
#            sequences = [uniprot_sequence_domain, ]

            # Load previously-calculated alignments
            alignment, __ = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)

            # Put the data into lists so that it's independent of whether we
            # are looking at a core or interface mutation
            pdb_id = uniprot_template.domain.pdb_id
            chain_ids = [uniprot_template.domain.pdb_chain, ]
            pdb_domain_definitions = [sql_db.decode_domain(uniprot_template.domain.pdb_domain_def, return_string=True), ]
            alignments = [alignment, ]
            pdb = PDBTemplate(self.pdb_path, pdb_id, chain_ids, pdb_domain_definitions, save_path, self.unique_temp_folder, self.log)
            pdb.extract()

        if type(uniprot_domain) == sql_db.UniprotDomainPair:
            self.log.debug('template pdb:')
            self.log.debug(uniprot_template.domain_1.pdb_id)
#            # Get sequences of the first and second protein
            uniprot_sequence_1 = self.db.get_uniprot_sequence(uniprot_domain.uniprot_domain_1.uniprot_id)
            uniprot_sequence_2 = self.db.get_uniprot_sequence(uniprot_domain.uniprot_domain_2.uniprot_id)
#            # Cut sequence to boundaries and set sequence ID
#            uniprot_sequence_1_domain = uniprot_sequence_1[
#                sql_db.decode_domain(uniprot_template.domain_def_1)[0]-1:
#                sql_db.decode_domain(uniprot_template.domain_def_1)[1]]
#            uniprot_sequence_2_domain = uniprot_sequence_2[
#                sql_db.decode_domain(uniprot_template.domain_def_2)[0]-1:
#                sql_db.decode_domain(uniprot_template.domain_def_2)[1]]
#            sequences = [uniprot_sequence_1_domain, uniprot_sequence_2_domain]

            # Load previously-calculated alignments
            alignment_1, alignment_2 = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)

            pdb_id = uniprot_template.domain_1.pdb_id
            chain_ids = [uniprot_template.domain_1.pdb_chain, uniprot_template.domain_2.pdb_chain]
            pdb_domain_definitions = [sql_db.decode_domain(uniprot_template.domain_1.pdb_domain_def, return_string=True),
                                           sql_db.decode_domain(uniprot_template.domain_2.pdb_domain_def, return_string=True)]
            alignments = [alignment_1, alignment_2]
            pdb = PDBTemplate(self.pdb_path, pdb_id, chain_ids, pdb_domain_definitions, save_path, self.unique_temp_folder, self.log)
            pdb.extract()

        alignments_pdb_order = [None] * len(alignments)
        for alignment in alignments:
            alignment_chain = alignment[1].id[-1]
            alignment_idx = chain_ids.index(alignment_chain)
            alignments_pdb_order[alignment_idx] = alignment
        alignments = alignments_pdb_order

        target_ids = []
        target_seqs = []
        template_ids = []
        template_seqs = []
        alignment_chain_ids = []
        for alignment in alignments:
            target_ids.append(alignment[0].id)
            target_seqs.append(str(alignment[0].seq))
            template_ids.append(alignment[1].id)
            template_seqs.append(str(alignment[1].seq))
            alignment_chain_ids.append(alignment[1].id[-1])
        for alignment, alignment_chain_id, chain_id in zip(alignments, alignment_chain_ids, chain_ids):
            self.log.debug('Alignment: {}'.format(alignment))
            self.log.debug('Chain alignment order: {}; Chain pdb order: {}'.format(alignment_chain_id, chain_id))
            assert alignment and alignment_chain_id == chain_id

        target_id = '_'.join(target_ids)
        target_seq = '/'.join(target_seqs)
        # Keep the pdb id from the first template together with all chains
        template_id = template_ids[0] + ''.join( [ pdb_id_and_chain[-1] for pdb_id_and_chain in template_ids[1:] ] )
        template_seq = '/'.join(template_seqs)

        hetatm_chains = []
        for chain in pdb.structure[0].child_list:
            if chain.id not in chain_ids:
                hetatm_chains.append(chain)
        self.log.debug('Number of extra chains: {}'.format(len(hetatm_chains)))
        if len(hetatm_chains) > 1:
            raise Exception('Too many extra chains!')
        if len(hetatm_chains) == 1:
            target_seq += '/' + '.' * len(hetatm_chains[0])
            template_seq += '/' + '.'* len(hetatm_chains[0])

        # Write the alignment file for modeller
        pir_alignment_filename = self.unique_temp_folder + 'modeller/outFile_wildtype'
        with open(pir_alignment_filename, 'w') as pir_alignment_filehandle:
            self._write_to_pir_alignment(pir_alignment_filehandle, 'sequence', target_id, target_seq)
            self._write_to_pir_alignment(pir_alignment_filehandle, 'structure', template_id, template_seq)

        # Make the homology model and check if it is knotted
        norm_dope_wt, pdbFile_wt, knotted = self._get_model(pir_alignment_filename, target_ids, template_ids, save_path)
        self.log.debug('model pdb file: %s, knotted: %s' % (pdbFile_wt, knotted,))
        if knotted:
            model_errors.append('knotted')

        # Save another copy of the modelled structure in the tmp export folder
        subprocess.check_call('cp ' + self.unique_temp_folder + 'modeller/' + pdbFile_wt + ' ' + save_path + pdbFile_wt, shell=True)
        structure = get_pdb_structure(self.unique_temp_folder + 'modeller/' + pdbFile_wt)
        model = structure[0]

        # Save the results to the database
        if isinstance(uniprot_domain, sql_db.UniprotDomain):
            uniprot_model = sql_db.UniprotDomainModel()
            uniprot_model.uniprot_domain_id = uniprot_template.uniprot_domain_id
            uniprot_model.model_filename = pdbFile_wt
            uniprot_model.norm_dope = norm_dope_wt
            uniprot_model.chain = 'A'

            # Get solvent-accessible surface area (SASA) for each residue
            analyze_structure_object = analyze_structure.AnalyzeStructure(
                self.unique_temp_folder + 'modeller/',
                self.unique_temp_folder + 'analyze_structure/',
                pdbFile_wt, [uniprot_model.chain], None, self.log)
            # Run the homology model through msms and get dataframes with all the
            # per atom and per residue SASA values
            (seasa_by_chain_together, seasa_by_chain_separately,
            seasa_by_residue_together, seasa_by_residue_separately) = analyze_structure_object.get_seasa()
            # Get SASA only for amino acids in the main chain of interest
            msms_length_mismatch = False
            rel_sasa_for_each_residue = []
            for res in model[uniprot_model.chain]:
                if res.resname in pdb_template.amino_acids:
                    rel_sasa = seasa_by_residue_separately[
                        (seasa_by_residue_separately['pdb_chain']==uniprot_model.chain) &
                        (seasa_by_residue_separately['res_name']==res.resname) &
                        (seasa_by_residue_separately['res_num']==(str(res.id[1]) + res.id[2].strip()))].iloc[0]['rel_sasa']
                    if rel_sasa is None or rel_sasa is np.nan:
                        msms_length_mismatch = True
                    rel_sasa_for_each_residue.append(rel_sasa)
            uniprot_model.sasa_score = ','.join([str(x) for x in rel_sasa_for_each_residue])
            if msms_length_mismatch:
                self.log.error('msms score length mismatch')
                model_errors.append('msms score length mismatch')

        if isinstance(uniprot_domain, sql_db.UniprotDomainPair):
            uniprot_model = sql_db.UniprotDomainPairModel()
            uniprot_model.uniprot_domain_pair_id = uniprot_template.uniprot_domain_pair_id
            uniprot_model.model_filename = pdbFile_wt
            uniprot_model.norm_dope = norm_dope_wt
            uniprot_model.chain_1 = 'A'
            uniprot_model.chain_2 = 'B'

            modeller_chains = [uniprot_model.chain_1, uniprot_model.chain_2]
            interactions_between_chains = analyze_structure.get_interactions_between_chains(model, modeller_chains[0], modeller_chains[1], 5)
            chain_1_interactions = list(set([key[:2] for key in interactions_between_chains.keys()]))
            chain_2_interactions = list(set([value[:2] for values in interactions_between_chains.values() for value in values]))
            if not chain_1_interactions and not chain_2_interactions:
                raise errors.ChainsNotInteracting(
                    'Chains %s and %s are not interacting! chain_1_interactions: %s, chain_2_interactions: %s' %
                    (modeller_chains[0], modeller_chains[1], ','.join(list(chain_1_interactions)), ','.join(list(chain_2_interactions)),))

            chain_1_interactions.sort(key=lambda x: x[0]) # Sort by residue for easy reading
            chain_2_interactions.sort(key=lambda x: x[0]) # Sort by residue for easy reading

            chain_1_interacting_resnum, chain_1_interacting_aa = get_unique_resnum_and_sequence(chain_1_interactions)
            chain_1_numbering = pdb_template.get_chain_numbering(model[modeller_chains[0]], return_extended=True)
            chain_2_interacting_resnum, chain_2_interacting_aa = get_unique_resnum_and_sequence(chain_2_interactions)
            chain_2_numbering = pdb_template.get_chain_numbering(model[modeller_chains[1]], return_extended=True)

            self.log.debug('chain_1_interactions:')
            self.log.debug(chain_1_interactions)
            self.log.debug('chain_1_interacting_resnum:')
            self.log.debug(chain_1_interacting_resnum)
            self.log.debug('chain_1_interacting_aa:')
            self.log.debug(chain_1_interacting_aa)
#            self.log.debug('chain_1_numbering:')
#            self.log.debug(chain_1_numbering)

            self.log.debug('chain_2_interactions:')
            self.log.debug(chain_2_interactions)
            self.log.debug('chain_2_interacting_resnum:')
            self.log.debug(chain_2_interacting_resnum)
            self.log.debug('chain_2_interacting_aa:')
            self.log.debug(chain_2_interacting_aa)
#            self.log.debug('chain_2_numbering:')
#            self.log.debug(chain_2_numbering)

            chain_1_interacting_uninum = [
                sql_db.decode_domain(uniprot_template.domain_def_1)[0] + chain_1_numbering.index(resnum)
                for resnum in chain_1_interacting_resnum]
            chain_2_interacting_uninum = [
                sql_db.decode_domain(uniprot_template.domain_def_2)[0] + chain_2_numbering.index(resnum)
                for resnum in chain_2_interacting_resnum]

            chain_1_interactions_uniprot = zip(chain_1_interacting_uninum, chain_1_interacting_aa)
            chain_1_interactions_uniprot.sort(key=lambda x: x[0])
            chain_1_interacting_uninum, chain_1_interacting_aa = zip(*chain_1_interactions_uniprot)
            chain_1_interacting_aa = ''.join(chain_1_interacting_aa)
            chain_1_interacting_aa_from_uniprot = ''
            for uninum, aa in chain_1_interactions_uniprot:
                uniprot_idx = uninum - 1
                chain_1_interacting_aa_from_uniprot += uniprot_sequence_1[uniprot_idx]

            chain_2_interactions_uniprot = zip(chain_2_interacting_uninum, chain_2_interacting_aa)
            chain_2_interactions_uniprot.sort(key=lambda x: x[0])
            chain_2_interacting_uninum, chain_2_interacting_aa = zip(*chain_2_interactions_uniprot)
            chain_2_interacting_aa = ''.join(chain_2_interacting_aa)
            chain_2_interacting_aa_from_uniprot = ''
            for uninum, aa in chain_2_interactions_uniprot:
                uniprot_idx = uninum - 1
                chain_2_interacting_aa_from_uniprot += uniprot_sequence_2[uniprot_idx]

            self.log.debug('domain_def_1: {}'.format(uniprot_template.domain_def_1))
            self.log.debug('chain_1_interactions_uniprot: {}'.format(chain_1_interactions_uniprot))
            self.log.debug('uniprot_sequence_1: {}'.format(uniprot_sequence_1.seq))
            self.log.debug('chain_1_interacting_aa: {}'.format(chain_1_interacting_aa))
            self.log.debug('chain_1_interacting_aa_from_uniprot: {}'.format(chain_1_interacting_aa_from_uniprot))

            self.log.debug('domain_def_2: {}'.format(uniprot_template.domain_def_2))
            self.log.debug('chain_2_interactions_uniprot: {}'.format(chain_2_interactions_uniprot))
            self.log.debug('uniprot_sequence_2: {}'.format(uniprot_sequence_2.seq))
            self.log.debug('chain_2_interacting_aa: {}'.format(chain_2_interacting_aa))
            self.log.debug('chain_2_interacting_aa_from_uniprot: {}'.format(chain_2_interacting_aa_from_uniprot))

            if ((len(chain_1_interacting_uninum) == 0 or
                    len(chain_2_interacting_uninum) == 0)):
                raise errors.ChainsNotInteracting(
                    'Chains %s and %s are not interacting! chain_1_interactions: %s, chain_2_interactions: %s' %
                    (modeller_chains[0], modeller_chains[1], ','.join(list(chain_1_interactions)), ','.join(list(chain_2_interactions)),))

            if ((chain_1_interacting_aa_from_uniprot != chain_1_interacting_aa) or
                    (chain_2_interacting_aa_from_uniprot != chain_2_interacting_aa)):
                raise Exception('Interacting amino acids do not match the uniprot sequence!')

            uniprot_model.interacting_aa_1 = ','.join([str(uniprot_num) for uniprot_num in chain_1_interacting_uninum])
            uniprot_model.interacting_aa_2 = ','.join([str(uniprot_num) for uniprot_num in chain_2_interacting_uninum])

            # Get interacting amino acids and interface area
            analyze_structure_object = analyze_structure.AnalyzeStructure(
                self.unique_temp_folder + 'modeller/',
                self.unique_temp_folder + 'analyze_structure/',
                pdbFile_wt, modeller_chains, None, self.log)
            interface_area = analyze_structure_object.get_interface_area()
            uniprot_model.interface_area_hydrophobic = interface_area[0]
            uniprot_model.interface_area_hydrophilic = interface_area[1]
            uniprot_model.interface_area_total = interface_area[2]

        # Values common for single domains and interactions
        model_errors =', '.join(model_errors)
        if model_errors != '':
            uniprot_model.model_errors = model_errors

        return uniprot_model

#    def prepareInput(self, pdbCode, chain_ids, domains, sequences, alignments, save_path):
#
#        pdb = PDBTemplate(self.pdb_path, pdbCode, chain_ids, domains, save_path, self.unique_temp_folder, self.log)
#        hetatm_positions, het_flag, chain_ids_pdb_order, chain_sequence, chain_numbering_extended = pdb.extract()
#
#        switch_chain = False
#        if chain_ids != chain_ids_pdb_order:
#            switch_chain = True
#            sequence_new_order = list()
#            alignments_new_order = list()
#            for item in chain_ids_pdb_order:
#                i = chain_ids.index(item)
#                sequence_new_order.append(sequences[i])
#                alignments_new_order.append(alignments[i])
#            sequences = sequence_new_order
#            alignments = alignments_new_order
#            chain_ids = chain_ids_pdb_order
#
#        # to get the range of the domain make sure to call extract before get_chain_numbering!!
#        chain_numbering_w_hetatms = dict()
#        for chain_id in chain_ids:
#            chain_numbering_w_hetatms[chain_id] = pdb.get_chain_numbering(chain_id, include_hetatms=True)
#
#        HETATMsInChain_SEQnumbering = [{}] * len(chain_ids)
#        for i, chain_id in enumerate(chain_ids):
#            HETATMsInChain_SEQnumbering[i][chain_id] = list()
#            for item in hetatm_positions[chain_id]:
#                try:
#                    HETATMsInChain_SEQnumbering[i][chain_id].append(chain_numbering_w_hetatms[chain_id].index(item))
#                except ValueError:
#                    HETATMsInChain_SEQnumbering[i][chain_id].append(10001) # add one add the end that will be outside of the sequence
#
#        return sequences, alignments, chain_ids, switch_chain, het_flag, HETATMsInChain_SEQnumbering


    def _write_to_pir_alignment(self, pir_alignment_filehandle, seq_type, seq_name, seq):
        pir_alignment_filehandle.write('>P1;' + seq_name + '\n')
        pir_alignment_filehandle.write(seq_type + ':' + seq_name + ':.:.:.:.::::\n')
        pir_alignment_filehandle.write(seq + '*')
        pir_alignment_filehandle.write('\n\n')


    def _get_model(self, pir_alignment_filename, target_ids, template_ids, save_path):
        """
        """
        modeller_target_id = '_'.join(target_ids)

        if len(template_ids) == 1:
            # If you're only modelling one domain
            modeller_template_id = template_ids[0]
        elif len(template_ids) == 2:
            if template_ids[0][-1] == template_ids[1][-1]:
                # Modelling two domains using the same chain in the pdb
                modeller_template_id = template_ids[0]
            else:
                # Modelling two domains using two chains in the pdb
                # if more than two chains are used this has to be adjusted
                modeller_template_id = template_ids[0] + template_ids[1][-1]

        modeller_path = self.unique_temp_folder + 'modeller/'
#        os.chdir(modeller_path)
        modeller = call_modeller.modeller(
            [pir_alignment_filename],
            modeller_target_id,
            modeller_template_id,
            save_path, # path_to_pdb_for_modeller
            self.unique_temp_folder, # path to folders with executables
            self.log,
            self.modeller_runs,
            self.subprocess_ids,
            loopRefinement=True)
        normDOPE, pdbFile, knotted = modeller.run()

        # If there is only one chain in the pdb, label that chain 'A'
        io = PDBIO()
        structure = get_pdb_structure(modeller_path + pdbFile)
        model = structure[0]
        chains = model.child_list
        self.log.debug(', '.join(['chain id: %s' % chain.id for chain in chains]))
        new_chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        for i in range(len(chains)):
            chains[i].id = new_chains[i]
        self.log.debug(', '.join(['chain id: %s' % chain.id for chain in chains]))
        io.set_structure(structure)
        io.save(modeller_path + pdbFile)
#        os.chdir(self.unique_temp_folder) # go back to the working directory

        return normDOPE, pdbFile, knotted