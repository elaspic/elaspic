# -*- coding: utf-8 -*-

import numpy as np
import subprocess

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBIO

import helper_functions as hf
import pdb_template
import errors
import call_modeller
import call_tcoffee
import sql_db
import analyze_structure


class GetModel(object):
    """
    Attributes
    ----------
    """
    def __init__(
            self, global_temp_path, temp_path, unique, pdb_path,
            db, logger, n_cores, modeller_runs):
        """
        """
        self.global_temp_path = global_temp_path
        self.temp_path = temp_path
        self.unique = unique
        self.unique_temp_folder = temp_path + unique + '/'
        self.pdb_path = pdb_path
        self.db = db
        self.logger = logger
        self.modeller_runs = modeller_runs
        self.n_cores = n_cores


    def __call__(self, d):
        """
        """
        if isinstance(d, sql_db.UniprotDomain):
            return self.get_domain_model(d)
        elif isinstance(d, sql_db.UniprotDomainPair):
            return self.get_domain_pair_model(d)
        else:
            raise Exception("Wrong object type for parameter 'd'")


    def get_domain_model(self, d):
        """
        """
        self.logger.debug(
            'Aligning: {}/{}*{}:{}'
            .format(d.uniprot_id, d.pdbfam_name, d.template.domain_def.replace(':', '-'), d.template.cath_id))

        alignmnets, alignment_filenames, norm_dope_wt, pdb_filename_wt, knotted, model_errors = \
            self.perform_alignments_and_modelling(
                [d.uniprot_id], [d.uniprot_sequence.uniprot_sequence], [d.template.domain_def],
                d.template.domain.pdb_id, [d.template.domain.pdb_chain], [d.template.domain.pdb_domain_def], d.path_to_data)

        structure = hf.get_pdb_structure(self.unique_temp_folder + 'modeller/' + pdb_filename_wt)
        model = structure[0]

        ###
        uniprot_model = sql_db.UniprotDomainModel()
        uniprot_model.uniprot_domain_id = d.uniprot_domain_id
        uniprot_model.model_filename = pdb_filename_wt
        uniprot_model.norm_dope = norm_dope_wt
        uniprot_model.chain = 'A'
        uniprot_model.alignment_filename = alignment_filenames[0]
        # Run the homology model through msms and get dataframes with all the
        # per atom and per residue SASA values
        analyze_structure_object = analyze_structure.AnalyzeStructure(
            self.unique_temp_folder + 'modeller/',
            self.unique_temp_folder + 'analyze_structure/',
            pdb_filename_wt, [uniprot_model.chain], None, self.logger)
        (seasa_by_chain_together, seasa_by_chain_separately,
        seasa_by_residue_together, seasa_by_residue_separately) = analyze_structure_object.get_seasa()
        # Get SASA only for amino acids in the chain of interest
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
            self.logger.error('msms score length mismatch')
            model_errors.append('msms score length mismatch')

        # Values common for single domains and interactions
        model_errors =', '.join(model_errors)
        if model_errors != '':
            uniprot_model.model_errors = model_errors
        return uniprot_model


    def get_domain_pair_model(self, d):
        """
        """
        self.logger.debug(
            'Aligning: {}/{}*{}:{}/{}*{}:{}/{}'
            .format(
                d.uniprot_domain_1.uniprot_id,
                d.uniprot_domain_1.pdbfam_name,
                d.uniprot_domain_1.template.domain_def.replace(':', '-'),
                d.template.cath_id_1,
                d.uniprot_domain_2.pdbfam_name,
                d.uniprot_domain_2.template.domain_def.replace(':', '-'),
                d.template.cath_id_2,
                d.uniprot_domain_2.uniprot_id))

        alignmnets, alignment_filenames, norm_dope_wt, pdb_filename_wt, knotted, model_errors = \
            self.perform_alignments_and_modelling(
                [d.uniprot_domain_1.uniprot_id, d.uniprot_domain_2.uniprot_id],
                [d.uniprot_domain_1.uniprot_sequence.uniprot_sequence,
                 d.uniprot_domain_2.uniprot_sequence.uniprot_sequence],
                [d.uniprot_domain_1.template.domain_def,
                 d.uniprot_domain_2.template.domain_def],
                d.template.domain_1.pdb_id,
                [d.template.domain_1.pdb_chain, d.template.domain_2.pdb_chain],
                [d.template.domain_1.pdb_domain_def, d.template.domain_2.pdb_domain_def],
                d.path_to_data)

        structure = hf.get_pdb_structure(self.unique_temp_folder + 'modeller/' + pdb_filename_wt)
        model = structure[0]

        ###
        uniprot_model = sql_db.UniprotDomainPairModel()
        uniprot_model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
        uniprot_model.model_filename = pdb_filename_wt
        uniprot_model.norm_dope = norm_dope_wt
        uniprot_model.chain_1 = 'A'
        uniprot_model.chain_2 = 'B'
        uniprot_model.alignment_filename_1 = alignment_filenames[0]
        uniprot_model.alignment_filename_2 = alignment_filenames[1]

        interactions_between_chains = (
            analyze_structure.get_interactions_between_chains(
                model, uniprot_model.chain_1, uniprot_model.chain_2, 5) )
        chain_1_interactions = list(set([key[:2] for key in interactions_between_chains.keys()]))
        chain_2_interactions = list(set([value[:2] for values in interactions_between_chains.values() for value in values]))
        if not chain_1_interactions and not chain_2_interactions:
            raise errors.ChainsNotInteracting(
                'Chains %s and %s are not interacting! chain_1_interactions: %s, chain_2_interactions: %s' %
                (uniprot_model.chain_1, uniprot_model.chain_2, ','.join(list(chain_1_interactions)), ','.join(list(chain_2_interactions)),))

        chain_1_interactions.sort(key=lambda x: x[0]) # Sort by residue for easy reading
        chain_2_interactions.sort(key=lambda x: x[0]) # Sort by residue for easy reading

        chain_1_interacting_resnum, chain_1_interacting_aa = self._get_unique_resnum_and_sequence(chain_1_interactions)
        chain_1_numbering = pdb_template.get_chain_numbering(model[uniprot_model.chain_1], return_extended=True)
        chain_2_interacting_resnum, chain_2_interacting_aa = self._get_unique_resnum_and_sequence(chain_2_interactions)
        chain_2_numbering = pdb_template.get_chain_numbering(model[uniprot_model.chain_2], return_extended=True)

        self.logger.debug('chain_1_interactions:')
        self.logger.debug(chain_1_interactions)
        self.logger.debug('chain_1_interacting_resnum:')
        self.logger.debug(chain_1_interacting_resnum)
        self.logger.debug('chain_1_interacting_aa:')
        self.logger.debug(chain_1_interacting_aa)
#            self.logger.debug('chain_1_numbering:')
#            self.logger.debug(chain_1_numbering)

        self.logger.debug('chain_2_interactions:')
        self.logger.debug(chain_2_interactions)
        self.logger.debug('chain_2_interacting_resnum:')
        self.logger.debug(chain_2_interacting_resnum)
        self.logger.debug('chain_2_interacting_aa:')
        self.logger.debug(chain_2_interacting_aa)
#            self.logger.debug('chain_2_numbering:')
#            self.logger.debug(chain_2_numbering)

        chain_1_interacting_uninum = [
            hf.decode_domain_def(d.uniprot_domain_1.template.domain_def)[0] + chain_1_numbering.index(resnum)
            for resnum in chain_1_interacting_resnum]
        chain_2_interacting_uninum = [
            hf.decode_domain_def(d.uniprot_domain_2.template.domain_def)[0] + chain_2_numbering.index(resnum)
            for resnum in chain_2_interacting_resnum]

        chain_1_interactions_uniprot = zip(chain_1_interacting_uninum, chain_1_interacting_aa)
        chain_1_interactions_uniprot.sort(key=lambda x: x[0])
        chain_1_interacting_uninum, chain_1_interacting_aa = zip(*chain_1_interactions_uniprot)
        chain_1_interacting_aa = ''.join(chain_1_interacting_aa)
        chain_1_interacting_aa_from_uniprot = ''
        for uninum, aa in chain_1_interactions_uniprot:
            uniprot_idx = uninum - 1
            chain_1_interacting_aa_from_uniprot += d.uniprot_domain_1.uniprot_sequence.uniprot_sequence[uniprot_idx]

        chain_2_interactions_uniprot = zip(chain_2_interacting_uninum, chain_2_interacting_aa)
        chain_2_interactions_uniprot.sort(key=lambda x: x[0])
        chain_2_interacting_uninum, chain_2_interacting_aa = zip(*chain_2_interactions_uniprot)
        chain_2_interacting_aa = ''.join(chain_2_interacting_aa)
        chain_2_interacting_aa_from_uniprot = ''
        for uninum, aa in chain_2_interactions_uniprot:
            uniprot_idx = uninum - 1
            chain_2_interacting_aa_from_uniprot += d.uniprot_domain_2.uniprot_sequence.uniprot_sequence[uniprot_idx]

        self.logger.debug('domain_def_1: {}'.format(d.uniprot_domain_1.template.domain_def))
        self.logger.debug('chain_1_interactions_uniprot: {}'.format(chain_1_interactions_uniprot))
        self.logger.debug('uniprot_sequence_1: {}'.format(d.uniprot_domain_1.uniprot_sequence.uniprot_sequence))
        self.logger.debug('chain_1_interacting_aa: {}'.format(chain_1_interacting_aa))
        self.logger.debug('chain_1_interacting_aa_from_uniprot: {}'.format(chain_1_interacting_aa_from_uniprot))

        self.logger.debug('domain_def_2: {}'.format(d.uniprot_domain_2.template.domain_def))
        self.logger.debug('chain_2_interactions_uniprot: {}'.format(chain_2_interactions_uniprot))
        self.logger.debug('uniprot_sequence_2: {}'.format(d.uniprot_domain_2.uniprot_sequence.uniprot_sequence))
        self.logger.debug('chain_2_interacting_aa: {}'.format(chain_2_interacting_aa))
        self.logger.debug('chain_2_interacting_aa_from_uniprot: {}'.format(chain_2_interacting_aa_from_uniprot))

        if ((len(chain_1_interacting_uninum) == 0 or
                len(chain_2_interacting_uninum) == 0)):
            raise errors.ChainsNotInteracting(
                'Chains %s and %s are not interacting! chain_1_interactions: %s, chain_2_interactions: %s' %
                (uniprot_model.chain_1, uniprot_model.chain_2, ','.join(list(chain_1_interactions)), ','.join(list(chain_2_interactions)),))

        if ((chain_1_interacting_aa_from_uniprot != chain_1_interacting_aa) or
                (chain_2_interacting_aa_from_uniprot != chain_2_interacting_aa)):
            raise Exception('Interacting amino acids do not match the uniprot sequence!')

        uniprot_model.interacting_aa_1 = ','.join([str(uniprot_num) for uniprot_num in chain_1_interacting_uninum])
        uniprot_model.interacting_aa_2 = ','.join([str(uniprot_num) for uniprot_num in chain_2_interacting_uninum])

        # Get interacting amino acids and interface area
        analyze_structure_object = analyze_structure.AnalyzeStructure(
            self.unique_temp_folder + 'modeller/',
            self.unique_temp_folder + 'analyze_structure/',
            pdb_filename_wt, [uniprot_model.chain_1, uniprot_model.chain_2], None, self.logger)
        interface_area = analyze_structure_object.get_interface_area()
        uniprot_model.interface_area_hydrophobic = interface_area[0] if not np.isnan(interface_area[0]) else None
        uniprot_model.interface_area_hydrophilic = interface_area[1] if not np.isnan(interface_area[1]) else None
        uniprot_model.interface_area_total = interface_area[2] if not np.isnan(interface_area[2]) else None

        # Values common for single domains and interactions
        model_errors =', '.join(model_errors)
        if model_errors != '':
            uniprot_model.model_errors = model_errors

        return uniprot_model


    def _get_unique_resnum_and_sequence(self, interactions):
        """
        """
        interacting_resnum = []
        interacting_aa = ''
        unique_resnums = set()
        for interaction in interactions:
            if interaction[0] not in unique_resnums:
                unique_resnums.add(interaction[0])
                interacting_resnum.append(interaction[0])
                interacting_aa += interaction[1]
        return interacting_resnum, interacting_aa


    def perform_alignments_and_modelling(
            self, uniprot_ids, uniprot_sequences, uniprot_domain_defs,
            pdb_id, pdb_chains, pdb_domain_defs, path_to_data):
        """
        """
        model_errors = []
        # Folder for storing files for export to output
        save_path = self.temp_path + path_to_data
        self.logger.debug('Template pdb: {}'.format(pdb_id))
        self.logger.debug('save path: {}'.format(save_path))
        uniprot_domain_seqrecs = self.get_uniprot_domain_seqrecs(uniprot_ids, uniprot_sequences, uniprot_domain_defs)
        pdb, pdb_domain_seqrecs = self.get_pdb_domain_seqrecs(pdb_id, pdb_chains, pdb_domain_defs)

        # Perform the alignments
        alignmnets = []
        alignment_filenames = []
        for uniprot_domain_seqrec, pdb_domain_seqrec in zip(uniprot_domain_seqrecs, pdb_domain_seqrecs):
            alignment, alignment_filename = self.perform_alignment(uniprot_domain_seqrec, pdb_domain_seqrec, 'expresso', path_to_data)
            alignmnets.append(alignment)
            alignment_filenames.append(alignment_filename)

        # Join alignments for different chains
        target_ids = [al[0].id for al in alignmnets]
        target_id = '_'.join(target_ids)
        target_seq = '/'.join([str(alignmnet[0].seq) for alignmnet in alignmnets])

        template_ids = [al[1].id for al in alignmnets]
        template_id = pdb_id + ''.join(pdb_chains)
        template_seq = '/'.join([str(alignmnet[1].seq) for alignmnet in alignmnets])

        # Add '.' in place of every heteroatom in the pdb
        hetatm_chains = []
        for chain in pdb.structure[0].child_list:
            if chain.id not in pdb_chains:
                hetatm_chains.append(chain)
        self.logger.debug('Number of extra chains: {}'.format(len(hetatm_chains)))
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
        norm_dope_wt, pdb_filename_wt, knotted = self._run_modeller(pir_alignment_filename, target_ids, template_ids, self.unique_temp_folder)
        self.logger.debug('model pdb file: %s, knotted: %s' % (pdb_filename_wt, knotted,))
        if knotted:
            model_errors.append('knotted')

        # Save another copy of the modelled structure in the tmp export folder
        subprocess.check_call('cp ' + self.unique_temp_folder + 'modeller/' + pdb_filename_wt + ' ' + save_path + pdb_filename_wt, shell=True)

        return alignmnets, alignment_filenames, norm_dope_wt, pdb_filename_wt, knotted, model_errors


    def get_uniprot_domain_seqrecs(self, uniprot_ids, uniprot_sequences, uniprot_domain_defs):
        """
        """
        uniprot_domain_seqrecs = []
        for uniprot_id, uniprot_sequence, uniprot_domain_def in zip(
                uniprot_ids, uniprot_sequences, uniprot_domain_defs):
            uniprot_seqrec = SeqRecord(seq=Seq(uniprot_sequence), id=uniprot_id)
            domain_def = hf.decode_domain_def(uniprot_domain_def, merge=True, return_string=False)
            uniprot_domain_seqrecs.append( uniprot_seqrec[domain_def[0]-1:domain_def[1]] )
        return uniprot_domain_seqrecs


    def get_pdb_domain_seqrecs(self, pdb_id, pdb_chains, pdb_domain_defs):
        """
        """
        pdb_domain_defs = [
            hf.decode_domain_def(pdb_domain_def, merge=False, return_string=True)
            for pdb_domain_def in pdb_domain_defs]
        pdb = pdb_template.PDBTemplate(
            self.pdb_path, pdb_id, pdb_chains, pdb_domain_defs,
            self.unique_temp_folder, self.unique_temp_folder, self.logger)
        pdb.extract()
        pdb_domain_seqrecs = []
        for pdb_chain in pdb_chains:
            chain_sequence = pdb.chain_sequence_dict[pdb_chain]
            pdb_domain_seqrecs.append(SeqRecord(seq=Seq(chain_sequence), id=pdb_id+pdb_chain))
        return pdb, pdb_domain_seqrecs


    def perform_alignment(self, uniprot_seqrecord, pdb_seqrecord, mode, path_to_data):
        """
        """
        # Write both sequences to a single file
        with open(self.unique_temp_folder + 'seqfiles.fasta', 'w') as seqfiles:
            SeqIO.write([uniprot_seqrecord, pdb_seqrecord], seqfiles, 'fasta')
        # Perform the alignment
        sequence_ids = [uniprot_seqrecord.id, pdb_seqrecord.id]
        tcoffee = call_tcoffee.tcoffee_alignment(
            self.global_temp_path,
            self.unique_temp_folder,
            [self.unique_temp_folder + 'seqfiles.fasta'],
            sequence_ids,
            self.n_cores,
            self.pdb_path,
            mode,
            self.logger)
        alignments = tcoffee.align()
        assert len(alignments) == 1
        alignment = alignments[0]
        # Save the alignment
        alignment_filename = alignment[0].id + '_' + alignment[1].id + '.aln'
        try:
            AlignIO.write(alignment, self.unique_temp_folder + 'tcoffee/' + alignment_filename, 'clustal')
        except IndexError as e:
            raise errors.EmptyPDBSequenceError('{}: {}'.format(type(e), e))
        temp_save_path = self.temp_path + path_to_data
        subprocess.check_call('mkdir -p ' + temp_save_path, shell=True)
        subprocess.check_call('cp ' + self.unique_temp_folder + 'tcoffee/' + alignment_filename +
                                ' ' + temp_save_path + alignment_filename, shell=True)
        return alignment, alignment_filename


    def analyze_alignment(self, alignment, pdb_contact_idxs=[]):
        """
        """
        pdb_aa_idx = -1
        sequence_1_length = 0
        sequence_1_identity = 0
        sequence_1_coverage = 0
        interface_1_identity = 0
        interface_1_coverage = 0

        for aa_1, aa_2 in zip(alignment):
            is_interface = False
            # Check if the amino acid falls in a gap
            if aa_1 == '-':
                continue
            sequence_1_length += 1
            # Check if the template is in a gap
            if aa_2 == '-':
                continue
            pdb_aa_idx += 1
            if pdb_aa_idx in pdb_contact_idxs:
                is_interface = True # This is an interface amino acid
            sequence_1_coverage += 1 # Count as coverage
            if is_interface:
                interface_1_coverage += 1 # Count as coverage
            # Check if the template is identical
            if aa_1 != aa_2:
                continue
            sequence_1_identity += 1 # Count as identity
            if is_interface:
                interface_1_identity += 1 # Count as identity

        identity = sequence_1_identity / float(sequence_1_length) * 100
        coverage = sequence_1_coverage / float(sequence_1_length) * 100
        if_identity = sequence_1_identity / float(len(pdb_contact_idxs)) * 100 if pdb_contact_idxs else None
        if_coverage = interface_1_coverage / float(len(pdb_contact_idxs)) * 100 if pdb_contact_idxs else None

        return identity, coverage, if_identity, if_coverage


    def score_alignment(self, identity, coverage, alpha=0.95):
        """
        """
        return alpha * (identity) * (coverage) + (1.0 - alpha) * (coverage)


    def _write_to_pir_alignment(self, pir_alignment_filehandle, seq_type, seq_name, seq):
        pir_alignment_filehandle.write('>P1;' + seq_name + '\n')
        pir_alignment_filehandle.write(seq_type + ':' + seq_name + ':.:.:.:.::::\n')
        pir_alignment_filehandle.write(seq + '*')
        pir_alignment_filehandle.write('\n\n')


    def _run_modeller(self, pir_alignment_filename, target_ids, template_ids, save_path):
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
        modeller = call_modeller.modeller(
            [pir_alignment_filename],
            modeller_target_id,
            modeller_template_id,
            save_path, # path_to_pdb_for_modeller
            self.unique_temp_folder, # path to folders with executables
            self.logger,
            self.modeller_runs,
            loopRefinement=True)
        with hf.switch_paths(modeller_path):
            normDOPE, pdbFile, knotted = modeller.run()

        # If there is only one chain in the pdb, label that chain 'A'
        io = PDBIO()
        structure = hf.get_pdb_structure(modeller_path + pdbFile)
        model = structure[0]
        chains = model.child_list
        self.logger.debug(', '.join(['chain id: %s' % chain.id for chain in chains]))
        new_chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        for i in range(len(chains)):
            chains[i].id = new_chains[i]
        self.logger.debug(', '.join(['chain id: %s' % chain.id for chain in chains]))
        io.set_structure(structure)
        io.save(modeller_path + pdbFile)

        return normDOPE, pdbFile, knotted





