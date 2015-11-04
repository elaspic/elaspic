# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 01:08:59 2015

@author: strokach
"""

from __future__ import absolute_import
from __future__ import unicode_literals
from builtins import zip
from builtins import range
from builtins import object

import os.path as op
import logging

import numpy as np
import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import conf, helper, errors, structure_tools, structure_analysis, database, model
configs = conf.Configs()

logger = logging.getLogger(__name__)



# %%
class GetModel(object):
    """
    Attributes
    ----------
    """
    def __init__(self, unique_temp_folder, db):
        """
        """
        self.global_temp_path = configs['global_temp_path']
        self.temp_archive_path = configs['temp_archive_path']
        self.unique_temp_folder = unique_temp_folder
        self.pdb_path = configs['pdb_path']
        self.db = db
        self.modeller_runs = configs['modeller_runs']
        self.n_cores = configs['n_cores']
        self._prepare_temp_folder()


    def _prepare_temp_folder(self):
        """
        Create temporary folders for `t_coffee` and `modeller`,
        required for making the homology models.
        """
        # t_coffee
        if not op.isdir(op.join(self.unique_temp_folder + 'tcoffee')):
            mkdir_command = (
                "mkdir -p '{0}/tcoffee' && "
                "mkdir -p '{0}/tcoffee/tmp' && "
                "mkdir -p '{0}/tcoffee/lck' && "
                "mkdir -p '{0}/tcoffee/cache' "
            ).format(self.unique_temp_folder)
            subprocess.check_call(mkdir_command, shell=True)

        # modeller
        if not op.isdir(op.join(self.unique_temp_folder + 'modeller')):
            mkdir_command = "mkdir -p '{}/modeller' ".format(self.unique_temp_folder)
            subprocess.check_call(mkdir_command, shell=True)


    def __call__(self, d):
        """
        """
        if isinstance(d, database.UniprotDomain):
            self.get_domain_model(d)
        elif isinstance(d, database.UniprotDomainPair):
            self.get_domain_pair_model(d)
        else:
            raise Exception("Wrong object type for parameter 'd'")


    def get_domain_model(self, d):
        """
        """
        logger.debug(
            'Aligning: {}/{}*{}:{}'
            .format(d.uniprot_id, d.pdbfam_name, d.template.domain_def.replace(':', '-'), d.template.cath_id))

        alignmnets, alignment_filenames, norm_dope_wt, pdb_filename_wt, knotted, model_errors, domain_def_offsets = (
            self.perform_alignments_and_modelling(
                [d.uniprot_id], [d.uniprot_sequence.uniprot_sequence], [d.template.domain_def],
                d.template.domain.pdb_id, [d.template.domain.pdb_chain],
                [d.template.domain.pdb_domain_def], d.path_to_data)
        )
        logger.debug('Finished performing alignments!')
        structure = structure_tools.get_pdb_structure(op.join(self.unique_temp_folder, 'modeller', pdb_filename_wt))
        model = structure[0]

        ###
        if d.template.model == None:
            d.template.model = database.UniprotDomainModel()
            d.template.model.uniprot_domain_id = d.uniprot_domain_id
        d.template.model.model_filename = pdb_filename_wt
        d.template.model.norm_dope = norm_dope_wt
        d.template.model.chain = 'A'
        d.template.model.alignment_filename = alignment_filenames[0]
        # Run the homology model through msms and get dataframes with all the
        # per atom and per residue SASA values
        analyze_structure_object = structure_analysis.AnalyzeStructure(
            self.unique_temp_folder + 'modeller/',
            self.unique_temp_folder + 'analyze_structure/',
            pdb_filename_wt, [d.template.model.chain], None)
        (seasa_by_chain_together, seasa_by_chain_separately,
        seasa_by_residue_together, seasa_by_residue_separately) = analyze_structure_object.get_seasa()
        # Get SASA only for amino acids in the chain of interest
        msms_length_mismatch = False
        rel_sasa_for_each_residue = []
        for res in model[d.template.model.chain]:
            if res.resname in structure.amino_acids:
                rel_sasa = seasa_by_residue_separately[
                    (seasa_by_residue_separately['pdb_chain']==d.template.model.chain) &
                    (seasa_by_residue_separately['res_name']==res.resname) &
                    (seasa_by_residue_separately['res_num']==(str(res.id[1]) + res.id[2].strip()))].iloc[0]['rel_sasa']
                if rel_sasa is None or rel_sasa is np.nan:
                    msms_length_mismatch = True
                rel_sasa_for_each_residue.append(rel_sasa)
        d.template.model.sasa_score = (
            ','.join(['{:.2f}'.format(x) for x in model.relative_sasa_scores[model.chain_ids]])
        )

        d.template.model.model_domain_def = self._truncate_domain_defs(
            d.template.domain_def,
            domain_def_offsets[0])

        if msms_length_mismatch:
            logger.error('msms score length mismatch')
            model_errors.append('msms score length mismatch')

        # Values common for single domains and interactions
        model_errors =', '.join('{}'.format(e) for e in model.errors)
        if model_errors != '':
            d.template.model.model_errors = model_errors



    def get_domain_pair_model(self, d):
        """
        """
        logger.debug(
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

        alignmnets, alignment_filenames, norm_dope_wt, pdb_filename_wt, knotted, model_errors, domain_def_offsets = \
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

        structure = structure_tools.get_pdb_structure(self.unique_temp_folder + 'modeller/' + pdb_filename_wt)
        model = structure[0]

        model_domain_def_1 = self._truncate_domain_defs(
            d.uniprot_domain_1.template.domain_def,
            domain_def_offsets[0])

        model_domain_def_2 = self._truncate_domain_defs(
            d.uniprot_domain_2.template.domain_def,
            domain_def_offsets[1])


        ###
        if d.template.model == None:
            d.template.model = database.UniprotDomainPairModel()
            d.template.model.uniprot_domain_pair_id = d.uniprot_domain_pair_id
        d.template.model.model_filename = pdb_filename_wt
        d.template.model.norm_dope = norm_dope_wt
        d.template.model.chain_1 = 'A'
        d.template.model.chain_2 = 'B'
        d.template.model.alignment_filename_1 = alignment_filenames[0]
        d.template.model.alignment_filename_2 = alignment_filenames[1]

        interactions_between_chains = (
            structure_analysis.get_interactions_between_chains(
                model, d.template.model.chain_1, d.template.model.chain_2, 5) )
        chain_1_interactions = list(set([key[:2] for key in list(interactions_between_chains.keys())]))
        chain_2_interactions = list(set([value[:2] for values in list(interactions_between_chains.values()) for value in values]))
        if not chain_1_interactions and not chain_2_interactions:
            raise errors.ChainsNotInteractingError(
                'Chains %s and %s are not interacting! chain_1_interactions: %s, chain_2_interactions: %s' %
                (d.template.model.chain_1, d.template.model.chain_2,
                 ','.join(list(chain_1_interactions)), ','.join(list(chain_2_interactions)),))

        chain_1_interactions.sort(key=lambda x: int(''.join([c for c in x[0] if c.isdigit()]))) # Sort by residue for easy reading
        chain_2_interactions.sort(key=lambda x: int(''.join([c for c in x[0] if c.isdigit()]))) # Sort by residue for easy reading

        chain_1_interacting_resnum, chain_1_interacting_aa = self._get_unique_resnum_and_sequence(chain_1_interactions)
        __, chain_1_numbering = structure.get_chain_sequence_and_numbering(model[d.template.model.chain_1])
        chain_2_interacting_resnum, chain_2_interacting_aa = self._get_unique_resnum_and_sequence(chain_2_interactions)
        __, chain_2_numbering = structure.get_chain_sequence_and_numbering(model[d.template.model.chain_2])

        logger.debug('chain_1_interactions:')
        logger.debug(chain_1_interactions)
        logger.debug('chain_1_interacting_resnum:')
        logger.debug(chain_1_interacting_resnum)
        logger.debug('chain_1_interacting_aa:')
        logger.debug(chain_1_interacting_aa)
#            logger.debug('chain_1_numbering:')
#            logger.debug(chain_1_numbering)

        logger.debug('chain_2_interactions:')
        logger.debug(chain_2_interactions)
        logger.debug('chain_2_interacting_resnum:')
        logger.debug(chain_2_interacting_resnum)
        logger.debug('chain_2_interacting_aa:')
        logger.debug(chain_2_interacting_aa)
#            logger.debug('chain_2_numbering:')
#            logger.debug(chain_2_numbering)

        model_domain_1_start = helper.decode_domain_def(model_domain_def_1)[0]
        chain_1_interacting_uninum = [
            model_domain_1_start + chain_1_numbering.index(resnum)
            for resnum in chain_1_interacting_resnum]

        model_domain_2_start = helper.decode_domain_def(model_domain_def_2)[0]
        chain_2_interacting_uninum = [
            model_domain_2_start + chain_2_numbering.index(resnum)
            for resnum in chain_2_interacting_resnum]

        chain_1_interactions_uniprot = list(zip(chain_1_interacting_uninum, chain_1_interacting_aa))
        chain_1_interactions_uniprot.sort(key=lambda x: x[0])
        chain_1_interacting_uninum, chain_1_interacting_aa = list(zip(*chain_1_interactions_uniprot))
        chain_1_interacting_aa = ''.join(chain_1_interacting_aa)
        chain_1_interacting_aa_from_uniprot = ''
        for uninum, aa in chain_1_interactions_uniprot:
            uniprot_idx = uninum - 1
            chain_1_interacting_aa_from_uniprot += d.uniprot_domain_1.uniprot_sequence.uniprot_sequence[uniprot_idx]

        chain_2_interactions_uniprot = list(zip(chain_2_interacting_uninum, chain_2_interacting_aa))
        chain_2_interactions_uniprot.sort(key=lambda x: x[0])
        chain_2_interacting_uninum, chain_2_interacting_aa = list(zip(*chain_2_interactions_uniprot))
        chain_2_interacting_aa = ''.join(chain_2_interacting_aa)
        chain_2_interacting_aa_from_uniprot = ''
        for uninum, aa in chain_2_interactions_uniprot:
            uniprot_idx = uninum - 1
            chain_2_interacting_aa_from_uniprot += d.uniprot_domain_2.uniprot_sequence.uniprot_sequence[uniprot_idx]

        logger.debug('domain_def_1: {}'.format(d.uniprot_domain_1.template.domain_def))
        logger.debug('model_domain_def_1: {}'.format(model_domain_def_1))
        logger.debug('chain_1_interactions_uniprot: {}'.format(chain_1_interactions_uniprot))
        logger.debug('uniprot_sequence_1: {}'.format(d.uniprot_domain_1.uniprot_sequence.uniprot_sequence))
        logger.debug('chain_1_interacting_aa: {}'.format(chain_1_interacting_aa))
        logger.debug('chain_1_interacting_aa_from_uniprot: {}'.format(chain_1_interacting_aa_from_uniprot))

        logger.debug('domain_def_2: {}'.format(d.uniprot_domain_2.template.domain_def))
        logger.debug('model_domain_def_2: {}'.format(model_domain_def_2))
        logger.debug('chain_2_interactions_uniprot: {}'.format(chain_2_interactions_uniprot))
        logger.debug('uniprot_sequence_2: {}'.format(d.uniprot_domain_2.uniprot_sequence.uniprot_sequence))
        logger.debug('chain_2_interacting_aa: {}'.format(chain_2_interacting_aa))
        logger.debug('chain_2_interacting_aa_from_uniprot: {}'.format(chain_2_interacting_aa_from_uniprot))

        if ((len(chain_1_interacting_uninum) == 0 or
                len(chain_2_interacting_uninum) == 0)):
            raise errors.ChainsNotInteracting(
                'Chains %s and %s are not interacting! chain_1_interactions: %s, chain_2_interactions: %s' %
                (d.template.model.chain_1, d.template.model.chain_2,
                ','.join(list(chain_1_interactions)), ','.join(list(chain_2_interactions)),))

        if ((chain_1_interacting_aa_from_uniprot != chain_1_interacting_aa) or
                (chain_2_interacting_aa_from_uniprot != chain_2_interacting_aa)):
            raise Exception('Interacting amino acids do not match the uniprot sequence!')

        d.template.model.interacting_aa_1 = ','.join([str(uniprot_num) for uniprot_num in chain_1_interacting_uninum])
        d.template.model.interacting_aa_2 = ','.join([str(uniprot_num) for uniprot_num in chain_2_interacting_uninum])

        # Get interacting amino acids and interface area
        analyze_structure_object = structure_analysis.AnalyzeStructure(
            self.unique_temp_folder + 'modeller/',
            self.unique_temp_folder + 'analyze_structure/',
            pdb_filename_wt, [d.template.model.chain_1, d.template.model.chain_2], None)
        interface_area = analyze_structure_object.get_interface_area()
        logger.debug('interface_area: {}'.format(interface_area))
        d.template.model.interface_area_hydrophobic = interface_area[0]
        d.template.model.interface_area_hydrophilic = interface_area[1]
        d.template.model.interface_area_total = interface_area[2]

        # Save model_domain_defs, which might be truncated compared to uniprot_domain_template domain defs
        d.template.model.model_domain_def_1 = model_domain_def_1
        d.template.model.model_domain_def_2 = model_domain_def_2

        # Values common for single domains and interactions
        model_errors =', '.join(model_errors)
        if model_errors != '':
            d.template.model.model_errors = model_errors

        # return uniprot_model


    def _truncate_domain_defs(self, domain_def, domain_def_offset):
        n_gaps_start, n_gaps_end = domain_def_offset
        domain_def_new = (
            [str(int(domain_def.split(':')[0]) + n_gaps_start)] +
            domain_def.split(':')[1:-1] +
            [str(int(domain_def.split(':')[-1]) - n_gaps_end)])
        domain_def_new = ':'.join(domain_def_new)
        return domain_def_new


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
        save_path = self.temp_archive_path + path_to_data
        logger.debug('Template pdb: {}'.format(pdb_id))
        logger.debug('save path: {}'.format(save_path))
        uniprot_domain_seqrecs = self.get_uniprot_domain_seqrecs(uniprot_ids, uniprot_sequences, uniprot_domain_defs)
        pdb, pdb_domain_seqrecs = self.get_pdb_domain_seqrecs(pdb_id, pdb_chains, pdb_domain_defs)

        # Perform the alignments
        logger.debug('Performing alignments...')
        alignmnets = []
        alignment_filenames = []
        domain_def_offsets = []
        for i in range(len(uniprot_domain_seqrecs)):
            alignment, alignment_filename = self.perform_alignment(
                uniprot_domain_seqrecs[i], pdb_domain_seqrecs[i], '3dcoffee', path_to_data)
            domain_def_offset = self._get_alignment_overhangs(alignment)
            if any(domain_def_offset):
                logger.debug('Shortening uniprot domain sequence because the alignment had large overhangs...')
                cut_from_start = domain_def_offset[0] if domain_def_offset[0] else None
                cut_from_end = -domain_def_offset[1] if domain_def_offset[1] else None
                uniprot_domain_seqrecs[i] = uniprot_domain_seqrecs[i][cut_from_start:cut_from_end]
                alignment, alignment_filename = self.perform_alignment(
                    uniprot_domain_seqrecs[i], pdb_domain_seqrecs[i], '3dcoffee', path_to_data)
            alignmnets.append(alignment)
            alignment_filenames.append(alignment_filename)
            domain_def_offsets.append(domain_def_offset)

        # Join alignments for different chains
        logger.debug('Joining alignments for different chains...')
        target_ids = [al[0].id for al in alignmnets]
        target_id = '_'.join(target_ids)
        target_seq = '/'.join([str(alignmnet[0].seq) for alignmnet in alignmnets])

        template_ids = [al[1].id for al in alignmnets]
        template_id = pdb_id + ''.join(pdb_chains)
        template_seq = '/'.join([str(alignmnet[1].seq) for alignmnet in alignmnets])

        # Add '.' in place of every heteroatom in the pdb
        logger.debug('Adding hetatms to alignment...')
        hetatm_chains = []
        for chain in pdb.structure[0].child_list:
            if chain.id not in pdb_chains:
                hetatm_chains.append(chain)
        logger.debug('Number of extra chains: {}'.format(len(hetatm_chains)))
        if len(hetatm_chains) > 1:
            raise Exception('Too many extra chains!')
        if len(hetatm_chains) == 1:
            number_of_hetatms_excluding_water = len([res for res in hetatm_chains[0] if res.id[0] != 'W'])
            target_seq += '/' + '.' * number_of_hetatms_excluding_water
            template_seq += '/' + '.'* number_of_hetatms_excluding_water

        # Write the alignment file for modeller
        pir_alignment_filename = self.unique_temp_folder + 'modeller/outFile_wildtype'
        with open(pir_alignment_filename, 'w') as pir_alignment_filehandle:
            model.write_to_pir_alignment(pir_alignment_filehandle, 'sequence', target_id, target_seq)
            model.write_to_pir_alignment(pir_alignment_filehandle, 'structure', template_id, template_seq)

        # Make the homology model and check if it is knotted
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
                # If more than two chains are used this has to be adjusted
                modeller_template_id = template_ids[0] + template_ids[1][-1]

        norm_dope_wt, pdb_filename_wt, knotted = model.run_modeller(
            pir_alignment_filename, modeller_target_id, modeller_template_id
        )
        logger.debug('model pdb file: %s, knotted: %s' % (pdb_filename_wt, knotted,))
        if knotted:
            model_errors.append('knotted')

        # Rename the pdb file in order to prevent overwrites for cases where you wish
        # to use multiple models (for example when making the training set)
        pdb_filename_wt_new = '_'.join(uniprot_ids) + '_' + pdb_id + ''.join(pdb_chains) + '.pdb'
        subprocess.check_call("cp -f '{}' '{}'".format(
            self.unique_temp_folder + 'modeller/' + pdb_filename_wt,
            self.unique_temp_folder + 'modeller/' + pdb_filename_wt_new), shell=True)
        # Save another copy of the modelled structure in the tmp export folder
        subprocess.check_call("cp -f '{}' '{}'".format(
            self.unique_temp_folder + 'modeller/' + pdb_filename_wt_new,
            save_path + pdb_filename_wt_new), shell=True)

        return alignmnets, alignment_filenames, norm_dope_wt, pdb_filename_wt_new, knotted, model_errors, domain_def_offsets


    def _get_alignment_overhangs(self, aln):
        """ Remove gap overhangs from the alignments.
        There are cases where no template sequence is availible for a big chunk
        of the protein. Return the number of amino acids that should be removed
        from the start and end of the query sequence in order to match the template.
        """
        n_gaps_start = 0
        n_gaps_end = 0
        for aa_query, aa_template in zip(*aln):
            if aa_query != '-' and aa_template == '-':
                n_gaps_start += 1
            else:
                break
        for aa_query, aa_template in reversed(list(zip(*aln))):
            if aa_query != '-' and aa_template == '-':
                n_gaps_end += 1
            else:
                break
        return n_gaps_start, n_gaps_end


    def get_uniprot_domain_seqrecs(self, uniprot_ids, uniprot_sequences, uniprot_domain_defs):
        """
        """
        uniprot_domain_seqrecs = []
        for uniprot_id, uniprot_sequence, uniprot_domain_def in zip(
                uniprot_ids, uniprot_sequences, uniprot_domain_defs):
            uniprot_seqrec = SeqRecord(id=uniprot_id, seq=Seq(uniprot_sequence))
            domain_def = helper.decode_domain_def(uniprot_domain_def, merge=True, return_string=False)
            uniprot_domain_seqrecs.append( uniprot_seqrec[domain_def[0]-1:domain_def[1]] )
        return uniprot_domain_seqrecs


    def get_pdb_domain_seqrecs(self, pdb_id, pdb_chains, pdb_domain_defs):
        """
        """
        logger.debug('Obtaining protein sequence for pdb: {}, chains: {}, domain_defs: {}'
            .format(pdb_id, pdb_chains, pdb_domain_defs))
        pdb = structure_tools.PDBTemplate(
            self.pdb_path, pdb_id, pdb_chains, pdb_domain_defs,
            self.unique_temp_folder, self.unique_temp_folder)
        pdb.extract()
        pdb.save_structure()
        pdb.save_sequences()
        pdb_domain_seqrecs = []
        for pdb_chain in pdb_chains:
            chain_sequence = pdb.chain_sequence_dict[pdb_chain]
            pdb_domain_seqrecs.append(SeqRecord(seq=Seq(chain_sequence), id=pdb_id+pdb_chain))
        return pdb, pdb_domain_seqrecs

