# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:25:45 2013

@author: niklas
"""
import os
import psutil
import time
import subprocess
import random

import numpy as np

import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sql_db
import errors
import call_tcoffee
import pdb_template
import helper_functions as hf


class GetTemplate():
    """
    Parent class holding functions for finding the correct template given a
    uniprot sequence
    """
    def __init__(
            self, global_temp_path, tmpPath, unique, pdb_path, db, log, n_cores,
            provean_temp_path, refine=False):
        """
        input:
        tmpPath             type 'str'
        unique              type 'str'
        pdb_path             type 'str'
        savePDB             type 'str'
        saveAlignments      type 'str'
        """
        self.global_temp_path = global_temp_path
        self.tmpPath = tmpPath
        self.unique = unique + '/'
        self.unique_temp_folder = tmpPath + unique + '/'
        self.pdb_path = pdb_path
        self.db = db
        self.log = log
        self.refine = refine
        self.n_cores = n_cores
        self.provean_temp_path = provean_temp_path
        self.alignment_identity_cutoff = 0.25
        self.bad_pdbs = ['3C4D', '3LB7', '3NSV', '2NP8', '2WN0']
        self.possible_template_expansion_errors = (
            errors.PDBDomainDefsError,
            errors.PDBChainError,
            errors.NoPDBFound, #PDBNotFoundError
            errors.EmptyPDBSequenceError,
            errors.PDBError,
            errors.TcoffeeBlastError,
            errors.TcoffeePDBidError,
            errors.TcoffeeError,
            Bio.PDB.PDBExceptions.PDBConstructionException,)


    def __call__(self, d):
        """
        """
        list_of_templates = self.run(d, None)
        if not list_of_templates:
            raise errors.NoTemplatesFound('Templates present in PDBfam were not usable')
        self.log.debug('Number of templates: %i' % len(list_of_templates))

        best_template = self.chose_best_template(list_of_templates)
        self.log.debug('The best template: ')
        self.log.debug(best_template)

        # Don't need refinement anymore, each template is refined already
        if self.refine:
            self.log.debug('Choosing the best templates withing a cluster...')
            list_of_templates = self.run(d, best_template)
            best_template = self.chose_best_template(list_of_templates)
            self.log.debug('The best template in the cluster is: ')
            self.log.debug(best_template)

        # Set up the paths and exporting the alignments
        if type(d) == sql_db.UniprotDomain:
            # Save one copy of the allignment for immediate usage in the tmp folder
            tmp_save_path = self.tmpPath + d.path_to_data
            subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
            subprocess.check_call('cp ' + self.unique_temp_folder + 'tcoffee/' + best_template.alignment_filename +
                                    ' ' + tmp_save_path + best_template.alignment_filename, shell=True)
#            (best_template.provean_supset_filename,
#            best_template.provean_supset_length) = self.build_provean_supporting_set(d, best_template)

        elif type(d) == sql_db.UniprotDomainPair:
            # Save one copy of the allignment for immediate usage in the tmp folder
            tmp_save_path = self.tmpPath + d.path_to_data
            subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
            subprocess.check_call('cp ' + self.unique_temp_folder + 'tcoffee/' + best_template.alignment_filename_1 +
                                    ' ' + tmp_save_path + best_template.alignment_filename_1, shell=True)
            subprocess.check_call('cp ' + self.unique_temp_folder + 'tcoffee/' + best_template.alignment_filename_2 +
                                    ' ' + tmp_save_path + best_template.alignment_filename_2, shell=True)
        return best_template


    def run(self, d, canonical_domain):
        """
        d: sql_db.UniprotDomain or sql_db.UniprotDomainPair object
        domain: sql_db.Domain or sql_db.DomainContact object
        """
        #######################################################################
        # Obtain a list of template domains and calculate some of their properties
        if isinstance(d, sql_db.UniprotDomain):
            # Add all partial combinations of superdomains to the set of superdomains
            self.log.debug('Finding templates for pfam domain {}...'.format(d.pfam_name))
            domain_list = self.db.get_domain([d.pfam_name])
            self._remove_bad_domains(domain_list)
            if not domain_list:
                self.log.debug('Did not find any structural templates for {}, trying subdomains...'.format(d.pfam_name))
                pfam_names = self._split_superdomains(d.pfam_name)
                self.log.debug('Pfam domain {} has subdomains {}'.format(d.pfam_name, str(pfam_names)))
                domain_list = self.db.get_domain(pfam_names, subdomains=True)
                self._remove_bad_domains(domain_list)
            if not domain_list:
                raise errors.NoStructuralTemplates('No templates found')

#            global_cdhit_clusters = set()
#            local_cdhit_clusters = set()
#            local_2_global = {}
#            for domain in domain_list:
#                # Domains might have different pfam names with the same cd-hit id
#                if domain.cdhit_cluster not in global_cdhit_clusters:
#                    # cdhit cluster for this domain is
#                    local_cdhit_clusters.add(domain.cdhit_cluster)
#                elif local_2_global.has_key(domain.cdhit_cluster):
#                    domain.cdhit_cluster = local_2_global[domain.cdhit_cluster]
#                    self.db.add_domain(domain)
#                else:
#                    for global_id in range(1000,2000):
#                        if global_id not in (global_cdhit_clusters | local_cdhit_clusters):
#                            break
#                    local_2_global[domain.cdhit_cluster] = global_id
#                    domain.cdhit_cluster = global_id
#                    self.db.add_domain(domain)
#                    local_cdhit_clusters.add(domain.cdhit_cluster)
#                global_cdhit_clusters.update(local_cdhit_clusters)
            # Get the length of the longest domain
            domain_list.sort(key=lambda d: d.cdhit_cluster_idx)
            max_domain_length = 0
            for domain in domain_list:
                domain_length = \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(sql_db.decode_domain(domain.pdb_domain_def, return_string=True)[1]) - \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(sql_db.decode_domain(domain.pdb_domain_def, return_string=True)[0]) + 1
                if domain_length > max_domain_length:
                    max_domain_length = domain_length

        elif isinstance(d, sql_db.UniprotDomainPair):
            self.log.debug('Finding templates for pfam domains {}, {}...'
                .format(d.uniprot_domain_1.pfam_name, d.uniprot_domain_2.pfam_name))
            domain_list_1, domain_list_2 =  self.db.get_domain_contact([d.uniprot_domain_1.pfam_name], [d.uniprot_domain_2.pfam_name])
            self._remove_bad_domain_pairs(domain_list_1)
            self._remove_bad_domain_pairs(domain_list_2)
            if not domain_list_1 and not domain_list_2:
                self.log.debug('Did not find any templates for domains {} and {}, trying their subdomains...'
                    .format(d.uniprot_domain_1.pfam_name, d.uniprot_domain_2.pfam_name))
                pfam_names_1 = self._split_superdomains(d.uniprot_domain_1.pfam_name)
                pfam_names_2 = self._split_superdomains(d.uniprot_domain_2.pfam_name)
                self.log.debug('Pfam domains {} and {} have subdomains {} and {}'
                    .format(d.uniprot_domain_1.pfam_name, d.uniprot_domain_2.pfam_name, pfam_names_1, pfam_names_2))
                domain_list_1, domain_list_2 =  self.db.get_domain_contact(pfam_names_1, pfam_names_2, subdomains=True)
                self._remove_bad_domain_pairs(domain_list_1)
                self._remove_bad_domain_pairs(domain_list_2)
            if not domain_list_1 and not domain_list_2:
                raise errors.NoStructuralTemplates('No templates found')

            # In the second list, domains are in opposite order relative to the query uniprots
            for idx in range(len(domain_list_2)):
                domain_contact = domain_list_2[idx]
                domain_contact.cath_id_1, domain_contact.cath_id_2 = \
                    domain_contact.cath_id_2, domain_contact.cath_id_1
                domain_contact.atom_count_1, domain_contact.atom_count_2 = \
                    domain_contact.atom_count_2, domain_contact.atom_count_1
                domain_contact.contact_residues_1, domain_contact.contact_residues_2 = \
                    domain_contact.contact_residues_2, domain_contact.contact_residues_1
                domain_contact.domain_1, domain_contact.domain_2 = \
                    domain_contact.domain_2, domain_contact.domain_1
                domain_list_2[idx] = domain_contact
            domain_list = list(set(domain_list_1 + domain_list_2))
            domain_list.sort(key=lambda dc: (dc.domain_1.cdhit_cluster_idx, dc.domain_2.cdhit_cluster_idx,))
            # Get the length of the longest first and second domains
            max_domain_length_1 = 0
            max_domain_length_2 = 0
            max_contact_length_1 = 0
            max_contact_length_2 = 0
            for domain in domain_list:
                domain_length_1 = \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(sql_db.decode_domain(domain.domain_1.pdb_domain_def, return_string=True)[1]) - \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(sql_db.decode_domain(domain.domain_1.pdb_domain_def, return_string=True)[0]) + 1
                domain_length_2 = \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(sql_db.decode_domain(domain.domain_2.pdb_domain_def, return_string=True)[1]) - \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(sql_db.decode_domain(domain.domain_2.pdb_domain_def, return_string=True)[0]) + 1
                if domain_length_1 > max_domain_length_1:
                    max_domain_length_1 = domain_length_1
                if domain_length_2 > max_domain_length_2:
                    max_domain_length_2 = domain_length_2
                contact_length_1 = len(sql_db.decode_aa_list(domain.contact_residues_1))
                contact_length_2 = len(sql_db.decode_aa_list(domain.contact_residues_2))
                if contact_length_1 > max_contact_length_1:
                    max_contact_length_1 = contact_length_1
                if contact_length_2 > max_contact_length_2:
                    max_contact_length_2 = contact_length_2

        #######################################################################
        # Filter out some of the templates that cannot be better than the cannonical template
        if canonical_domain:
            domain_list_new = []
            if isinstance(d, sql_db.UniprotDomain):
                for domain in domain_list:
                    if domain.cdhit_cluster != canonical_domain.domain.cdhit_cluster:
                        continue
                    if domain.pdb_type != 'X-ray' and canonical_domain.domain.pdb_type == 'X-ray':
                        continue
                    if domain.pdb_resolution >= canonical_domain.domain.pdb_resolution \
                    and domain.cdhit_cluster_identity > 99.9:
                        continue
                    domain_list_new.append(domain)
                    self.log.debug('keeping domain with cluster id: %i, cluster idx: %i' % (domain.cdhit_cluster, domain.cdhit_cluster_idx))
                domain_list = domain_list_new
            if isinstance(d, sql_db.UniprotDomainPair):
                for domain in domain_list:
                    if domain.domain_1.cdhit_cluster != canonical_domain.domain_1.cdhit_cluster:
                        continue
                    if domain.domain_2.cdhit_cluster != canonical_domain.domain_2.cdhit_cluster:
                        continue
                    if domain.domain_1.pdb_resolution >= canonical_domain.domain_1.pdb_resolution \
                    and domain.domain_1.cdhit_cluster_identity > 99.9 \
                    and domain.domain_2.cdhit_cluster_identity > 99.9:
                        continue
                    domain_list_new.append(domain)
                domain_list_new.append(domain)

        #######################################################################
        # Make and expand alignments for each of the obtained structural domain templates
        list_of_templates = []
        pfam_cluster_ids_visited = dict()
        self.log.debug('Found {} structural templates.'.format(len(domain_list)))
        for domain in domain_list:
            self.log.debug('-' * 80)
            if isinstance(d, sql_db.UniprotDomain):
                # Check if we need to evaluate this template
                pfam_cluster_id = (domain.pfam_name, domain.cdhit_cluster,)
                self.log.debug('Pfam cluster id: {}'.format(pfam_cluster_id))
                self.log.debug(
                    'Is {}a canonical domain!'
                    .format('' if canonical_domain else 'not '))
                self.log.debug(
                    'Has {}been visited!'
                    .format('' if pfam_cluster_ids_visited.has_key(pfam_cluster_id) else 'not '))
                if (not canonical_domain and
                        pfam_cluster_ids_visited.has_key(pfam_cluster_id)):
                    continue

                # Get the parameters required to make alignments
                uniprot_sequence = self.db.get_uniprot_sequence(d.uniprot_id)
                pdb_domain_def = sql_db.decode_domain(domain.pdb_domain_def, return_string=True)
                try:
                    pdb = pdb_template.PDBTemplate(
                        self.pdb_path, domain.pdb_id, [domain.pdb_chain], [pdb_domain_def, ],
                        self.unique_temp_folder, self.unique_temp_folder, self.log)
                except errors.NoPDBFoundError as e:
                    self.log.error(str(type(e)) + ': ' + e.__str__())
                    self.log.error("Didn't find the pdb file? Check if it is correct. Skipping...")
                    continue
                pdb.extract()
                __, chain_sequence = pdb.get_chain_numbering(domain.pdb_chain, return_sequence=True, return_extended=True)
                chain_sequence = SeqRecord(seq=Seq(chain_sequence), id=domain.pdb_id+domain.pdb_chain)
                domain_def = sql_db.decode_domain(d.envelope_def)
                template = sql_db.UniprotDomainTemplate()
                template.domain = domain
                template.cath_id = domain.cath_id
                template.uniprot_domain_id = d.uniprot_domain_id
                score_align = lambda alignment: self.score_align(alignment, max_domain_length, None, None)

                if not len(chain_sequence.seq):
                    self.log.error('PB chain is empty!')
                    self.log.error('PDB chain 1: {}'.format(chain_sequence.seq))
                    self.log.debug('Skipping...')
                    continue

                try: # Do iterative alignments and catch errors
                    self.log.debug(
                        'Aligning: {}/{}*{}:{}{}'
                        .format(d.uniprot_id, domain.pfam_name, domain.pdb_domain_def, domain.pdb_id, domain.pdb_chain))
                    (template.domain_def,
                    template.alignment_id,
                    template.alignment_score,
                    template.alignment_identity,
                    template.alignment_filename) = \
                        self.calculate_alignment(
                            uniprot_sequence, domain_def, chain_sequence,
                            score_align, refine=True)
                    self.log.debug(
                        'Done aligning: {}:{}'
                        .format(d.uniprot_id, template.alignment_id))
                except (errors.LowIdentity,
                            errors.EmptyPDBSequenceError) as e:
                        self.log.error(e.__str__())
                        self.log.debug('Skipping...')
                        continue
                except self.possible_template_expansion_errors as e:
                    raise e
                    domain.domain_errors = str(type(e)) + ': ' + e.__str__()
                    self.log.error(domain.domain_errors)
                    if 'result' in dir(e):
                        self.log.error(e.result)
                    self.db.add_domain(domain)
                    continue
                if template.alignment_identity < self.alignment_identity_cutoff:
                    continue

                # Save alignment results
                list_of_templates.append(template)
                pfam_cluster_ids_visited[pfam_cluster_id] = (domain.cdhit_cluster_idx,)
                self.log.debug('cdhit cluster: %i, cdhit cluster idx: %i' % (domain.cdhit_cluster, domain.cdhit_cluster_idx) )


            elif isinstance(d, sql_db.UniprotDomainPair):
                # Check if we need to evaluate this template
                pfam_cluster_id = (
                    domain.domain_1.pfam_name, domain.domain_1.cdhit_cluster,
                    domain.domain_2.pfam_name, domain.domain_2.cdhit_cluster,)
                self.log.debug('Pfam cluster id: {}'.format(pfam_cluster_id))
                self.log.debug(
                    'Is {}a canonical domain!'
                    .format('' if canonical_domain else 'not '))
                self.log.debug(
                    'Has {}been visited!'
                    .format('' if pfam_cluster_ids_visited.has_key(pfam_cluster_id) else 'not '))
                self.log.debug(
                    'Pfam cluster visited previously had cdhit cluster idxs: {}, {}'
                    .format(*pfam_cluster_ids_visited.get(pfam_cluster_id, [None, None])))
                if (not canonical_domain and
                        pfam_cluster_ids_visited.has_key(pfam_cluster_id) and
                        pfam_cluster_ids_visited[pfam_cluster_id][0] <= domain.domain_1.cdhit_cluster_idx and
                        pfam_cluster_ids_visited[pfam_cluster_id][1] <= domain.domain_2.cdhit_cluster_idx):
                    self.log.debug('Already covered this cluster pair. Skipping...')
                    continue

                # Get the parameters required to make alignments
                uniprot_sequence_1 = self.db.get_uniprot_sequence(d.uniprot_domain_1.uniprot_id)
                uniprot_sequence_2 = self.db.get_uniprot_sequence(d.uniprot_domain_2.uniprot_id)
                pdb_domain_def_1 = sql_db.decode_domain(domain.domain_1.pdb_domain_def, return_string=True)
                pdb_domain_def_2 = sql_db.decode_domain(domain.domain_2.pdb_domain_def, return_string=True)
                try:
                    pdb = pdb_template.PDBTemplate(
                        self.pdb_path,  domain.domain_1.pdb_id,
                        [domain.domain_1.pdb_chain, domain.domain_2.pdb_chain],
                        [pdb_domain_def_1, pdb_domain_def_2],
                        self.unique_temp_folder, self.unique_temp_folder, self.log)
                except errors.NoPDBFoundError as e:
                    self.log.error(str(type(e)) + ': ' + e.__str__())
                    self.log.error("Didn't find the pdb file? Check if it is correct. Skipping...")
                    continue
                pdb.extract()
                chain_numbering_1, chain_sequence_1 = pdb.get_chain_numbering(domain.domain_1.pdb_chain, return_sequence=True, return_extended=True)
                chain_sequence_1 = SeqRecord(seq=Seq(chain_sequence_1), id=domain.domain_1.pdb_id+domain.domain_1.pdb_chain)
                chain_numbering_2, chain_sequence_2 = pdb.get_chain_numbering(domain.domain_2.pdb_chain, return_sequence=True, return_extended=True)
                chain_sequence_2 = SeqRecord(seq=Seq(chain_sequence_2), id=domain.domain_2.pdb_id+domain.domain_2.pdb_chain)
                domain_def_1 = sql_db.decode_domain(d.uniprot_domain_1.envelope_def)
                domain_def_2 = sql_db.decode_domain(d.uniprot_domain_2.envelope_def)
                template = sql_db.UniprotDomainPairTemplate()
                template.uniprot_domain_pair_id = d.uniprot_domain_pair_id
                template.cath_id_1, template.cath_id_2 = domain.cath_id_1, domain.cath_id_2
                template.domain_1, template.domain_2 = domain.domain_1, domain.domain_2
                contact_residue_idxs_1 = [
                    chain_numbering_1.index(resid) for resid
                    in sql_db.decode_aa_list(domain.contact_residues_1)
                    if resid in chain_numbering_1]
                contact_residue_idxs_2 = [
                    chain_numbering_2.index(resid) for resid
                    in sql_db.decode_aa_list(domain.contact_residues_2)
                    if resid in chain_numbering_2]
                score_align_1 = lambda alignment: self.score_align(alignment, max_domain_length_1, contact_residue_idxs_1, max_contact_length_1)
                score_align_2 = lambda alignment: self.score_align(alignment, max_domain_length_2, contact_residue_idxs_2, max_contact_length_2)

                if not len(chain_sequence_1.seq) or not len(chain_sequence_2.seq):
                    self.log.error('At least one of the pdb chains is empty!')
                    self.log.error('PDB chain 1: {}'.format(chain_sequence_1.seq))
                    self.log.error('PDB chain 2: {}'.format(chain_sequence_2.seq))
                    self.log.debug('Skipping...')
                    continue

                try: # Do iterative alignments and catch errors
                    self.log.debug(
                        'Aligning partner 1: {}/{}*{}:{}{}'.format(
                            d.uniprot_domain_1.uniprot_id, domain.domain_1.pfam_name,
                            domain.domain_1.pdb_domain_def,
                            domain.domain_1.pdb_id, domain.domain_1.pdb_chain))
                    (template.domain_def_1,
                    template.alignment_id_1,
                    template.alignment_score_1,
                    template.alignment_identity_1,
                    template.alignment_filename_1) = \
                        self.calculate_alignment(
                            uniprot_sequence_1, domain_def_1, chain_sequence_1,
                            score_align_1, refine=True)
                    self.log.debug(
                        'Done aligning partner 1: {}:{}'.format(
                            d.uniprot_domain_1.uniprot_id, template.alignment_id_1))
                    self.log.debug(
                        'Aligning partner 2: {}/{}*{}:{}{}'.format(
                            d.uniprot_domain_2.uniprot_id, domain.domain_2.pfam_name,
                            domain.domain_2.pdb_domain_def,
                            domain.domain_2.pdb_id, domain.domain_2.pdb_chain))
                    (template.domain_def_2,
                    template.alignment_id_2,
                    template.alignment_score_2,
                    template.alignment_identity_2,
                    template.alignment_filename_2) = \
                        self.calculate_alignment(
                            uniprot_sequence_2, domain_def_2, chain_sequence_2,
                            score_align_2, refine=True)
                    self.log.debug(
                        'Done aligning partner 2: {}:{}'.format(
                            d.uniprot_domain_2.uniprot_id, template.alignment_id_2))
                except (errors.LowIdentity,
                            errors.EmptyPDBSequenceError) as e:
                        self.log.error(e.__str__())
                        self.log.debug('Skipping...')
                        continue
                except self.possible_template_expansion_errors as e:
                    raise e
                    domain.domain_contact_errors = str(type(e)) + ': ' + e.__str__()
                    self.log.error(domain.domain_contact_errors)
                    if 'result' in dir(e):
                        self.log.error(e.result)
                    self.db.add_domain(domain)
                    continue
                if (template.alignment_identity_1 < self.alignment_identity_cutoff or
                        template.alignment_identity_2 < self.alignment_identity_cutoff):
                    self.log.debug('Skipping...\n\n')
                    continue

                # Save alignment results
                pfam_cluster_ids_visited[pfam_cluster_id] = (
                    domain.domain_1.cdhit_cluster_idx, domain.domain_2.cdhit_cluster_idx,)
                list_of_templates.append(template)
                self.log.debug('Adding alignmnet to list...\n\n')

        # Return templates for either a single domain or a domain pair
        return list_of_templates


    def _split_superdomains(self, superdomain):
        domains = [d.split('_')[0] for d in superdomain.split('+')]
        return domains


    def _remove_bad_domains(self, domain_list):
        d_idx = 0
        while d_idx < len(domain_list):
            if domain_list[d_idx].domain_errors != None:
                self.log.debug('Domain has errors: {}. Skipping...'.format(domain_list[d_idx].domain_errors))
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].pdb_id in self.bad_pdbs:
                self.log.debug('Domain pdb is a known bad pdb. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].cdhit_cluster == None:
                self.log.debug('No cdhit cluster information. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].cdhit_cluster == -1:
                self.log.debug('Bad cdhit cluster (-1). Skipping...')
                del domain_list[d_idx]
                continue
            d_idx += 1


    def _remove_bad_domain_pairs(self, domain_list):
        d_idx = 0
        while d_idx < len(domain_list):
            if domain_list[d_idx].domain_contact_errors != None:
                self.log.debug('Domain has errors: {}. Skipping...'.format(domain_list[d_idx].domain_contact_errors))
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.pdb_id in self.bad_pdbs:
                self.log.debug('Domain pdb is a known bad pdb. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.pdb_chain == domain_list[d_idx].domain_2.pdb_chain:
                # For now, we are just focusing on interactions between different chains
                # in the pdb. This may be changed in the future.
                self.log.debug('Interacting domains are on the same chain in the pdb. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.cdhit_cluster == None\
            or domain_list[d_idx].domain_2.cdhit_cluster == None:
                self.log.debug(domain_list[d_idx])
                self.log.debug(domain_list[d_idx].domain_1)
                self.log.debug(domain_list[d_idx].domain_2)
                self.log.debug(domain_list[d_idx].domain_1.cath_id)
                self.log.debug(domain_list[d_idx].domain_2.cath_id)
                self.log.debug('No cdhit cluster information. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.cdhit_cluster == -1 \
            or domain_list[d_idx].domain_2.cdhit_cluster == -1:
                self.log.debug('Bad cdhit cluster (-1). Skipping...')
                del domain_list[d_idx]
                continue
            d_idx += 1


    def calculate_alignment(self, uniprot_sequence, domain_def, chain_sequence, score_align, refine=True):
        """
        """
        # get the uniprot sequence and align it to the pdb sequence
        uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
        alignment, alignment_id = self.do_align(uniprot_sequence_domain, chain_sequence, 'quick')
        #----------------------------------------------------------------------
        # Expanding domain boundaries
        uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
        self.log.debug(alignment_id)
        self.log.debug(uniprot_alignment_sequence.seq)
        self.log.debug(pdb_alignment_sequence.seq)

        left_extend_length, right_extend_length = [
            int((random.random() * 2 + 1) * extend_length)
            for extend_length
            in self.count_overhangs_and_gaps(uniprot_alignment_sequence)]
        alignment_score, alignment_identity, interface_score, global_coverage, local_coverage = score_align(alignment)
        self.log.debug('Extend left by: %i' % left_extend_length)
        self.log.debug('Extend right by: %i' % right_extend_length)
        if alignment_identity < 0.25:
            raise errors.LowIdentity(
                'Initial alignment identity using {} as template is too low ({:.2f}).'
                .format(alignment_id, alignment_identity))
            return domain_def, alignment_id, alignment_score, alignment_identity, None

        loop_counter = 0
        while (refine and
                ((abs(left_extend_length) > 2 and domain_def[0] > 1) or
                (abs(right_extend_length) > 2 and domain_def[1] < len(uniprot_sequence))) and
                (loop_counter < 5)):
            loop_counter += 1
            #------------------------------------------------------------------
            # Expand left domain boundaries
            if (domain_def[0] - 1 - left_extend_length) >= 0:
                domain_def[0] -= left_extend_length
            else:
                domain_def[0] = 1
            # Expand right domain boundaries
            if (domain_def[1] + right_extend_length) <= len(uniprot_sequence):
                domain_def[1] += right_extend_length
            else:
                domain_def[1] = len(uniprot_sequence)
            #------------------------------------------------------------------
            # get the uniprot sequence and align it to the pdb sequence
            uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
            alignment, alignment_id = self.do_align(uniprot_sequence_domain, chain_sequence, 'expresso')
            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
            self.log.debug(uniprot_alignment_sequence.seq)
            self.log.debug(pdb_alignment_sequence.seq)

            left_extend_length, right_extend_length = [-1 * extend_length for extend_length in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
            self.log.debug('Extend left by: %i' % left_extend_length)
            self.log.debug('Extend right by: %i' % right_extend_length)

            domain_def[0] -= left_extend_length
            domain_def[1] += right_extend_length
            #------------------------------------------------------------------
            # get the uniprot sequence and align it to the pdb sequence
            uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
            alignment, alignment_id = self.do_align(uniprot_sequence_domain, chain_sequence, 'expresso')
            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
            self.log.debug(uniprot_alignment_sequence.seq)
            self.log.debug(pdb_alignment_sequence.seq)

            left_extend_length, right_extend_length = \
                [int((random.random() * 2 + 1) * extend_length) for extend_length in self.count_overhangs_and_gaps(uniprot_alignment_sequence)]
            self.log.debug('Extend left by: %i' % left_extend_length)
            self.log.debug('Extend right by: %i' % right_extend_length)

#        # Removed this part
#        #----------------------------------------------------------------------
#        self.log.debug('Removing final overhangs...')
#        left_extend_length, right_extend_length = [
#            -1 * extend_length
#            for extend_length
#            in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
#        alignment_score, alignment_identity, interface_score = score_align(alignment)
#        self.log.debug('Extend left by: %i' % left_extend_length)
#        self.log.debug('Extend right by: %i' % right_extend_length)
#
#        domain_def[0] -= left_extend_length
#        domain_def[1] += right_extend_length
#
#        loop_counter = 0
#        while ((abs(left_extend_length) > 0 or abs(right_extend_length) > 0)
#                and (domain_def[0] < domain_def[1])):
#            loop_counter += 1
#            #------------------------------------------------------------------
#            # get the uniprot sequence and align it to the pdb sequence
#            uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
#            alignment, alignment_id = self.do_align(uniprot_sequence_domain, chain_sequence)
#            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
#            self.log.debug(uniprot_alignment_sequence.seq)
#            self.log.debug(pdb_alignment_sequence.seq)
#
#            left_extend_length, right_extend_length = [
#                -1 * extend_length
#                for extend_length
#                in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
#            alignment_score, alignment_identity, interface_score = score_align(alignment)
#            self.log.debug('Extend left by: %i' % left_extend_length)
#            self.log.debug('Extend right by: %i' % right_extend_length)
#
#            domain_def[0] -= left_extend_length
#            domain_def[1] += right_extend_length
#        #----------------------------------------------------------------------

        # Added this part
        #----------------------------------------------------------------------
        alignment_score, alignment_identity, interface_score, global_coverage, local_coverage = score_align(alignment)
#        if (alignment_identity < 0.95) and (local_coverage < 0.95):
#            self.log.debug('Performing the final, structural alignment...')
#            alignment, alignment_id = self.do_align(uniprot_sequence_domain, chain_sequence, 'expresso')
#            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
#            self.log.debug(uniprot_alignment_sequence.seq)
#            self.log.debug(pdb_alignment_sequence.seq)


        def get_overhangs(seqrec_1, seqrec_2, do_reversed=False):
            top_overhang = 0
            bottom_overhang = 0
            if do_reversed:
                custom_iterator = reversed(zip(str(seqrec_1.seq), str(seqrec_2.seq)))
            else:
                custom_iterator = zip(str(seqrec_1.seq), str(seqrec_2.seq))
            for aa_1, aa_2 in custom_iterator:
                if (aa_1 != '-') and (aa_2 != '-'):
                    break
                elif (aa_1 == '-') and (aa_2 == '-'):
                    bottom_overhang += 1
                elif (aa_1 == '-') and (aa_2 != '-'):
                    bottom_overhang += 1
                elif (aa_1 != '-') and (aa_2 == '-'):
                    top_overhang += 1
                else:
                    raise Exception("Didn't take something into account!")
            return top_overhang, bottom_overhang


        # Remove overhangs in the alignment
        seqrec_1, seqrec_2 = alignment
        left_uniprot_overhang, left_pdb_overhang = get_overhangs(seqrec_1, seqrec_2, False)
        right_uniprot_overhang, right_pdb_overhang = get_overhangs(seqrec_1, seqrec_2, True)
        left_overhang = left_uniprot_overhang + left_pdb_overhang
        right_overhang = (
            -(right_uniprot_overhang + right_pdb_overhang)
            if (right_uniprot_overhang + right_pdb_overhang) > 0
            else None)
        seqrec_1.seq = seqrec_1.seq[left_overhang:right_overhang]
        seqrec_2.seq = seqrec_2.seq[left_overhang:right_overhang]
        domain_def[0] += left_uniprot_overhang
        domain_def[1] -= right_uniprot_overhang
        domain_def = sql_db.encode_domain(domain_def) # turn it into a string object to be saved in the database
        alignment_score, alignment_identity, interface_score, global_coverage, local_coverage = score_align(alignment)
        alignment_filename = alignment[0].id + '_' + alignment[1].id + '.aln'
        try:
            AlignIO.write(alignment, self.unique_temp_folder + 'tcoffee/' + alignment_filename, 'clustal')
        except IndexError as e:
            raise errors.EmptyPDBSequenceError(str(type(e)) + ': ' + e.__str__())
        return domain_def, alignment_id, alignment_score, alignment_identity, alignment_filename


    def count_overhangs_and_gaps(self, alignment_sequence):
        overhang_length = [0, 0]
        loner_started = [False, False]
        loner_length = [0, 0]
        gap_started = [False, False]
        gap_length = [0, 0]
        second_loner_started = [False, False]
        second_loner_length = [0, 0]
        second_gap_started = [False, False]
        second_gap_length = [0, 0]
        second_gap_ended = [False, False]

        for i in range(len(alignment_sequence)):

            # Stop after reaching the half-way point
            if i >= len(alignment_sequence)-1-i:
                break
            if second_gap_ended[0] and second_gap_ended[1]:
                break

            # One step forward, one step reverse
            for direction_flag in [0, 1]:
                if direction_flag:
                    i = len(alignment_sequence)-1-i

                if alignment_sequence[i] == '-' and not loner_started[direction_flag]:
                    overhang_length[direction_flag] += 1
                elif alignment_sequence[i] == '-' and loner_started[direction_flag] and not second_loner_started[direction_flag]:
                    gap_started[direction_flag] = True
                    gap_length[direction_flag] += 1
                elif alignment_sequence[i] == '-' and second_loner_started[direction_flag] and not second_gap_ended[direction_flag]:
                    second_gap_started[direction_flag] = True
                    second_gap_length[direction_flag] += 1

                elif alignment_sequence[i] != '-' and not gap_started[direction_flag]:
                    loner_started[direction_flag] = True
                    loner_length[direction_flag] += 1
                elif alignment_sequence[i] != '-' and gap_started[direction_flag] and not second_gap_started[direction_flag]:
                    second_loner_started[direction_flag] = True
                    second_loner_length[direction_flag] += 1
                elif alignment_sequence[i] != '-' and second_gap_started[direction_flag]:
                    second_gap_ended[direction_flag] = True

                elif second_gap_ended[direction_flag]:
                    continue

                else:
                    raise Exception('Didn\'t take into account all possibilities!')

        extend_length = [None,None]
        for direction_flag in [0, 1]:
            if second_gap_length[direction_flag] * 2 > (loner_length[direction_flag] + second_loner_length[direction_flag]):
                extend_length[direction_flag] = second_gap_length[direction_flag] + gap_length[direction_flag] + overhang_length[direction_flag]
            elif gap_length[direction_flag] * 2 > loner_length[direction_flag]:
                extend_length[direction_flag] = gap_length[direction_flag] + overhang_length[direction_flag]
            else:
                extend_length[direction_flag] = overhang_length[direction_flag]

#        self.log.debug('Left overhang:\t %i,\t gap:\t %i,\t gap^2:\t %i,\t total:\t %i' % (overhang_length[0], gap_length[0], second_gap_length[0], extend_length[0]))
#        self.log.debug('Right overhang:\t %i,\t gap:\t %i,\t gap^2:\t %i,\t total:\t %i' % (overhang_length[1], gap_length[1], second_gap_length[1], extend_length[1]))

        return (extend_length[0], extend_length[1])


    def do_align(self, uniprot_sequence, pdb_sequence, mode):
        """
        Align the sequences in the seqFile.fasta file and return the alignment
        and the percentage identity

        input
        seqIDs              type 'list'     ;look like ['P01112', '1FOEB']


        alignments[0]       type class 'Bio.Align.MultipleSeqAlignment'
        score               type 'float'
        """
        # write both sequences to one file
        with open(self.unique_temp_folder + 'seqfiles.fasta', 'w') as seqFiles:
            SeqIO.write([uniprot_sequence, pdb_sequence], seqFiles, 'fasta')

        seqIDs = [uniprot_sequence.id, pdb_sequence.id]

#        self.log.debug('Calling tcoffee with parameters:')
#        self.log.debug('global_temp_path: {}'.format(self.global_temp_path))
#        self.log.debug('unique_temp_path: {}'.format(self.unique_temp_folder))
        tcoffee = call_tcoffee.tcoffee_alignment(
            self.global_temp_path,
            self.unique_temp_folder,
            [self.unique_temp_folder + 'seqfiles.fasta', ],
            seqIDs,
            self.n_cores,
            self.pdb_path,
            mode,
            self.log)
        alignments = tcoffee.align()

        return alignments[0], pdb_sequence.id


    def score_align(self, alignment, max_domain_length, contact_residue_idxs, max_contact_length):
        """
        """
        multiplier = 10000
        seq_identity = self.get_identity(alignment)
        global_coverage, local_coverage = self.get_coverage(alignment, max_domain_length)

        # New way to discourage seq identity < 40%
        a = 0.95
        if seq_identity < 0.40:
            score = a * (seq_identity)**2/0.40 * (global_coverage) + (1.0 - a) * (global_coverage)
        else:
            score = a * (seq_identity) * (global_coverage) + (1.0 - a) * (global_coverage)
        score = int(score * multiplier)

        # Below are the scoring metrics that were used previously. I left
        # them here for comparison.
        # see http://www.nature.com/nmeth/journal/v10/n1/full/nmeth.2289.html#methods
        # getting the score like Aloy did, based on sequence identity and coverage
        a = 0.95
        score2 = a * seq_identity * global_coverage + (1.0 - a) * global_coverage
        score2 = int(score2 * multiplier)

        # Another posibility...
        score3 = seq_identity**2 * global_coverage
        score3 = int(score3 * multiplier)

        self.log.debug(
            'Identity: {:.3f}; coverage: {:.3f}; score: {:n}; score2: {:n}; score3: {:n}'
            .format(seq_identity, global_coverage, score, score2, score3))

        interface_identity = None
        interface_score = None
        if contact_residue_idxs:
            interface_identity = self.get_interacting_identity(alignment, contact_residue_idxs, max_contact_length)
            interface_score = interface_identity * multiplier
            score = (score + interface_score) / 2
            self.log.debug(
                'Interface identity: {:.3f}; interface score: {:n}; final score: {:n}'
                .format(interface_identity, interface_score, score))

        return score, seq_identity, interface_score, global_coverage, local_coverage


    def get_coverage(self, alignment, global_max_domain_length=None):
        """
        Returns the coverage of the alginment in %

        input
        alignment       type class 'Bio.Align.MultipleSeqAlignment'>
        """

        seqrec_1, seqrec_2 = alignment
        len_seq_1 = len([aa for aa in str(seqrec_1.seq) if aa != '-'])
        len_seq_2 = len([aa for aa in str(seqrec_2.seq) if aa != '-'])
        domain_length = min([len_seq_1, len_seq_2])

        max_domain_length = max([len_seq_1, len_seq_2])
        if global_max_domain_length is None:
            global_max_domain_length = max_domain_length

        self.log.debug(
            'Length of template: {}. Length of longest template: {}.'
            .format(domain_length, global_max_domain_length))

        global_coverage = float(domain_length) / global_max_domain_length
        local_coverage = float(domain_length) / max_domain_length

        return global_coverage, local_coverage


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
        seqrec_1, seqrec_2 = alignment
        len_seq_1 = len([aa for aa in str(seqrec_1.seq) if aa != '-'])
        len_seq_2 = len([aa for aa in str(seqrec_2.seq) if aa != '-'])
        max_domain_length = max([len_seq_1, len_seq_2])
        longer_seqrec, shorter_seqrec = (seqrec_1, seqrec_2) if len_seq_1 >= len_seq_2 else (seqrec_2, seqrec_1)

        num_identical = 0
        for aa_1, aa_2 in zip(str(longer_seqrec.seq), str(shorter_seqrec.seq)):
            if aa_1 != '-':
                if aa_1 == aa_2:
                    num_identical += 1

        return float(num_identical) / max_domain_length


    def get_interacting_identity(self, alignment, contact_residue_idxs, max_contact_length):
        uniprot_alignment_sequence = alignment[0].seq.tostring()
        pdb_alignment_sequence = alignment[1].seq.tostring()
        num_of_same = 0
        num_of_different = 0
#        self.log.debug('Uniprot domain length: %i, pdb domain length: %i'
#            % (len(uniprot_alignment_sequence), len(pdb_alignment_sequence),))
        idx = -1
        for uniprot_aa, pdb_aa in zip(uniprot_alignment_sequence, pdb_alignment_sequence):
            if pdb_aa == '-':
                continue
            else:
                idx += 1
            if idx in contact_residue_idxs:
#                self.log.debug('Contact idx: %i%s' % (idx, pdb_aa))
                if uniprot_aa == pdb_aa:
                    num_of_same += 1
                else:
                    num_of_different += 1
        assert num_of_same + num_of_different == len(contact_residue_idxs)
        return float(num_of_same) / max_contact_length


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


    ###########################################################################

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
        if isinstance(domain_template[0], sql_db.UniprotDomainTemplate):
            domain_template.sort(key=lambda k: k.alignment_score, reverse=True)
            max_score = domain_template[0].alignment_score
            self.log.debug('max alignment score: %f' % max_score)
            # Collect all templates with the highest alignment score
            best_domain_templates = []
            for t in domain_template:
                self.log.debug('alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score, t.domain.cdhit_cluster, t.domain.cdhit_cluster_idx))
                if t.alignment_score == max_score:
                    best_domain_templates.append(t)
            best_domain_templates.sort(key=lambda k: k.domain.pdb_resolution, reverse=False)
            for t in best_domain_templates:
                self.log.debug('Best alignments:')
                self.log.debug('alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score, t.domain.cdhit_cluster, t.domain.cdhit_cluster_idx))

        if isinstance(domain_template[0], sql_db.UniprotDomainPairTemplate):
            domain_template.sort(key=lambda k: k.alignment_score_1 + k.alignment_score_2, reverse=True)
            max_score = domain_template[0].alignment_score_1 + domain_template[0].alignment_score_2
            self.log.debug('max alignment score: %f' % max_score)
            # Collect all templates with the highest alignment score
            best_domain_templates = []
            for t in domain_template:
                self.log.debug('domain 1. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_1, t.domain_1.cdhit_cluster, t.domain_1.cdhit_cluster_idx))
                self.log.debug('domain 2. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_2, t.domain_2.cdhit_cluster, t.domain_2.cdhit_cluster_idx))
                if t.alignment_score_1 + t.alignment_score_2 == max_score:
                    best_domain_templates.append(t)
            best_domain_templates.sort(key=lambda k: k.domain_1.pdb_resolution, reverse=False)
            for t in best_domain_templates:
                self.log.debug('Best alignments:')
                self.log.debug('domain 1. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_1, t.domain_1.cdhit_cluster, t.domain_1.cdhit_cluster_idx))
                self.log.debug('domain 2. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_2, t.domain_2.cdhit_cluster, t.domain_2.cdhit_cluster_idx))
        return best_domain_templates[0]


    ###########################################################################

    def build_provean_supporting_set(self, d, t):
        current_path = os.getcwd()

        uniprot_sequence = self.db.get_uniprot_sequence(d.uniprot_id).seq.tostring()
        domain_def = sql_db.decode_domain(t.domain_def)
        uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
        first_aa = uniprot_sequence_domain[0]
        uniprot_sequence_domain = SeqRecord(
            seq=Seq(uniprot_sequence_domain), id=str(d.uniprot_domain_id),
            description=d.uniprot_id + '_' + t.domain_def)
        SeqIO.write(uniprot_sequence_domain, self.unique_temp_folder + 'sequence_conservation/sequence.fasta', 'fasta')
        provean_supset_filename = t.alignment_filename.replace('.aln', '.supset')

        # Get initial measurements of how much virtual memory and disk space is availible
        disk_space_availible = psutil.disk_usage(self.provean_temp_path).free / float(1024)**3
        self.log.debug('Disk space availible: {}'.format(disk_space_availible))
        if disk_space_availible < 5:
            raise errors.ProveanError('Not enough disk space (%i GB) to run provean' % disk_space_availible)
        memory_availible = psutil.virtual_memory().available / float(1024)**3

        # Run provean
        subprocess.check_call('echo {0}1{0} > {1}sequence_conservation/decoy.var'.format(first_aa, self.unique_temp_folder), shell=True)
        system_command = (
            './provean ' +
            ' -q ' + './sequence.fasta ' +
            ' -v ' + './decoy.var ' +
            ' -d ' + self.global_temp_path + 'blast/db/nr ' +
            ' --tmp_dir ' + self.provean_temp_path +
            ' --num_threads ' + '{}'.format(self.n_cores) +
            ' --psiblast ' + hf.get_which('psiblast') +
            ' --blastdbcmd ' + hf.get_which('blastdbcmd') +
            ' --cdhit ' + hf.get_which('cd-hit') +
            ' --save_supporting_set ' + self.tmpPath + d.path_to_data + provean_supset_filename)
        self.log.debug(system_command)
        child_process = hf.run_subprocess_locally(
            self.unique_temp_folder + 'sequence_conservation/',
            system_command)
        self.log.debug('Parent group id: {}'.format(os.getpgrp()))
        child_process_group_id = os.getpgid(child_process.pid)
        self.log.debug('Child group id: {}'.format(child_process_group_id))

        # Keep an eye on provean to make sure it doesn't do anything crazy
        while True:
            if child_process.poll() is not None:
                break
            disk_space_availible_now = psutil.disk_usage(self.provean_temp_path).free / float(1024)**3
            if disk_space_availible_now < 5: # less than 5 GB of free disk space left
                raise errors.ProveanResourceError(
                    'Ran out of disk space and provean had to be terminated ({} GB used)'
                    .format(disk_space_availible-disk_space_availible_now),
                    child_process_group_id)
            memory_availible_now = psutil.virtual_memory().available / float(1024)**3
            if memory_availible_now < 0.5:
                raise errors.ProveanResourceError(
                    'Ran out of RAM ({} GB left)'
                    .format(memory_availible - memory_availible_now),
                    child_process_group_id)
            time.sleep(60) # Wait for 1 minute before checking again

        # Collect the results and check for errors
        result, error_message = child_process.communicate()
        return_code = child_process.returncode
        if return_code != 0:
            self.log.error(error_message)
            raise errors.ProveanError(error_message)

        ### Results look something like this:
        #[23:28:34] clustering subject sequences...
        #[23:28:34] selecting clusters...
        #[23:28:34] 0 subject sequences in 0 clusters were selected for supporting sequences.
        #[23:28:34] use the query itself as a supporting sequence
        #[23:28:34] loading subject sequences from a FASTA file...
        #[23:28:34] scores were computed based on the query sequence itself.
        ## Number of clusters:	1
        ## Number of supporting sequences used:	1
        #[23:28:34] computing delta alignment scores...
        #[23:28:34] printing PROVEAN scores...
        ### PROVEAN scores ##
        ## VARIATION	SCORE
        #M1A	-6.000

        # Commented out for now because we just need to get the supporting sequences
        for line in result.split('\n'):
            if 'Number of supporting sequences used:' in line:
                provean_supset_length = int(line.split()[-1])

        os.chdir(current_path)
        return provean_supset_filename, provean_supset_length


###############################################################################
# Obsolete

    def get_pdb_sequence(self, pdb_code, chain_id, pdb_domain_def):
        """
        Return the pdb file sequence (ATOM)
        """
        domains = [pdb_domain_def, ]
        pdb = pdb_template.PDBTemplate(self.pdb_path, pdb_code, chain_id, domains, self.unique_temp_folder, self.unique_temp_folder, self.log)
        pdb.extract()
        chain_numbering_extended, chain_sequence = pdb.get_chain_numbering(chain_id, return_sequence=True, return_extended=True)
        chain_seqrecord = SeqRecord(seq=Seq(chain_sequence), id=pdb_code+chain_id)

        # It happend that the given domain boundaries are larger than the
        # chain found in the pdb file. In this case, set the boundaries to the
        # maximum of the pdb chain
        if pdb_domain_def[0] not in chain_numbering_extended:
            pdb_domain_def[0] = chain_numbering_extended[0]
        if pdb_domain_def[1] not in chain_numbering_extended:
            pdb_domain_def[1] = chain_numbering_extended[-1]

        # Translate the pdb_domain from pdb numbering to sequence numbering
        try:
            pdb_domain_def_index = chain_numbering_extended.index(pdb_domain_def[0])+1, chain_numbering_extended.index(pdb_domain_def[1])+1
        except ValueError:
            raise errors.PDBError('ValueError when mapping domain boundaries to sequence numbering: ' + pdb_code + '_' + chain_id)

        return chain_seqrecord, pdb_domain_def_index, chain_numbering_extended, pdb_domain_def


def split_superdomains_1(superdomain):
    domains = [domain.split('_') for domain in superdomain.split('+')]
    print domains
    all_superdomains = set()
    for d1_idx in range(len(domains)):
        for d2_idx in range(d1_idx+1, len(domains)+1):
            all_superdomains.add('+'.join(domains[d1_idx:d2_idx]))
    for domain in list(all_superdomains):
        subdomains = domain.split('_')
        for i in range(2, len(subdomains)+1):
            all_superdomains.add('_'.join(subdomains[0:i]))
    return all_superdomains


def split_superdomains_2(superdomain):
    domains = [d.split('_') for d in superdomain.split('+')]
    data = np.empty( (len(domains), len(domains)), dtype=list)
    for row_idx in range(len(domains)):
        for col_idx in range(len(domains) - row_idx):
            data_sublist = data[row_idx, col_idx] = []
            if col_idx == 0:
                for d_1 in domains[row_idx]:
                    data_sublist.append( data_sublist[-1] + '_' + d_1 if data_sublist else d_1 )
            else:
                d_1_list = data[row_idx, col_idx - 1]
                d_2_list = []
                for d_2 in domains[row_idx + col_idx]:
                    d_2_list.append( d_2_list[-1] + '_' + d_2 if d_2_list else d_2 )
                for d_1 in d_1_list:
                    for d_2 in d_2_list:
                        data_sublist.append(d_1 + '+' + d_2)
    return [d for ddd in data[:] for dd in ddd if dd for d in dd]



if __name__ == '__main__':
    db = sql_db.MyDatabase('', path_to_archive='/home/kimlab1/database_data/elaspic/')
    tmp_path = '/tmp/test_template/'
    unique = 'Consumer-1'
    pdb_path = '/home/kimlab1/database_data/pdb/data/data/structures/divided/pdb/'
    saveAlignments = '/tmp/test_template/Consumer-1/'
    subprocess.check_call('mkdir -p ' + tmp_path + unique, shell=True)
    subprocess.check_call('mkdir -p ' + tmp_path + unique + '/sequence_conservation/', shell=True)
    subprocess.check_call('mkdir -p ' + tmp_path + 'blast/', shell=True)
    subprocess.check_call('cd ' + tmp_path + 'blast/ && ln -sf /home/kimlab1/strokach/ncbi-blast-2.2.28+/db', shell=True)
    subprocess.check_call('cp ' + '/home/kimlab1/strokach/working/pipeline/bin/provean ' + tmp_path + unique + '/sequence_conservation/', shell=True)

    ###########################################################################
    import logging
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
#    handler = logging.FileHandler(tmp_path + 'templates.log', mode='w', delay=True)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    ###########################################################################

    get_template = GetTemplate(tmp_path, unique, pdb_path, db, logger)
    temp = Bio.Align.MultipleSeqAlignment([Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq('--AGGA-')), Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq('MMMGGMM'))])
    print get_template.get_coverage(temp, len(temp[0]))
    print get_template.get_identity(temp)
    print get_template.get_interacting_identity(temp, [4,5], 2)

    p = db.get_uniprot_domain('Q8NEU8') + db.get_uniprot_domain_pair('Q8NEU8')
    p = [db.get_uniprot_domain_pair('Q8NEU8')[0]]
#    protein_domains = db.get_uniprot_domain_pair('Q8NEU8')
    for d, t, m in p:
        template = get_template(d)
#        if isinstance(template, sql_db.UniprotDomainTemplate):
#            print get_template.build_provean_supporting_set(d, template)

