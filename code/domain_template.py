# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:25:45 2013

@author: niklas
"""
import os
import psutil
import signal
import time
import subprocess
from math import fabs

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from pdb_template import pdbTemplate
import call_tcoffee
import errors
import sql_db
import pdb_template

import random

#from Bio.Align import MultipleSeqAlignment


def get_subdomains(superdomain):
    subdomains = superdomain.split('+')
    all_superdomains = set(subdomains)
    for i in range(2, len(subdomains)):
        for j in range(0, len(subdomains) + 1 - i):
            all_superdomains.add('+'.join(subdomains[j:i+j]))
    return all_superdomains


class GetTemplate():
    """
    Parent class holding functions for finding the correct template given a
    uniprot sequence
    """
    def __init__(self, global_temp_path, tmpPath, unique, pdbPath, db, log, n_cores, refine=False):
        """
        input:
        tmpPath             type 'str'
        unique              type 'str'
        pdbPath             type 'str'
        savePDB             type 'str'
        saveAlignments      type 'str'
        """
        self.global_temp_path = global_temp_path
        self.tmpPath = tmpPath
        self.unique = unique + '/'
        self.saveAlignments = self.tmpPath + self.unique + 'tcoffee/'
        self.pdbPath = pdbPath
        self.db = db

        # get the logger from the parent and add a handler
        self.log = log
        self.refine = refine
        self.bad_pdbs = ['3C4D', '3LB7', '3NSV', '2NP8', '2WN0']
        self.n_cores = n_cores


    def __call__(self, d):
        """
        """
        list_of_templates = self.run(d, None)
        if not list_of_templates:
            raise errors.NoTemplatesFound('Templates present in PDBfam were not useable')
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

        #######################################################################
        # Set up the paths and exporting the alignments
        if type(d) == sql_db.UniprotDomain:

            # Save one copy of the allignment for immediate usage in the tmp folder
            tmp_save_path = self.tmpPath + d.path_to_data
            subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
            subprocess.check_call('cp ' + self.saveAlignments + best_template.alignment_filename +
                                    ' ' + tmp_save_path + best_template.alignment_filename, shell=True)
#            (best_template.provean_supset_filename,
#            best_template.provean_supset_length) = self.build_provean_supporting_set(d, best_template)

        elif type(d) == sql_db.UniprotDomainPair:

            # Save one copy of the allignment for immediate usage in the tmp folder
            tmp_save_path = self.tmpPath + d.path_to_data
            subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
            subprocess.check_call('cp ' + self.saveAlignments + best_template.alignment_filename_1 +
                                    ' ' + tmp_save_path + best_template.alignment_filename_1, shell=True)
            subprocess.check_call('cp ' + self.saveAlignments + best_template.alignment_filename_2 +
                                    ' ' + tmp_save_path + best_template.alignment_filename_2, shell=True)

        #######################################################################

        return best_template


    def build_provean_supporting_set(self, d, t):

        current_path = os.getcwd()
        os.chdir(self.tmpPath + self.unique + 'sequence_conservation/')

        uniprot_sequence = self.db.get_uniprot_sequence(d.uniprot_id).seq.tostring()
        domain_def = sql_db.decode_domain(t.domain_def)
        uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
        first_aa = uniprot_sequence_domain[0]
        uniprot_sequence_domain = SeqRecord(seq=Seq(uniprot_sequence_domain), id=str(d.uniprot_domain_id), description=d.uniprot_id + '_' + t.domain_def)
        SeqIO.write(uniprot_sequence_domain, self.tmpPath + self.unique + 'sequence_conservation/sequence.fasta', 'fasta')

        disk_space_availible = psutil.disk_usage('/tmp').free / float(1024)**3
        if disk_space_availible < 5:
            raise errors.ProveanError('Not enough disk space (%i GB) to run provean' % disk_space_availible)
        memory_availible = psutil.virtual_memory().available / float(1024)**3
        provean_supset_filename = t.alignment_filename.replace('.aln', '.supset')
        system_command = (
            'echo ' + first_aa + '1' + first_aa + ' > ' + self.tmpPath + self.unique + 'sequence_conservation/decoy.var && ' +
            './provean ' +
            ' -q ' + self.tmpPath + self.unique + 'sequence_conservation/sequence.fasta ' +
            ' -v ' + self.tmpPath + self.unique + 'sequence_conservation/decoy.var ' +
            ' -d ' + self.global_temp_path + 'blast/db/nr ' +
            ' --tmp_dir ' + self.tmpPath + self.unique +
            ' --num_threads ' + '{}'.format(self.n_cores) +
            ' --psiblast `which psiblast` ' +
            ' --blastdbcmd `which blastdbcmd` ' +
            ' --cdhit `which cd-hit` ' +
            ' --save_supporting_set ' + self.tmpPath + d.path_to_data + provean_supset_filename)
        child_process = subprocess.Popen(system_command,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, preexec_fn=os.setsid)

        while True:
            if child_process.poll() is not None:
                break
            disk_space_availible_now = psutil.disk_usage('/tmp').free / float(1024)**3
            if disk_space_availible_now < 5: # less than 5 GB of free disk space left
                os.killpg(child_process.pid, signal.SIGTERM) # Send the signal to all the process groups
                raise errors.ProveanError(
                    'Ran out of disk space and provean had to be terminated ({} GB used)'
                    .format(disk_space_availible-disk_space_availible_now) )
                break
            memory_availible_now = psutil.virtual_memory().available / float(1024)**3
            if memory_availible_now < 1.5:
                os.killpg(child_process.pid, signal.SIGTERM) # Send the signal to all the process groups
                raise errors.ProveanError(
                    'Ran out of RAM ({} GB left)'
                    .format(memory_availible - memory_availible_now) )
                break
            time.sleep(60) # Wait for 1 minute before checking again

        result, error = child_process.communicate()
        if child_process.returncode != 0:
            self.log.error(error)
            raise errors.ProveanError(error)

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


    def run(self, d, canonical_domain, blacklisted_domain_ids=[]):
        """
        """
        if isinstance(d, sql_db.UniprotDomain):

            # Add all partial combinations of superdomains to the set of superdomains
            pfam_names = get_subdomains(d.pfam_name)
            global_cdhit_clusters = set()
            domain_list_all = []
            for pfam_name in pfam_names:
                local_cdhit_clusters = set()
                local_2_global = {}
                domain_list = self.db.get_domain(pfam_name)
                if not domain_list:
                    continue
                for domain in domain_list:
                    if domain.cath_id in blacklisted_domain_ids:
                        continue
                    elif domain.cdhit_cluster not in global_cdhit_clusters:
                        local_cdhit_clusters.add(domain.cdhit_cluster)
                    elif local_2_global.has_key(domain.cdhit_cluster):
                        domain.cdhit_cluster = local_2_global[domain.cdhit_cluster]
                    else:
                        for global_id in range(1000,2000):
                            if global_id not in (global_cdhit_clusters | local_cdhit_clusters):
                                break
                        local_2_global[domain.cdhit_cluster] = global_id
                        domain.cdhit_cluster = global_id
                        local_cdhit_clusters.add(domain.cdhit_cluster)
                global_cdhit_clusters.update(local_cdhit_clusters)
                domain_list_all.extend(domain_list)
            domain_list = domain_list_all

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
            pfam_names_1 = get_subdomains(d.uniprot_domain_1.pfam_name)
            pfam_names_2 = get_subdomains(d.uniprot_domain_2.pfam_name)
            domain_list_1_all = []
            domain_list_2_all = []
            for pfam_name_1 in pfam_names_1:
                for pfam_name_2 in pfam_names_2:
                    domain_list_1, domain_list_2 = self.db.get_domain_contact(pfam_name_1, pfam_name_2)
                    domain_list_1_all.extend(domain_list_1)
                    domain_list_2_all.extend(domain_list_2)
            domain_list_1, domain_list_2 = domain_list_1_all, domain_list_2_all
            # In the second list, domains are in opposite order relative to the query uniprots
            for idx, domain_contact in enumerate(domain_list_2):
                domain_contact.cath_id_1, domain_contact.cath_id_2 = \
                    domain_contact.cath_id_2, domain_contact.cath_id_1
                domain_contact.contact_residues_1, domain_contact.contact_residues_2 = \
                    domain_contact.contact_residues_2, domain_contact.contact_residues_1
                domain_contact.domain_1, domain_contact.domain_2 = \
                    domain_contact.domain_2, domain_contact.domain_1
                domain_list_2[idx] = domain_contact

            domain_list = list(set(domain_list_1 + domain_list_2))
            domain_list.sort(key=lambda dc: (dc.domain_1.cdhit_cluster_idx, dc.domain_2.cdhit_cluster_idx,))
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

        if len(domain_list) == 0:
            raise errors.NoStructuralTemplates('No templates found')

        #######################################################################
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
        cdhit_clusters_visited = dict()
        list_of_templates = []
        for domain in domain_list:
            if isinstance(d, sql_db.UniprotDomain):
                # there are some obsolete pdbs, ignore them... or 2NP8 has only one chain
                if domain.pdb_id in self.bad_pdbs:
                    continue
                if not domain.cdhit_cluster:
                    continue
                if domain.cdhit_cluster == -1:
                    continue
                if not canonical_domain:
                    if cdhit_clusters_visited.has_key(domain.cdhit_cluster):
                        continue
                template = sql_db.UniprotDomainTemplate()
                template.domain = domain
                template.cath_id = domain.cath_id
                template.uniprot_domain_id = d.uniprot_domain_id
                try:
                    (template.domain_def, template.alignment_id,
                    template.alignment_score, template.alignment_identity, template.alignment_filename) = \
                    self.calculate_alignment(d, domain, max_domain_length, None, None, refine=True, second_domain=False)
                except (
                errors.PDBDomainDefsError,
                errors.PDBChainError,
                errors.NoPDBFound, #PDBNotFoundError
                errors.EmptyPDBSequenceError, #PDBEmptySequenceError
                errors.pdbError, #PDBError
                errors.TcoffeeBlastError,
                errors.TcoffeePDBidError) as e:
                    self.log.error(str(type(e)))
                    self.log.error(e.message)
                    continue
                list_of_templates.append(template)
                cdhit_clusters_visited[domain.cdhit_cluster] = domain.cdhit_cluster_idx
                self.log.debug('cdhit cluster: %i, cdhit cluster idx: %i' % (domain.cdhit_cluster, domain.cdhit_cluster_idx) )

            elif isinstance(d, sql_db.UniprotDomainPair):
                # there are some obsolete pdbs, ignore them... or 2NP8 has only one chain
                if domain.domain_1.pdb_id in self.bad_pdbs:
                    continue
                if domain.domain_1.pdb_chain == domain.domain_2.pdb_chain:
                    # For now, we are just focusing on interactions between different chains
                    # in the pdb. This may be changed in the future.
                    continue
                if not domain.domain_1.cdhit_cluster \
                or not domain.domain_2.cdhit_cluster:
                    continue
                if domain.domain_1.cdhit_cluster == -1 \
                or domain.domain_2.cdhit_cluster == -1:
                    continue

                cdhit_key = (domain.domain_1.cdhit_cluster, domain.domain_2.cdhit_cluster,)
                if not canonical_domain \
                and cdhit_clusters_visited.has_key(cdhit_key) \
                and cdhit_clusters_visited[cdhit_key][0] <= domain.domain_1.cdhit_cluster_idx \
                and cdhit_clusters_visited[cdhit_key][1] <= domain.domain_2.cdhit_cluster_idx:
                    continue
                self.log.debug(cdhit_key)
                self.log.debug(not canonical_domain)
                self.log.debug(cdhit_clusters_visited.has_key(cdhit_key))
                self.log.debug(cdhit_clusters_visited.get(cdhit_key))

                template = sql_db.UniprotDomainPairTemplate()
                template.uniprot_domain_pair_id = d.uniprot_domain_pair_id
                template.cath_id_1, template.cath_id_2 = domain.cath_id_1, domain.cath_id_2
                template.domain_1, template.domain_2 = domain.domain_1, domain.domain_2
                self.log.debug('partner 1. pfam name: %s, cdhit cluster: %i, cdhit cluster idx: %i' \
                    % (domain.domain_1.pfam_name, domain.domain_1.cdhit_cluster, domain.domain_1.cdhit_cluster_idx) )
                self.log.debug('partner 2. pfam name: %s, cdhit cluster: %i, cdhit cluster idx: %i' \
                    % (domain.domain_2.pfam_name, domain.domain_2.cdhit_cluster, domain.domain_2.cdhit_cluster_idx) )

                try:
                    contact_residue_ressids_1 = pdb_template.get_pdb(domain.domain_1.pdb_id, self.pdbPath)[0]
                    self.log.debug('contact residue resids 1: ')
                    self.log.debug(contact_residue_ressids_1)
                    contact_residue_idxs_1 = pdb_template.convert_resid_to_position(
                        contact_residue_ressids_1,
                        domain.domain_1.pdb_chain,
                        sql_db.decode_aa_list(domain.contact_residues_1),
                        sql_db.decode_domain(domain.domain_1.pdb_domain_def, return_string=True)[0],
                        sql_db.decode_domain(domain.domain_1.pdb_domain_def, return_string=True)[1])
                    self.log.debug('contact residue idxs 1: ')
                    self.log.debug(contact_residue_idxs_1)

                    contact_residue_ressids_2 = pdb_template.get_pdb(domain.domain_2.pdb_id, self.pdbPath)[0]
                    self.log.debug('contact residue resids 2: ')
                    self.log.debug(contact_residue_ressids_2)
                    contact_residue_idxs_2 = pdb_template.convert_resid_to_position(
                        contact_residue_ressids_2,
                        domain.domain_2.pdb_chain,
                        sql_db.decode_aa_list(domain.contact_residues_2),
                        sql_db.decode_domain(domain.domain_2.pdb_domain_def, return_string=True)[0],
                        sql_db.decode_domain(domain.domain_2.pdb_domain_def, return_string=True)[1])
                    self.log.debug('contact residue idxs 2: ')
                    self.log.debug(contact_residue_idxs_2)

                    template.domain_def_1, template.alignment_id_1, \
                    template.alignment_score_1, template.alignment_identity_1, template.alignment_filename_1 = \
                        self.calculate_alignment(d, domain, max_domain_length_1, contact_residue_idxs_1,
                                                 max_contact_length_1, refine=True, second_domain=False)
                    template.domain_def_2, template.alignment_id_2, \
                    template.alignment_score_2, template.alignment_identity_2, template.alignment_filename_2 = \
                        self.calculate_alignment(d, domain, max_domain_length_2, contact_residue_idxs_2,
                                                 max_contact_length_2, refine=True, second_domain=True)
                except (
                errors.PDBDomainDefsError,
                errors.PDBChainError,
                errors.NoPDBFound, #PDBNotFoundError
                errors.EmptyPDBSequenceError, #PDBEmptySequenceError
                errors.pdbError, #PDBError
                errors.TcoffeeBlastError,
                errors.TcoffeePDBidError) as e:
                    self.log.error(str(type(e)))
                    self.log.error(e.message)
                    continue
                list_of_templates.append(template)
                cdhit_clusters_visited[(domain.domain_1.cdhit_cluster, domain.domain_2.cdhit_cluster)] = \
                    (domain.domain_1.cdhit_cluster_idx, domain.domain_2.cdhit_cluster_idx,)
                self.log.debug('cdhit cluster 1: %i, cdhit cluster 2: %i, cdhit cluster idx 1: %i, cdhit cluster idx 2: %i' % \
                (domain.domain_1.cdhit_cluster, domain.domain_2.cdhit_cluster, domain.domain_1.cdhit_cluster_idx, domain.domain_2.cdhit_cluster_idx,))

            return list_of_templates


    def calculate_alignment(self, d, domain, max_domain_length, contact_residue_idxs, max_contact_length, refine=True, second_domain=False):

        if isinstance(d, sql_db.UniprotDomainPair):
            if not second_domain:
                d = d.uniprot_domain_1
                domain = domain.domain_1
            else:
                d = d.uniprot_domain_2
                domain = domain.domain_2

        uniprot_id = d.uniprot_id

        domain_def = sql_db.decode_domain(d.envelope_def)
        pdb_id = domain.pdb_id
        pdb_chain = domain.pdb_chain
        pdb_domain_def = sql_db.decode_domain(domain.pdb_domain_def, return_string=True)

        uniprot_sequence = self.db.get_uniprot_sequence(uniprot_id)

        # get the sequence from the pdb
        pdb_sequence, pdb_domain_def, chainNumberingDomain, pdb_domain_def = self.get_pdb_sequence(pdb_id, pdb_chain, pdb_domain_def)
        domain.pdb_domain_def = sql_db.encode_domain(pdb_domain_def)

        self.log.debug("Aligning: " + uniprot_id + ':' + pdb_id + pdb_chain)

        # get the uniprot sequence and align it to the pdb sequence
        uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]

        alignment, alignment_id = self.map_to_uniprot_helper(uniprot_sequence_domain, pdb_sequence, self.saveAlignments)

        #----------------------------------------------------------------------
        # Expanding domain boundaries
        uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
        self.log.debug(alignment_id)
        self.log.debug(alignment)
        self.log.debug(uniprot_alignment_sequence.seq)
        self.log.debug(pdb_alignment_sequence.seq)

        left_extend_length, right_extend_length = \
            [int((random.random() * 2 + 1) * extend_length) for extend_length in self.count_overhangs_and_gaps(uniprot_alignment_sequence)]
        alignment_score, alignment_identity, interface_score = self.score_align(alignment, max_domain_length, contact_residue_idxs, max_contact_length)
        self.log.debug('Extend left by: %i' % left_extend_length)
        self.log.debug('Extend right by: %i' % right_extend_length)

        loop_counter = 0
        while ((abs(left_extend_length) > 0 and domain_def[0] > 1) \
            or (abs(right_extend_length) > 0 and domain_def[1] < len(uniprot_sequence))) \
        and loop_counter <= 6 \
        and refine:
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
            alignment, alignment_id = self.map_to_uniprot_helper(uniprot_sequence_domain, pdb_sequence, self.saveAlignments)
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
            alignment, alignment_id = self.map_to_uniprot_helper(uniprot_sequence_domain, pdb_sequence, self.saveAlignments)
            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
            self.log.debug(uniprot_alignment_sequence.seq)
            self.log.debug(pdb_alignment_sequence.seq)

            left_extend_length, right_extend_length = \
                [int((random.random() * 2 + 1) * extend_length) for extend_length in self.count_overhangs_and_gaps(uniprot_alignment_sequence)]
            alignment_score, alignment_identity, interface_score = self.score_align(alignment, max_domain_length, contact_residue_idxs, max_contact_length)
            self.log.debug('Extend left by: %i' % left_extend_length)
            self.log.debug('Extend right by: %i' % right_extend_length)


        #----------------------------------------------------------------------
        self.log.debug('Removing final overhangs...')
        left_extend_length, right_extend_length = [-1 * extend_length for extend_length in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
        alignment_score, alignment_identity, interface_score = self.score_align(alignment, max_domain_length, contact_residue_idxs, max_contact_length)
        self.log.debug('Extend left by: %i' % left_extend_length)
        self.log.debug('Extend right by: %i' % right_extend_length)

        domain_def[0] -= left_extend_length
        domain_def[1] += right_extend_length

        loop_counter = 0
        while (abs(left_extend_length) > 0 or abs(right_extend_length) > 0) \
        and (domain_def[0] < domain_def[1]):
            loop_counter += 1
            #------------------------------------------------------------------
            # get the uniprot sequence and align it to the pdb sequence
            uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
            alignment, alignment_id = self.map_to_uniprot_helper(uniprot_sequence_domain, pdb_sequence, self.saveAlignments)
            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
            self.log.debug(uniprot_alignment_sequence.seq)
            self.log.debug(pdb_alignment_sequence.seq)

            left_extend_length, right_extend_length = [-1 * extend_length for extend_length in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
            alignment_score, alignment_identity, interface_score = self.score_align(alignment, max_domain_length, contact_residue_idxs, max_contact_length)
            self.log.debug('Extend left by: %i' % left_extend_length)
            self.log.debug('Extend right by: %i' % right_extend_length)

            domain_def[0] -= left_extend_length
            domain_def[1] += right_extend_length

        self.log.debug("Done aligning: " + uniprot_id + ':' + alignment_id + '\n\n')

        domain_def = sql_db.encode_domain(domain_def) # turn it into a string object to be saved in the database
        alignment_filename = alignment[0].id + '_' + alignment[1].id + '.aln'

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
        with open(self.tmpPath + self.unique + 'seqfiles.fasta', 'w') as seqFiles:
            SeqIO.write([uniprot_sequence, pdb_sequence], seqFiles, 'fasta')

        seqIDs = [uniprot_sequence.id, pdb_sequence.id]

        # do the alignment and get the score
        alignment = self.do_align(seqIDs, saveAlignments)

        return alignment, pdb_sequence.id



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

        tcoffee = call_tcoffee.tcoffee_alignment(
            self.global_temp_path,
            self.tmpPath + self.unique,
            saveAlignments,
            [self.tmpPath + self.unique + 'seqfiles.fasta', ],
            seqIDs,
            self.n_cores)
        alignments = tcoffee.align()

        return alignments[0]


    def score_align(self, alignment, max_domain_length, contact_residue_idxs, max_contact_length):

        seq_identity = self.get_identity(alignment)/100.0  # percent -> decimal
        seq_coverage = self.get_coverage(alignment,max_domain_length)/100.0 # percent -> decimal

        interface_score = None
        if contact_residue_idxs:
            interface_score = self.get_interacting_identity(alignment, contact_residue_idxs, max_contact_length)/100.0 # percent -> decimal

        # New way to discourage seq identity < 40%
        a = 0.95
        if seq_identity < 0.40:
            score = a * (seq_identity)**2 / 0.40 * (seq_coverage) + (1.0 - a) * (seq_coverage)
        else:
            score = a * (seq_identity) * (seq_coverage) + (1.0 - a) * (seq_coverage)
        score = int(score*10000)

        # see http://www.nature.com/nmeth/journal/v10/n1/full/nmeth.2289.html#methods
        # getting the score like Aloy did, based on sequence identity and coverage
        a = 0.95
        score2 = a * seq_identity * seq_coverage + (1.0-a)*seq_coverage
        score2 = int(score2*10000)

        # another posibility...
        score3 = seq_identity**2*seq_coverage
        score3 = int(score3*10000)

        self.log.debug('Identity: %.3f, coverage: %.3f, score: %i, score2: %i, score3: %i, interface score: %s' \
            % (seq_identity, seq_coverage, score, score2, score3, interface_score))

        return score, seq_identity, interface_score


    def get_coverage(self, alignment, max_domain_length=None):
        """
        Returns the coverage of the alginment in %

        input
        alignment       type class 'Bio.Align.MultipleSeqAlignment'>
        """

        assert( len(alignment) == 2 )

        length_seq1 = len(alignment[0].seq.tostring().replace('-', ''))

        if not max_domain_length:
            # The old way of doing it, which assumes that the alignments are always
            # of the same length
            length_seq2 = len(alignment[1].seq.tostring().replace('-', ''))
        else:
            length_seq2 = max_domain_length

        self.log.debug('Length of target domain: %i; length of longest template: %i' % (length_seq1, length_seq2))

        if length_seq1 >= length_seq2:
            return 100.0 * float(length_seq2) / float(length_seq1)
        else:
            return 100.0 * float(length_seq1) / float(length_seq2)


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
        return 100.0 * float(num_of_same) / max_contact_length


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
                alignment, pdb_sequence.id = self.map_to_uniprot_helper(uniprot_sequence, pdb_sequence, saveAlignments)

        return alignment, [cut_left_all_uniprot, cut_right_all_uniprot], [cut_left_all_pdb, cut_right_all_pdb]


    def align_cut(self, alignment, pdb_sequence_id):
        """
        Checks for loners and shortens the sequence

        Parameters
        ----------
        alignment               type class 'Bio.Align.MultipleSeqAlignment'
        pdb_sequence_id         type 'str'

        Returns
        -------
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
        seq_left, uniprot_or_pdb_left = self.check_for_loners(uniprot_alignment, pdb_alignment, 'left')
        do_shorten_left, loners_left = self.check_loners(seq_left, 'left')
        # right
        seq_right, uniprot_or_pdb_right = self.check_for_loners(uniprot_alignment, pdb_alignment, 'right')
        do_shorten_right, loners_right = self.check_loners(seq_right, 'right')


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
                    self.log.error('You must specify uniprot or pdb!')
                    self.log.error('uniprot_or_pdb_right', uniprot_or_pdb_right)
                    self.log.error('uniprot_or_pdb_left', uniprot_or_pdb_left)
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
                    self.log.error('You must specify uniprot or pdb!')
                    self.log.error('uniprot_or_pdb_right', uniprot_or_pdb_right)
                    self.log.error('uniprot_or_pdb_left', uniprot_or_pdb_left)

        # since only the uniprot sequence is shortend, but above it works in
        # principle for both sequences simultaniously, a endless while loop can occur
        # to avoid this, check if the uniprot sequence is cut or not.
        if len(uniprot_sequence_old.seq) == len(uniprot_sequence.seq):
            cut = False

        return cut, cut_left, cut_right, uniprot_or_pdb_left, uniprot_or_pdb_right, uniprot_sequence, pdb_sequence


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
            self.log.error('You must specify left or right correctly!')
            return


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


    def get_pdb_sequence(self, pdbCode, chain, pdb_domain_def):
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
        domains = [pdb_domain_def, ]
        pdb = pdbTemplate(self.pdbPath, pdbCode, chain, domains, self.tmpPath + self.unique)

        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()

        __, chainNumberingDomain = pdb.getChainNumberingNOHETATMS(chain, return_extended=True)
        if chainNumberingDomain == []:
            raise errors.pdbError('Could not get the pdb numbering for ' + pdbCode + '_' + chain)

        # it happend that the given domain boundaries are larger than the
        # chain found in the pdb file. In this case, set the boundaries to the
        # maximum of the pdb chain
#        if domain_pdb[0] < chainNumberingDomain[0]:
#            domain_pdb[0] = chainNumberingDomain[0]
#        if domain_pdb[1] > chainNumberingDomain[-1]:
#            domain_pdb[1] = chainNumberingDomain[-1]

        ### AS changed
        if pdb_domain_def[0] not in chainNumberingDomain:
            pdb_domain_def[0] = chainNumberingDomain[0]
        if pdb_domain_def[1] not in chainNumberingDomain:
            pdb_domain_def[1] = chainNumberingDomain[-1]


        # translate the pdb_domain from pdb numbering to sequence numbering
        try:
            domain_pdb = chainNumberingDomain.index(pdb_domain_def[0])+1, chainNumberingDomain.index(pdb_domain_def[1])+1
        except ValueError:
            raise errors.pdbError('ValueError when mapping domain boundaries to sequence numbering: ' + pdbCode + '_' + chain)

        # seq.txt is the sequence used for crystalization and not the seqres sequence
        pdb_sequence = next(SeqIO.parse(self.tmpPath + self.unique + pdbCode + chain + '.seq.txt', 'fasta'))
        if str(pdb_sequence.seq) == '':
            raise errors.EmptyPDBSequenceError(pdbCode, chain)
        return pdb_sequence, domain_pdb, chainNumberingDomain, pdb_domain_def


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
#            self.log.error("PDB code, chain:", pdbCode, chain)
#            self.log.error("Domains:", domains)
#            self.log.error("Dobain_pdb 1/2:", domain_pdb[0], domain_pdb[1])
#            self.log.error("ChainNumberingDomain:")
#            self.log.error(chainNumberingDomain)
#            raise error.pdbError('ValueError when mapping domain boundaries to sequence numbering: ' + pdbCode + '_' + chain)
#
#        pdb_sequence = next(SeqIO.parse(self.tmpPath + self.unique + pdbCode + chain + '.seq.txt', 'fasta'))
#        return pdb_sequence, domain_pdb, chainNumberingDomain





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
            self.log.error('You must specify left or right!')
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


###############################################################################

if __name__ == '__main__':
    import Bio
    import logging

    db = sql_db.MyDatabase('', path_to_archive='/home/kimlab1/database_data/elaspic/')

    tmp_path = '/tmp/test_template/'
    unique = 'Consumer-1'
    pdbPath = '/home/kimlab1/database_data/pdb/data/data/structures/divided/pdb/'
    saveAlignments = '/tmp/test_template/Consumer-1/'
    subprocess.check_call('mkdir -p ' + tmp_path + unique, shell=True)
    subprocess.check_call('mkdir -p ' + tmp_path + unique + '/sequence_conservation/', shell=True)
    subprocess.check_call('mkdir -p ' + tmp_path + 'blast/', shell=True)
    subprocess.check_call('cd ' + tmp_path + 'blast/ && ln -sf /home/kimlab1/strokach/ncbi-blast-2.2.28+/db', shell=True)
    subprocess.check_call('cp ' + '/home/kimlab1/strokach/working/pipeline/bin/provean ' + tmp_path + unique + '/sequence_conservation/', shell=True)
    ###########################################################################

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

#    handler = logging.FileHandler(tmp_path + 'templates.log', mode='w', delay=True)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)

    ###########################################################################


    get_template = GetTemplate(tmp_path, unique, pdbPath, db, logger)

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

