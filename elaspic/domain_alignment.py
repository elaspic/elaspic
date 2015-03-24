# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import zip
from builtins import range
from builtins import object
from past.utils import old_div

import os
import psutil
import time
import subprocess
import random

import six
import numpy as np

import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import sql_db
from . import errors
from . import call_tcoffee
from . import pdb_template
from . import helper_functions as hf



def convert_basestring_to_seqrecord(sequence, sequence_id='id'):
    if any([isinstance(sequence, string_type) for string_type in six.string_types]):
        seqrec = SeqRecord(Seq(sequence), id=str(sequence_id))
    elif isinstance(sequence, Seq):
        seqrec = SeqRecord(sequence, id=str(sequence_id))
    elif isinstance(sequence, SeqRecord):
        seqrec = sequence
    else:
        raise Exception("Wrong class type %s for ``sequence``" % str(type(sequence)))

    return seqrec



def check_provean_supporting_set(
        domain_mutation, sequence, configs, unique_temp_folder, provean_temp_path, logger,
        sequence_id='id', path_to_provean_supset=None, save_supporting_set=False, check_mem_usage=False):
    """

    Provean results look something like this::

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

    Parameters
    ----------
    domain_mutation : string
        Mutation in domain coordinates (i.e. relative to the start of the domain)

    Returns
    -------
    list
        [result, error_message, return_code] -- The output from running a provean system command.

    Raises
    ------
    errors.ProveanError
        Can raise this exception only if ``check_mem_usage`` is set to ``True``.
    """

    if check_mem_usage:
        # Get initial measurements of how much virtual memory and disk space is availible
        disk_space_availible = psutil.disk_usage(provean_temp_path).free / (1024**3)
        logger.debug('Disk space availible: {:.2f} GB'.format(disk_space_availible))
        if disk_space_availible < 5:
            raise errors.ProveanError('Not enough disk space ({:.2f} GB) to run provean'.format(disk_space_availible))
        memory_availible = old_div(psutil.virtual_memory().available, float(1024)**3)
        logger.debug('Memory availible: {:.2f} GB'.format(memory_availible))
        if memory_availible < 0.5:
            raise errors.ProveanError('Not enough memory ({:.2f} GB) to run provean'.format(memory_availible))

    # Run provean
    seqrec = convert_basestring_to_seqrecord(sequence, sequence_id)
    SeqIO.write(seqrec, unique_temp_folder + 'sequence_conservation/sequence.fasta', 'fasta')

    subprocess.check_call(
        'echo {0} > {1}sequence_conservation/decoy.var'
        .format(domain_mutation, unique_temp_folder), shell=True)

    system_command = (
        './provean ' +
        ' -q ./sequence.fasta ' +
        ' -v ./decoy.var ' +
        ' -d ' + configs['global_temp_path'] + 'blast/db/nr ' +
        ' --tmp_dir ' + provean_temp_path +
        ' --num_threads ' + '{}'.format(configs['n_cores']) +
        ' --psiblast ' + hf.get_which('psiblast') +
        ' --blastdbcmd ' + hf.get_which('blastdbcmd') +
        ' --cdhit ' + hf.get_which('cd-hit')
    )
    if save_supporting_set:
        system_command += ' --save_supporting_set ' + path_to_provean_supset
    elif path_to_provean_supset: # use supporting set
        system_command += ' --supporting_set ' + path_to_provean_supset

    logger.debug(system_command)
    child_process = hf.run_subprocess_locally(
        unique_temp_folder + 'sequence_conservation/',
        system_command)

    logger.debug('Parent group id: {}'.format(os.getpgrp()))
    child_process_group_id = os.getpgid(child_process.pid)
    logger.debug('Child group id: {}'.format(child_process_group_id))

    # Keep an eye on provean to make sure it doesn't do anything crazy
    while check_mem_usage and child_process.poll() is None:
        disk_space_availible_now = old_div(psutil.disk_usage(provean_temp_path).free, float(1024)**3)
        if disk_space_availible_now < 5: # less than 5 GB of free disk space left
            raise errors.ProveanResourceError(
                'Ran out of disk space and provean had to be terminated ({} GB used)'
                .format(disk_space_availible-disk_space_availible_now),
                child_process_group_id)
        memory_availible_now = old_div(psutil.virtual_memory().available, float(1024)**3)
        if memory_availible_now < 0.5:
            raise errors.ProveanResourceError(
                'Ran out of RAM and provean had to be terminated ({} GB left)'
                .format(memory_availible - memory_availible_now),
                child_process_group_id)
        time.sleep(60) # Wait for 1 minute before checking again

    # Collect the results and check for errors
    result, error_message = child_process.communicate()
    if six.PY3:
        result = str(result, encoding='utf-8')
        error_message = str(error_message, encoding='utf-8')
    return_code = child_process.returncode

    return result, error_message, return_code



def build_provean_supporting_set(
        uniprot_id, uniprot_name, uniprot_sequence, 
        configs, unique_temp_folder, provean_temp_path, logger, 
        supset_version=0):
    """
    """
    # Get the required parameters
    first_aa = uniprot_sequence[0]
    domain_mutation = '{0}1{0}'.format(first_aa)

    uniprot_seqrecord = SeqRecord(
        seq=Seq(uniprot_sequence), id=str(uniprot_id), description=uniprot_name)

    provean_supset_path = (
        configs['temp_path'] + sql_db.get_uniprot_base_path(
            {'uniprot_id': uniprot_id, 'uniprot_name': uniprot_name}))
    provean_supset_filename = uniprot_id + '_provean_supset_{}'.format(supset_version)

    # Run provean
    result, error_message, return_code = check_provean_supporting_set(
        domain_mutation, uniprot_seqrecord, 
        configs, unique_temp_folder, provean_temp_path, logger,
        uniprot_id, provean_supset_path + provean_supset_filename,
        save_supporting_set=True, check_mem_usage=True)

    if return_code != 0:
        logger.error(error_message)
        raise errors.ProveanError(error_message)

    for line in result.split('\n'):
        if 'Number of supporting sequences used:' in line:
            provean_supset_length = int(line.split()[-1])

    return provean_supset_filename, provean_supset_length




class DomainPairTemplate(object):
    """
    R CODE PROTOTYPE:

    .. code-block:: r

        aligment.domains <- function(x,chain,dir_db,dir_query,dir_blastp,definition,libraries) {
            BE = definition[,x,drop=FALSE]
            family = gsub(pattern="\\|([[:digit:]])*$",replacement="",x=colnames(BE))
            B = BE["Begin",]
            E = BE["End",]
            if(family %in% libraries) {
              query = paste(dir_query,chain,".fasta",sep = "")
              db = paste(dir_db,family,"/",family,sep = "")
              cmd = paste("blastp -db",db,"-query",query,"-outfmt 10 -evalue 0.001")
              aligment = system(cmd,intern=TRUE)
              # Start inner if
              if (length(aligment) == 0) {
                templete = unique(x,chain,B,E,type=2)
              } else {
                aligment = lapply(aligment,formating)
                aligment = do.call(rbind,aligment)
                B2 = as.vector(aligment[,"q_start"])
                E2 = as.vector(aligment[,"q_end"])
                overlapings = mapply(overlap.internal,B1 = B,B2 = B2,E1 = E,E2 = E2,
                                     MoreArgs=list(type="percentage_1"),
                                     SIMPLIFY=TRUE)
                t_score = mapply(FUN = T_score,
                                 identity = as.vector(aligment[,"identity"]),
                                 coverage = overlapings,
                                 SIMPLIFY = TRUE)
                type = 1
                aligment = cbind(aligment,overlapings,t_score,type)
                dir_blastp_domain = paste(dir_blastp,x,".csv",sep="")
                write.csv(aligment,file=dir_blastp_domain)
                ind_templete = which.max(aligment[,"t_score"]) # criterion see t_score function
                templete = aligment[ind_templete,,drop=FALSE]
              }
            # End inner if
            } else {
              templete = unique(x,chain,B,E,type=3)
            }
            return(templete)
          }
    """

    def __init__(
            self, global_temp_path, temp_path, unique, pdb_path, db, logger,
            n_cores, provean_temp_path, refine=False):
        """
        """
        self.global_temp_path = global_temp_path
        self.temp_path = temp_path
        self.unique_temp_folder = temp_path + unique + '/'
        self.pdb_path = pdb_path
        self.db = db
        self.logger = logger
        self.refine = refine
        self.n_cores = n_cores





class GetTemplate(object):
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
        self.temp_path = tmpPath
        self.unique = unique
        self.unique_temp_folder = tmpPath + unique + '/'
        self.pdb_path = pdb_path
        self.db = db
        self.logger = log
        self.refine = refine
        self.n_cores = n_cores
        self.provean_temp_path = provean_temp_path


    def __call__(self, d):
        """
        """
        list_of_templates = self.run(d, None)
        if not list_of_templates:
            raise errors.NoTemplatesFound('Templates present in PDBfam were not usable')
        self.logger.debug('Number of templates: %i' % len(list_of_templates))

        best_template = self.chose_best_template(list_of_templates)
        self.logger.debug('The best template: ')
        self.logger.debug(best_template)

        # Don't need refinement anymore, each template is refined already
        if self.refine:
            self.logger.debug('Choosing the best templates withing a cluster...')
            list_of_templates = self.run(d, best_template)
            best_template = self.chose_best_template(list_of_templates)
            self.logger.debug('The best template in the cluster is: ')
            self.logger.debug(best_template)

        # Set up the paths and exporting the alignments
        if type(d) == sql_db.UniprotDomain:
            # Save one copy of the allignment for immediate usage in the tmp folder
            tmp_save_path = self.temp_path + d.path_to_data
            subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
            subprocess.check_call('cp ' + self.unique_temp_folder + 'tcoffee/' + best_template.alignment_filename +
                                    ' ' + tmp_save_path + best_template.alignment_filename, shell=True)
#            (best_template.provean_supset_filename,
#            best_template.provean_supset_length) = self.build_provean_supporting_set(d, best_template)

        elif type(d) == sql_db.UniprotDomainPair:
            # Save one copy of the allignment for immediate usage in the tmp folder
            tmp_save_path = self.temp_path + d.path_to_data
            subprocess.check_call('mkdir -p ' + tmp_save_path, shell=True)
            subprocess.check_call('cp ' + self.unique_temp_folder + 'tcoffee/' + best_template.alignment_filename_1 +
                                    ' ' + tmp_save_path + best_template.alignment_filename_1, shell=True)
            subprocess.check_call('cp ' + self.unique_temp_folder + 'tcoffee/' + best_template.alignment_filename_2 +
                                    ' ' + tmp_save_path + best_template.alignment_filename_2, shell=True)
        return best_template


    def run(self, d, canonical_domain):
        """
        Parameters
        ----------
        d : sql_db.UniprotDomain or sql_db.UniprotDomainPair
            Information about the given domain
        domain : sql_db.Domain or sql_db.DomainContact
            Information about the particular domain domain interaction
        """
        #######################################################################
        # Obtain a list of template domains and calculate some of their properties
        if isinstance(d, sql_db.UniprotDomain):
            # Add all partial combinations of superdomains to the set of superdomains
            self.logger.debug('Finding templates for pfam domain {}...'.format(d.pfam_name))
            domain_list = self.db.get_domain([d.pfam_name])
            self._remove_bad_domains(domain_list)
            if not domain_list:
                self.logger.debug('Did not find any structural templates for {}, trying subdomains...'.format(d.pfam_name))
                pfam_names = self._split_superdomains(d.pfam_name)
                self.logger.debug('Pfam domain {} has subdomains {}'.format(d.pfam_name, str(pfam_names)))
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
                    pdb_template.convert_resnum_alphanumeric_to_numeric(
                        sql_db.decode_domain(domain.pdb_domain_def, return_string=True)[1]) - \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(
                        sql_db.decode_domain(domain.pdb_domain_def, return_string=True)[0]) + 1
                if domain_length > max_domain_length:
                    max_domain_length = domain_length

        elif isinstance(d, sql_db.UniprotDomainPair):
            self.logger.debug('Finding templates for pfam domains {}, {}...'
                .format(d.uniprot_domain_1.pfam_name, d.uniprot_domain_2.pfam_name))
            domain_list_1, domain_list_2 =  self.db.get_domain_contact(
                [d.uniprot_domain_1.pfam_name], [d.uniprot_domain_2.pfam_name])
            self._remove_bad_domain_pairs(domain_list_1)
            self._remove_bad_domain_pairs(domain_list_2)
            if not domain_list_1 and not domain_list_2:
                self.logger.debug('Did not find any templates for domains {} and {}, trying their subdomains...'
                    .format(d.uniprot_domain_1.pfam_name, d.uniprot_domain_2.pfam_name))
                pfam_names_1 = self._split_superdomains(d.uniprot_domain_1.pfam_name)
                pfam_names_2 = self._split_superdomains(d.uniprot_domain_2.pfam_name)
                self.logger.debug('Pfam domains {} and {} have subdomains {} and {}'
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
                    pdb_template.convert_resnum_alphanumeric_to_numeric(
                        sql_db.decode_domain(domain.domain_1.pdb_domain_def, return_string=True)[1]) - \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(
                        sql_db.decode_domain(domain.domain_1.pdb_domain_def, return_string=True)[0]) + 1
                domain_length_2 = \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(
                        sql_db.decode_domain(domain.domain_2.pdb_domain_def, return_string=True)[1]) - \
                    pdb_template.convert_resnum_alphanumeric_to_numeric(
                        sql_db.decode_domain(domain.domain_2.pdb_domain_def, return_string=True)[0]) + 1
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
                    self.logger.debug(
                        'keeping domain with cluster id: %i, cluster idx: %i' %
                        (domain.cdhit_cluster, domain.cdhit_cluster_idx))
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
        self.logger.debug('Found {} structural templates.'.format(len(domain_list)))
        for domain in domain_list:
            self.logger.debug('-' * 80)
            if isinstance(d, sql_db.UniprotDomain):
                # Check if we need to evaluate this template
                pfam_cluster_id = (domain.pfam_name, domain.cdhit_cluster,)
                self.logger.debug('Pfam cluster id: {}'.format(pfam_cluster_id))
                self.logger.debug(
                    'Is {}a canonical domain!'
                    .format('' if canonical_domain else 'not '))
                self.logger.debug(
                    'Has {}been visited!'
                    .format('' if pfam_cluster_id in pfam_cluster_ids_visited else 'not '))
                if (not canonical_domain and
                        pfam_cluster_id in pfam_cluster_ids_visited):
                    continue

                # Get the parameters required to make alignments
                uniprot_sequence = self.db.get_uniprot_sequence(d.uniprot_id)
                pdb_domain_def = sql_db.decode_domain(domain.pdb_domain_def, return_string=True)
                try:
                    pdb = pdb_template.PDBTemplate(
                        self.pdb_path, domain.pdb_id, [domain.pdb_chain], [pdb_domain_def, ],
                        self.unique_temp_folder, self.unique_temp_folder, self.logger)
                except errors.NoPDBFoundError as e:
                    self.logger.error(str(type(e)) + ': ' + e.__str__())
                    self.logger.error("Didn't find the pdb file? Check if it is correct. Skipping...")
                    continue
                pdb.extract()
                chain_sequence, __ = pdb.get_chain_sequence_and_numbering(domain.pdb_chain)
                chain_sequence = SeqRecord(seq=Seq(chain_sequence), id=domain.pdb_id+domain.pdb_chain)
                domain_def = sql_db.decode_domain(d.envelope_def)
                template = sql_db.UniprotDomainTemplate()
                template.domain = domain
                template.cath_id = domain.cath_id
                template.uniprot_domain_id = d.uniprot_domain_id
                score_align = lambda alignment: self.score_align(alignment, max_domain_length, None, None)

                if not len(chain_sequence.seq):
                    self.logger.error('PB chain is empty!')
                    self.logger.error('PDB chain 1: {}'.format(chain_sequence.seq))
                    self.logger.debug('Skipping...')
                    continue

                try: # Do iterative alignments and catch errors
                    self.logger.debug(
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
                    self.logger.debug(
                        'Done aligning: {}:{}'
                        .format(d.uniprot_id, template.alignment_id))
                except (errors.LowIdentity,
                            errors.EmptyPDBSequenceError) as e:
                        self.logger.error(e.__str__())
                        self.logger.debug('Skipping...\n\n')
                        continue
                except self.possible_template_expansion_errors as e:
                    raise e
                    domain.domain_errors = str(type(e)) + ': ' + e.__str__()
                    self.logger.error(domain.domain_errors)
                    if 'result' in dir(e):
                        self.logger.error(e.result)
                    self.db.add_domain(domain)
                    self.logger.debug('Skipping...\n\n')
                    continue
                if template.alignment_identity < self.alignment_identity_cutoff:
                    self.logger.debug('Skipping...\n\n')
                    continue

                # Save alignment results
                list_of_templates.append(template)
                pfam_cluster_ids_visited[pfam_cluster_id] = (domain.cdhit_cluster_idx,)
                list_of_templates.append(template)
                self.logger.debug('Adding alignmnet to list...\n\n')


            elif isinstance(d, sql_db.UniprotDomainPair):
                # Check if we need to evaluate this template
                pfam_cluster_id = (
                    domain.domain_1.pfam_name, domain.domain_1.cdhit_cluster,
                    domain.domain_2.pfam_name, domain.domain_2.cdhit_cluster,)
                self.logger.debug('Pfam cluster id: {}'.format(pfam_cluster_id))
                self.logger.debug(
                    'Is {}a canonical domain!'
                    .format('' if canonical_domain else 'not '))
                self.logger.debug(
                    'Has {}been visited!'
                    .format('' if pfam_cluster_id in pfam_cluster_ids_visited else 'not '))
                self.logger.debug(
                    'Pfam cluster visited previously had cdhit cluster idxs: {}, {}'
                    .format(*pfam_cluster_ids_visited.get(pfam_cluster_id, [None, None])))
                if (not canonical_domain and
                        pfam_cluster_id in pfam_cluster_ids_visited and
                        pfam_cluster_ids_visited[pfam_cluster_id][0] <= domain.domain_1.cdhit_cluster_idx and
                        pfam_cluster_ids_visited[pfam_cluster_id][1] <= domain.domain_2.cdhit_cluster_idx):
                    self.logger.debug('Already covered this cluster pair. Skipping...')
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
                        self.unique_temp_folder, self.unique_temp_folder, self.logger)
                except errors.NoPDBFoundError as e:
                    self.logger.error(str(type(e)) + ': ' + e.__str__())
                    self.logger.error("Didn't find the pdb file? Check if it is correct. Skipping...")
                    continue
                pdb.extract()
                chain_sequence_1, chain_numbering_1 = pdb.get_chain_sequence_and_numbering(domain.domain_1.pdb_chain)
                chain_sequence_1 = SeqRecord(seq=Seq(chain_sequence_1), id=domain.domain_1.pdb_id+domain.domain_1.pdb_chain)
                chain_sequence_2, chain_numbering_2 = pdb.get_chain_sequence_and_numbering(domain.domain_2.pdb_chain)
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
                score_align_1 = lambda alignment: self.score_align(
                    alignment, max_domain_length_1, contact_residue_idxs_1, max_contact_length_1)
                score_align_2 = lambda alignment: self.score_align(
                    alignment, max_domain_length_2, contact_residue_idxs_2, max_contact_length_2)

                if not len(chain_sequence_1.seq) or not len(chain_sequence_2.seq):
                    self.logger.error('At least one of the pdb chains is empty!')
                    self.logger.error('PDB chain 1: {}'.format(chain_sequence_1.seq))
                    self.logger.error('PDB chain 2: {}'.format(chain_sequence_2.seq))
                    self.logger.debug('Skipping...')
                    continue

                try: # Do iterative alignments and catch errors
                    self.logger.debug(
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
                    self.logger.debug(
                        'Done aligning partner 1: {}:{}'.format(
                            d.uniprot_domain_1.uniprot_id, template.alignment_id_1))
                    self.logger.debug(
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
                    self.logger.debug(
                        'Done aligning partner 2: {}:{}'.format(
                            d.uniprot_domain_2.uniprot_id, template.alignment_id_2))
                except (errors.LowIdentity,
                            errors.EmptyPDBSequenceError) as e:
                        self.logger.error(e.__str__())
                        self.logger.debug('Skipping...\n\n')
                        continue
                except self.possible_template_expansion_errors as e:
                    raise e
                    domain.domain_contact_errors = str(type(e)) + ': ' + e.__str__()
                    self.logger.error(domain.domain_contact_errors)
                    if 'result' in dir(e):
                        self.logger.error(e.result)
                    self.db.add_domain(domain)
                    continue
                if (template.alignment_identity_1 < self.alignment_identity_cutoff or
                        template.alignment_identity_2 < self.alignment_identity_cutoff):
                    self.logger.debug('Skipping...\n\n')
                    continue

                # Save alignment results
                pfam_cluster_ids_visited[pfam_cluster_id] = (
                    domain.domain_1.cdhit_cluster_idx, domain.domain_2.cdhit_cluster_idx,)
                list_of_templates.append(template)
                self.logger.debug('Adding alignmnet to list...\n\n')

        # Return templates for either a single domain or a domain pair
        return list_of_templates


    def _split_superdomains(self, superdomain):
        domains = [d.split('_')[0] for d in superdomain.split('+')]
        return domains


    def _remove_bad_domains(self, domain_list):
        d_idx = 0
        while d_idx < len(domain_list):
            if domain_list[d_idx].domain_errors != None:
                self.logger.debug('Domain has errors: {}. Skipping...'.format(domain_list[d_idx].domain_errors))
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].pdb_id in self.bad_pdbs:
                self.logger.debug('Domain pdb is a known bad pdb. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].cdhit_cluster == None:
                self.logger.debug('No cdhit cluster information. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].cdhit_cluster == -1:
                self.logger.debug('Bad cdhit cluster (-1). Skipping...')
                del domain_list[d_idx]
                continue
            d_idx += 1


    def _remove_bad_domain_pairs(self, domain_list):
        d_idx = 0
        while d_idx < len(domain_list):
            if domain_list[d_idx].domain_contact_errors != None:
                self.logger.debug('Domain has errors: {}. Skipping...'.format(domain_list[d_idx].domain_contact_errors))
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.pdb_id in self.bad_pdbs:
                self.logger.debug('Domain pdb is a known bad pdb. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.pdb_chain == domain_list[d_idx].domain_2.pdb_chain:
                # For now, we are just focusing on interactions between different chains
                # in the pdb. This may be changed in the future.
                self.logger.debug('Interacting domains are on the same chain in the pdb. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.cdhit_cluster == None\
            or domain_list[d_idx].domain_2.cdhit_cluster == None:
                self.logger.debug(domain_list[d_idx])
                self.logger.debug(domain_list[d_idx].domain_1)
                self.logger.debug(domain_list[d_idx].domain_2)
                self.logger.debug(domain_list[d_idx].domain_1.cath_id)
                self.logger.debug(domain_list[d_idx].domain_2.cath_id)
                self.logger.debug('No cdhit cluster information. Skipping...')
                del domain_list[d_idx]
                continue
            if domain_list[d_idx].domain_1.cdhit_cluster == -1 \
            or domain_list[d_idx].domain_2.cdhit_cluster == -1:
                self.logger.debug('Bad cdhit cluster (-1). Skipping...')
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
        self.logger.debug(alignment_id)
        self.logger.debug(uniprot_alignment_sequence.seq)
        self.logger.debug(pdb_alignment_sequence.seq)

        left_extend_length, right_extend_length = [
            int((random.random() * 2 + 1) * extend_length)
            for extend_length
            in self.count_overhangs_and_gaps(uniprot_alignment_sequence)]
        alignment_score, alignment_identity, interface_score, global_coverage, local_coverage = score_align(alignment)
        self.logger.debug('Extend left by: %i' % left_extend_length)
        self.logger.debug('Extend right by: %i' % right_extend_length)
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
            self.logger.debug(uniprot_alignment_sequence.seq)
            self.logger.debug(pdb_alignment_sequence.seq)

            left_extend_length, right_extend_length = [
                -1 * extend_length for extend_length in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
            self.logger.debug('Extend left by: %i' % left_extend_length)
            self.logger.debug('Extend right by: %i' % right_extend_length)

            domain_def[0] -= left_extend_length
            domain_def[1] += right_extend_length
            #------------------------------------------------------------------
            # get the uniprot sequence and align it to the pdb sequence
            uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
            alignment, alignment_id = self.do_align(uniprot_sequence_domain, chain_sequence, 'expresso')
            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
            self.logger.debug(uniprot_alignment_sequence.seq)
            self.logger.debug(pdb_alignment_sequence.seq)

            left_extend_length, right_extend_length = \
                [int((random.random() * 2 + 1) * extend_length)
                for extend_length in self.count_overhangs_and_gaps(uniprot_alignment_sequence)]
            self.logger.debug('Extend left by: %i' % left_extend_length)
            self.logger.debug('Extend right by: %i' % right_extend_length)

#        # Removed this part
#        #----------------------------------------------------------------------
#        self.logger.debug('Removing final overhangs...')
#        left_extend_length, right_extend_length = [
#            -1 * extend_length
#            for extend_length
#            in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
#        alignment_score, alignment_identity, interface_score = score_align(alignment)
#        self.logger.debug('Extend left by: %i' % left_extend_length)
#        self.logger.debug('Extend right by: %i' % right_extend_length)
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
#            self.logger.debug(uniprot_alignment_sequence.seq)
#            self.logger.debug(pdb_alignment_sequence.seq)
#
#            left_extend_length, right_extend_length = [
#                -1 * extend_length
#                for extend_length
#                in self.count_overhangs_and_gaps(pdb_alignment_sequence)]
#            alignment_score, alignment_identity, interface_score = score_align(alignment)
#            self.logger.debug('Extend left by: %i' % left_extend_length)
#            self.logger.debug('Extend right by: %i' % right_extend_length)
#
#            domain_def[0] -= left_extend_length
#            domain_def[1] += right_extend_length
#        #----------------------------------------------------------------------

        # Added this part
        #----------------------------------------------------------------------
        alignment_score, alignment_identity, interface_score, global_coverage, local_coverage = score_align(alignment)
#        if (alignment_identity < 0.95) and (local_coverage < 0.95):
#            self.logger.debug('Performing the final, structural alignment...')
#            alignment, alignment_id = self.do_align(uniprot_sequence_domain, chain_sequence, 'expresso')
#            uniprot_alignment_sequence, pdb_alignment_sequence = self.pick_sequence(alignment, alignment_id)
#            self.logger.debug(uniprot_alignment_sequence.seq)
#            self.logger.debug(pdb_alignment_sequence.seq)


        def get_overhangs(seqrec_1, seqrec_2, do_reversed=False):
            top_overhang = 0
            bottom_overhang = 0
            if do_reversed:
                custom_iterator = reversed(list(zip(str(seqrec_1.seq), str(seqrec_2.seq))))
            else:
                custom_iterator = list(zip(str(seqrec_1.seq), str(seqrec_2.seq)))
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
                elif (alignment_sequence[i] == '-' and
                    second_loner_started[direction_flag] and
                    not second_gap_ended[direction_flag]):
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
                extend_length[direction_flag] = (
                    second_gap_length[direction_flag] + gap_length[direction_flag] + overhang_length[direction_flag]
                )
            elif gap_length[direction_flag] * 2 > loner_length[direction_flag]:
                extend_length[direction_flag] = gap_length[direction_flag] + overhang_length[direction_flag]
            else:
                extend_length[direction_flag] = overhang_length[direction_flag]

#        self.logger.debug('Left overhang:\t %i,\t gap:\t %i,\t gap^2:\t %i,\t total:\t %i' % (
#            overhang_length[0], gap_length[0], second_gap_length[0], extend_length[0]))
#        self.logger.debug('Right overhang:\t %i,\t gap:\t %i,\t gap^2:\t %i,\t total:\t %i' % (
#            overhang_length[1], gap_length[1], second_gap_length[1], extend_length[1]))

        return (extend_length[0], extend_length[1])


    def do_align(self, uniprot_sequence, pdb_sequence, mode):
        """
        Align the sequences in the seqFile.fasta file and return the alignment
        and the percentage identity.

        Parameters
        ----------
        uniprot_sequence : SeqRecord
            The sequence of the protein domain to be aligned
        pdb_sequence : SeqRecord
            The sequence of the PDB domain to be aligned
        mode : str
            The type of alignment to perform.
            Can be "3dcoffee", "expresso", "t_coffee", "quick".
            
        Returns
        -------
        Bio.Align.MultipleSeqAlignment
            The alignment of the uniprot domain with the PDB domain
        str
            PDB id of the protein
        """
        # write both sequences to one file
        with open(self.unique_temp_folder + 'seqfiles.fasta', 'w') as seqFiles:
            SeqIO.write([uniprot_sequence, pdb_sequence], seqFiles, 'fasta')

        seqIDs = [uniprot_sequence.id, pdb_sequence.id]

#        self.logger.debug('Calling tcoffee with parameters:')
#        self.logger.debug('global_temp_path: {}'.format(self.global_temp_path))
#        self.logger.debug('unique_temp_path: {}'.format(self.unique_temp_folder))
        tcoffee = call_tcoffee.tcoffee_alignment(
            self.global_temp_path,
            self.unique_temp_folder,
            [self.unique_temp_folder + 'seqfiles.fasta', ],
            seqIDs,
            self.n_cores,
            self.pdb_path,
            mode,
            self.logger)
        alignments = tcoffee.align()

        return alignments[0], pdb_sequence.id



    def pick_sequence(self, alignment, pdb_sequence_id):
        """
        Pick the uniprot and pdb sequences from the alignment

        Parameters
        ----------
        alignment : list
            A list of `Bio.Align.MultipleSeqAlignment` objects.
        pdb_sequence_id : str
            The ID of the PDB sequence

        Returns
        -------
        uniprot_alignment : Bio.SeqRecord.SeqRecord
        pdb_alignemnt : Bio.SeqRecord.SeqRecord
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

        Parameters
        ----------
        domain_template : list
            A list of structural templates.

        Returns
        -------
            The best structural template.
        """

        # First sort by identity score:
        if isinstance(domain_template[0], sql_db.UniprotDomainTemplate):
            domain_template.sort(key=lambda k: k.alignment_score, reverse=True)
            max_score = domain_template[0].alignment_score
            self.logger.debug('max alignment score: %f' % max_score)
            # Collect all templates with the highest alignment score
            best_domain_templates = []
            for t in domain_template:
                self.logger.debug('alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score, t.domain.cdhit_cluster, t.domain.cdhit_cluster_idx))
                if t.alignment_score == max_score:
                    best_domain_templates.append(t)
            best_domatemplatesin_templates.sort(key=lambda k: k.domain.pdb_resolution, reverse=False)
            for t in best_domain_templates:
                self.logger.debug('Best alignments:')
                self.logger.debug('alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score, t.domain.cdhit_cluster, t.domain.cdhit_cluster_idx))

        if isinstance(domain_template[0], sql_db.UniprotDomainPairTemplate):
            domain_template.sort(key=lambda k: k.alignment_score_1 + k.alignment_score_2, reverse=True)
            max_score = domain_template[0].alignment_score_1 + domain_template[0].alignment_score_2
            self.logger.debug('max alignment score: %f' % max_score)
            # Collect all templates with the highest alignment score
            best_domain_templates = []
            for t in domain_template:
                self.logger.debug('domain 1. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_1, t.domain_1.cdhit_cluster, t.domain_1.cdhit_cluster_idx))
                self.logger.debug('domain 2. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_2, t.domain_2.cdhit_cluster, t.domain_2.cdhit_cluster_idx))
                if t.alignment_score_1 + t.alignment_score_2 == max_score:
                    best_domain_templates.append(t)
            best_domain_templates.sort(key=lambda k: k.domain_1.pdb_resolution, reverse=False)
            for t in best_domain_templates:
                self.logger.debug('Best alignments:')
                self.logger.debug('domain 1. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_1, t.domain_1.cdhit_cluster, t.domain_1.cdhit_cluster_idx))
                self.logger.debug('domain 2. alignment score: %f, cluster id: %i, cluster idx: %i' \
                    % (t.alignment_score_2, t.domain_2.cdhit_cluster, t.domain_2.cdhit_cluster_idx))
        return best_domain_templates[0]


    ###########################################################################


###############################################################################
# Obsolete

    def get_pdb_sequence(self, pdb_code, chain_id, pdb_domain_def):
        """
        Return the pdb file sequence (ATOM)
        """
        domains = [pdb_domain_def, ]
        pdb = pdb_template.PDBTemplate(self.pdb_path, pdb_code, chain_id, domains,
                                       self.unique_temp_folder, self.unique_temp_folder, self.logger)
        pdb.extract()
        chain_sequence, chain_numbering_extended = pdb.get_chain_sequence_and_numbering(chain_id)
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
            pdb_domain_def_index = chain_numbering_extended.index(
                pdb_domain_def[0])+1, chain_numbering_extended.index(pdb_domain_def[1])+1
        except ValueError:
            raise errors.PDBError(
                'ValueError when mapping domain boundaries to sequence numbering: ' + pdb_code + '_' + chain_id)

        return chain_seqrecord, pdb_domain_def_index, chain_numbering_extended, pdb_domain_def


def split_superdomains_1(superdomain):
    domains = [domain.split('_') for domain in superdomain.split('+')]
    print(domains)
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
    subprocess.check_call(
        'cd ' + tmp_path + 'blast/ && ln -sf /home/kimlab1/strokach/ncbi-blast-2.2.28+/db', shell=True)
    subprocess.check_call('cp ' + '/home/kimlab1/strokach/working/pipeline/bin/provean ' +
        tmp_path + unique + '/sequence_conservation/', shell=True)

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
    temp = Bio.Align.MultipleSeqAlignment([
        Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq('--AGGA-')),
        Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq('MMMGGMM'))])
    print(get_template.get_coverage(temp, len(temp[0])))
    print(get_template.get_identity(temp))
    print(get_template.get_interacting_identity(temp, [4,5], 2))

    p = db.get_uniprot_domain('Q8NEU8') + db.get_uniprot_domain_pair('Q8NEU8')
    p = [db.get_uniprot_domain_pair('Q8NEU8')[0]]
#    protein_domains = db.get_uniprot_domain_pair('Q8NEU8')
    for d, t, m in p:
        template = get_template(d)
#        if isinstance(template, sql_db.UniprotDomainTemplate):
#            print get_template.build_provean_supporting_set(d, template)

