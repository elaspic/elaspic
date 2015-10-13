# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import unicode_literals
from future import standard_library
standard_library.install_aliases()
from builtins import next
from builtins import zip
from builtins import range
from builtins import object
import shutil

import os
import os.path as op
import subprocess
import pickle as pickle
import six
import logging

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser

from . import (
    conf, helper, errors, sequence, structure_tools, structure_analysis, call_foldx,
    database_tables,
)

logger = logging.getLogger(__name__)
configs = conf.Configs()


#%%
class GetMutation(object):

    def __init__(self, **kwargs):
        """
        uniprot_id_1 = None
        mutation = None
        domain_sequences = None
        chains_modeller = None
        uniprot_domain_id = None ## have to get rid of this 
        path_to_provean_supset = None # provean_supset_file
        save_path = None # archive_dir ? 
        pdbFile_wt = None # pdb_file ? 
        position_domain = None # mutation_domain_pos
        mutation_domain = None # mutation_domain
        structure = None
        position_modeller = None # modeller_mutation_pos
        mutation_modeller = None # modeller_mutation
        """
        self.__dict__ = kwargs
        self.is_interface = len(self.domain_sequences) > 1
        self.validate()


    def validate(self):
        """Run some sanity checks on the input data.
        """
        mutated_aa_domain = self.domain_sequences[0][int(self.mutation_domain[1:-1])-1]
        if mutated_aa_domain != self.mutation_domain[0]:
            logger.error(
                'Domain sequence: {}'
                .format(str(self.domain_sequences[0].seq)))
            logger.error(
                'Domain AA: {};\t Mutation AA: {}'
                .format(mutated_aa_domain, self.mutation_domain[0]))
            raise Exception('Mutated amino acid was not found inside the specified domain!')


    def __call__(mutation_features):
        """
        Parameters
        ----------
        mutation_features : 
            An object that can hold all mutation feautres for core or interface.
        """
        ...


    def evaluate_mutation(self, d, mut_data, uniprot_mutation):
        """
        """
        if (mut_data.path_to_provean_supset and
            uniprot_mutation.provean_score in [None, 0, 0.0, 1.0]): # TODO: change back to None
            logger.debug('Calculating the provean score for the mutation...')
            try:
                provean_mutation, provean_score = self.get_provean_score(
                    mut_data.uniprot_domain_id, mut_data.mutation_domain,
                    mut_data.domain_sequences[0], mut_data.path_to_provean_supset)
            except errors.ProveanError as e:
                logger.error(str(type(e)) + ': ' + e.__str__())
                provean_mutation, provean_score = None, None
            logger.debug('provean mutation:')
            logger.debug(provean_mutation)
            logger.debug('provean score:')
            logger.debug(provean_score)
            uniprot_mutation.provean_score = provean_score

        if not uniprot_mutation.stability_energy_wt:
            logger.debug('Evaluating the structural impact of the mutation...')
            uniprot_mutation = self.evaluate_structural_impact(d, mut_data, uniprot_mutation)

        if (uniprot_mutation.provean_score and
            uniprot_mutation.stability_energy_wt and
            uniprot_mutation.ddg in [None, 1.0]): # TODO: Change back to None
                logger.debug('Predicting the thermodynamic effect of the mutation...')
                uniprot_mutation = self.predict_thermodynamic_effect(d, uniprot_mutation)

        return uniprot_mutation


    def get_provean_score(self, uniprot_domain_id, domain_mutation, domain_sequence, path_to_provean_supset):
        """
        """
#        uniprot_sequence = self.db.get_uniprot_sequence(d.uniprot_id).seq.tostring()
#        domain_def = database_tables.decode_domain_def(t.domain_def)
#        uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
        if isinstance(domain_sequence, six.string_types):
            pass
        elif isinstance(domain_sequence, Seq):
            domain_sequence = str(domain_sequence)
        elif isinstance(domain_sequence, SeqRecord):
            domain_sequence = str(domain_sequence.seq)
        else:
            raise Exception('Wrong class type %s for domain_sequence' % str(type(domain_sequence)))

#        first_aa = uniprot_sequence_domain[0]
    #    provean_supset_filename = t.alignment_filename.replace('.aln', '.supset')
#        provean_supset_filename = t.provean_supset_filename

        path_to_provean_supset_local = (
            conf.configs['unique_temp_folder'] + 'sequence_conservation/' +
            path_to_provean_supset.split('/')[-1]
        )
        subprocess.check_call(
            'cp -f ' + path_to_provean_supset + ' ' + path_to_provean_supset_local, shell=True)
        subprocess.check_call(
            'cp -f ' + path_to_provean_supset + '.fasta ' + path_to_provean_supset_local + '.fasta', shell=True)

        conf.configs = {
            'provean_temp_path': self.provean_temp_path,
            'global_temp_path': conf.configs['global_temp_path'],
            'n_cores': conf.configs['n_cores'],
            'path_to_provean_supset': path_to_provean_supset,
            'path_to_provean_supset_local': path_to_provean_supset_local,

        }
        result, error_message, return_code = sequence.check_provean_supporting_set(
            domain_mutation, domain_sequence, conf.configs, conf.configs['unique_temp_folder'], self.provean_temp_path,
            str(uniprot_domain_id), path_to_provean_supset_local,
            save_supporting_set=False, check_mem_usage=False)

        logger.debug(result)
        while (return_code != 0 and
                ('IDs are not matched' in error_message or 'OID not found' in error_message)):
            logger.error(error_message)
            line_to_remove = error_message.split(':')[1].split(',')[0].strip()
            logger.error('Removing line with id: {} from the supporting set...'.format(line_to_remove))
            with open(path_to_provean_supset_local) as ifh, \
                    open(path_to_provean_supset_local + '.mod', 'w') as ofh:
                for line in ifh:
                    if line_to_remove not in line:
                        ofh.write(line)
            subprocess.check_call('mv -f ' + path_to_provean_supset_local + '.mod ' + path_to_provean_supset_local, shell=True)
            result, error_message, return_code = sequence.check_provean_supporting_set(
                self, domain_mutation, domain_sequence, str(uniprot_domain_id), path_to_provean_supset_local,
                save_supporting_set=False, check_mem_usage=False)
            logger.debug(result)

        if return_code != 0:
            logger.error(error_message)
            raise errors.ProveanError(error_message)

        ### Results look something like this:
        #[23:28:34] clustering subject domain_sequences...
        #[23:28:34] selecting clusters...
        #[23:28:34] 0 subject domain_sequences in 0 clusters were selected for supporting domain_sequences.
        #[23:28:34] use the query itself as a supporting sequence
        #[23:28:34] loading subject domain_sequences from a FASTA file...
        #[23:28:34] scores were computed based on the query sequence itself.
        ## Number of clusters:	1
        ## Number of supporting domain_sequences used:	1
        #[23:28:34] computing delta alignment scores...
        #[23:28:34] printing PROVEAN scores...
        ### PROVEAN scores ##
        ## VARIATION	SCORE
        #M1A	-6.000

        # Commented out for now because we just need to get the supporting domain_sequences
        variations_started = False
        provean_mutation, provean_score = None, np.nan
        for line in result.split('\n'):
            if 'VARIATION\tSCORE' in line:
                variations_started = True
                continue
            if variations_started:
                provean_mutation, provean_score = line.split()
                provean_score = np.float(provean_score)
                break
        return provean_mutation, provean_score


    def evaluate_structural_impact(self, d, mut_data, uniprot_mutation):

        uniprot_id_1 = mut_data.uniprot_id_1
        mutation = mut_data.mutation

        #######################################################################
        # Copy the model pdb to the foldx folder
        system_command = (
            'cp -u ' + mut_data.save_path + mut_data.pdbFile_wt + ' ' +
            conf.configs['unique_temp_folder'] + 'FoldX/' + mut_data.pdbFile_wt.split('/')[-1])
        childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell=True, universal_newlines=True)
        result, error_message = childProcess.communicate()
        if childProcess.returncode != 0:
            logger.error('cp result: {}'.format(result))
            logger.error('cp error: {}'.format(error_message))

        #######################################################################
        ## 2nd: use the 'Repair' feature of FoldX to optimise the structure
        fX = call_foldx.FoldX(conf.configs['unique_temp_folder'],
                   mut_data.save_path + mut_data.pdbFile_wt,
                   mut_data.chains_modeller[0],
                   conf.configs['foldx_num_of_runs'],
                   conf.configs['foldx_water'],
                   logger)
        repairedPDB_wt = fX('RepairPDB')

        #######################################################################
        ## 3rd: introduce the mutation using FoldX
        if len(mut_data.domain_sequences) == 1:
            logger.debug(mut_data.domain_sequences[0].seq)
        else:
            logger.debug(mut_data.domain_sequences[0].seq)
            logger.debug(mut_data.domain_sequences[1].seq)
#        mutations_foldX = [prepareMutationFoldX(mut_data.domain_sequences[0], mut_data.mutation_domain),]
#        logger.debug("mutations_foldX:")
#        logger.debug(mutations_foldX)

        # compile a list of mutations
        mutCodes = [mutation[0] + mut_data.chains_modeller[0] + mut_data.position_modeller[0] + mutation[-1], ]
        logger.debug('Mutcodes for foldx:')
        logger.debug(mutCodes)

        # Introduce the mutation using foldX
        fX_wt = call_foldx.FoldX(conf.configs['unique_temp_folder'],
                      repairedPDB_wt,
                      mut_data.chains_modeller[0],
                      conf.configs['foldx_num_of_runs'],
                      conf.configs['foldx_water'],
                      logger)
        repairedPDB_wt_list, repairedPDB_mut_list = fX_wt('BuildModel', mutCodes)

        wt_chain_sequences = structure_tools.get_structure_sequences(repairedPDB_wt_list[0])
        mut_chain_sequences = structure_tools.get_structure_sequences(repairedPDB_mut_list[0])
        logger.debug('repairedPDB_wt_list: %s' % str(repairedPDB_wt_list))
        logger.debug('wt_chain_sequences: %s' % str(wt_chain_sequences))
        logger.debug('repairedPDB_mut_list: %s' % str(repairedPDB_mut_list))
        logger.debug('mut_chain_sequences: %s' % str(mut_chain_sequences))
        self.__check_structure_match(repairedPDB_wt_list[0], mut_data, mutation[0])
        self.__check_structure_match(repairedPDB_mut_list[0], mut_data, mutation[-1])

        # Copy the foldX wildtype and mutant pdb files (use the first model if there are multiple)
        model_filename_wt = uniprot_id_1 + '_' + mutation + '/' +  repairedPDB_wt_list[0].split('/')[-1]
        model_filename_mut = uniprot_id_1 + '_' + mutation + '/MUT_' +  repairedPDB_mut_list[0].split('/')[-1]
        subprocess.check_call('mkdir -p ' + mut_data.save_path + uniprot_id_1 + '_' + mutation + '/', shell=True)
        subprocess.check_call('cp ' + repairedPDB_wt_list[0] + ' ' + mut_data.save_path + model_filename_wt, shell=True)
        subprocess.check_call('cp ' + repairedPDB_mut_list[0] + ' ' + mut_data.save_path + model_filename_mut, shell=True)

        #######################################################################
        ## 4th: set up the classes for the wildtype and the mutant structures
        fX_wt_list = list()
        for wPDB in repairedPDB_wt_list:
            fX_wt_list.append(
                call_foldx.FoldX(
                    conf.configs['unique_temp_folder'], wPDB, mut_data.chains_modeller[0],
                    conf.configs['foldx_num_of_runs'], conf.configs['foldx_water']))

        fX_mut_list = list()
        for mPDB in repairedPDB_mut_list:
            fX_mut_list.append(
                call_foldx.FoldX(
                    conf.configs['unique_temp_folder'], mPDB, mut_data.chains_modeller[0],
                    conf.configs['foldx_num_of_runs'], conf.configs['foldx_water']))

        #######################################################################
        ## 5th: calculate the energy for the wildtype
        #if isinstance(d, database_tables.UniprotDomain):
        stability_values_wt = helper.encode_list_as_text([foldx('Stability') for foldx in fX_wt_list])
        stability_values_mut = helper.encode_list_as_text([foldx('Stability') for foldx in fX_mut_list])

        if isinstance(d, database_tables.UniprotDomainPair):
            complex_stability_values_wt = helper.encode_list_as_text([ foldx('AnalyseComplex') for foldx in fX_wt_list ])
            complex_stability_values_mut = helper.encode_list_as_text([ foldx('AnalyseComplex') for foldx in fX_mut_list ])

        #######################################################################
        ## 9th: calculate the pysico-chemical properties
        def get_contact_vectors(physi_chem, pdb_filename_list, mutated_chain_id, mutation_domain):
            opposite_chain_contact_vector_all = []
            same_chain_contact_vector_all = []
            for pdb_filename in pdb_filename_list:
                opposite_chain_contact_vector, same_chain_contact_vector = physi_chem(
                    pdb_filename, mutated_chain_id, mutation_domain)
                opposite_chain_contact_vector_all.append(opposite_chain_contact_vector)
                same_chain_contact_vector_all.append(same_chain_contact_vector)
            return [helper.encode_list_as_text(opposite_chain_contact_vector_all),
                    helper.encode_list_as_text(same_chain_contact_vector_all)]

        physi_chem = structure_analysis.PhysiChem(5.0, 4.0, conf.configs['unique_temp_folder'])
        opposite_chain_contact_vector_all_wt, same_chain_contact_vector_all_wt = get_contact_vectors(
            physi_chem, repairedPDB_wt_list, mut_data.chains_modeller[0], mut_data.mutation_domain)
        opposite_chain_contact_vector_all_mut, same_chain_contact_vector_all_mut = get_contact_vectors(
            physi_chem, repairedPDB_wt_list, mut_data.chains_modeller[0], mut_data.mutation_domain)


        #######################################################################
        # Calculate secondary structure, sasa, and interchain distance
        def obtain_additional_mutation_properties(repaired_pdb_list, is_domain_pair=False):
            analyze_structure_instance = structure_analysis.AnalyzeStructure(
                conf.configs['unique_temp_folder'] + 'FoldX/',
                conf.configs['unique_temp_folder'] + 'analyze_structure/',
                repaired_pdb_list[0].split('/')[-1], # dssp file wildtype
                mut_data.chains_modeller, None)

            (seasa_by_chain_together, seasa_by_chain_separately,
            seasa_by_residue_together, seasa_by_residue_separately) = analyze_structure_instance.get_seasa()
            seasa_info = seasa_by_residue_separately[
                (seasa_by_residue_separately['pdb_chain']==mut_data.chains_modeller[0]) &
                (seasa_by_residue_separately['res_num']==mut_data.position_modeller[0])].iloc[0]
            if (structure_tools.convert_aa(seasa_info['res_name']) != mut_data.mutation_domain[0] and
                structure_tools.convert_aa(seasa_info['res_name']) != mut_data.mutation_domain[-1]):
                    logger.error('Wrong amino acid for msms mutant!')
                    logger.error(structure_tools.convert_aa(seasa_info['res_name']))
                    logger.error(mut_data.mutation_domain)
                    logger.error(seasa_info)
                    logger.error(seasa_by_residue_separately)
                    raise Exception('surface area calculated for the wrong atom!')
            solvent_accessibility = seasa_info['rel_sasa']

            # Secondary structure
            secondary_structure_df = analyze_structure_instance.get_secondary_structure()
            secondary_structure_df = secondary_structure_df[
                (secondary_structure_df.chain == mut_data.chains_modeller[0]) &
                (secondary_structure_df.idx == int(mut_data.mutation_domain[1:-1]))]
            assert len(secondary_structure_df) == 1
            secondary_structure_df = secondary_structure_df.iloc[0]
            if (secondary_structure_df.amino_acid != mut_data.mutation_domain[0] and
                secondary_structure_df.amino_acid != mut_data.mutation_domain[-1]):
                    logger.error('Wrong amino acid for stride output!')
                    logger.error(secondary_structure_df.amino_acid)
                    logger.error(mut_data.mutation_domain)
                    raise Exception('surface area calculated for the wrong atom!')
            secondary_structure = secondary_structure_df.ss_code

            contact_distance = None
            if is_domain_pair:
                try:
                    contact_distance = analyze_structure_instance.get_interchain_distances(mut_data.chains_modeller[0], mut_data.mutation_modeller)
                    contact_distance = contact_distance[mut_data.chains_modeller[0]][mut_data.chains_modeller[1]]
                    logger.debug(
                        'The shortest interchain distance between chain {} and chain {} is {}'
                        .format(mut_data.chains_modeller[0], mut_data.chains_modeller[1], contact_distance))
                    if not contact_distance:
                        raise ValueError
                except (IndexError, KeyError, ValueError) as e:
                    logger.error('Could not calculate the shortest contact distance between two chains!')
                    logger.error(str(e))
                    logger.error(contact_distance)
                    raise e
            return solvent_accessibility, secondary_structure, contact_distance

        solvent_accessibility_wt, secondary_structure_wt, contact_distance_wt = \
            obtain_additional_mutation_properties(repairedPDB_wt_list, isinstance(d, database_tables.UniprotDomainPair))
        solvent_accessibility_mut, secondary_structure_mut, contact_distance_mut = \
            obtain_additional_mutation_properties(repairedPDB_mut_list, isinstance(d, database_tables.UniprotDomainPair))

        logger.debug('solvent_accessibility (wt/mut): ({}/{})'.format(solvent_accessibility_wt, solvent_accessibility_mut))
        logger.debug('secondary_structure (wt/mut): ({}/{})'.format(secondary_structure_wt, secondary_structure_mut))
        logger.debug('contact_distance (wt/mut): ({}/{})'.format(contact_distance_wt, contact_distance_mut))


        #######################################################################
        ## 5th: calculate the energy for the wildtype
        uniprot_mutation.uniprot_id = uniprot_id_1
        uniprot_mutation.mutation = mutation

        uniprot_mutation.chain_modeller = mut_data.chains_modeller[0]
        uniprot_mutation.mutation_modeller = mut_data.mutation_modeller

        uniprot_mutation.model_filename_wt = model_filename_wt
        uniprot_mutation.model_filename_mut = model_filename_mut

        uniprot_mutation.stability_energy_wt = stability_values_wt
        uniprot_mutation.stability_energy_mut = stability_values_mut

        uniprot_mutation.physchem_wt = opposite_chain_contact_vector_all_wt
        uniprot_mutation.physchem_wt_ownchain = same_chain_contact_vector_all_wt
        uniprot_mutation.physchem_mut = opposite_chain_contact_vector_all_mut
        uniprot_mutation.physchem_mut_ownchain = same_chain_contact_vector_all_mut

        uniprot_mutation.secondary_structure_wt = secondary_structure_wt
        uniprot_mutation.solvent_accessibility_wt = solvent_accessibility_wt
        uniprot_mutation.secondary_structure_mut = secondary_structure_mut
        uniprot_mutation.solvent_accessibility_mut = solvent_accessibility_mut

        if isinstance(d, database_tables.UniprotDomainPair):
            uniprot_mutation.analyse_complex_energy_wt = complex_stability_values_wt
            uniprot_mutation.analyse_complex_energy_mut = complex_stability_values_mut
            uniprot_mutation.contact_distance_wt = contact_distance_wt
            uniprot_mutation.contact_distance_mut = contact_distance_mut
        #######################################################################
        # Save alignments and modeller models to output database for storage
        # Move template files to the output folder as a tar archives
        # self.make_tarfile(self.HOME + self.outputPath + save_path + '_' + mutation + '.tar.bz2', save_path[:-1])

        #######################################################################
        logger.info('Finished processing template:')
        logger.info(mut_data.save_path.split('/')[-2])

        return uniprot_mutation


    def __check_structure_match(self, pdb_filename, mut_data, expecte_aa):
        """
        """
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        structure = parser.get_structure('ID', pdb_filename)
        model = structure[0]
        chain = model[mut_data.chains_modeller[0]]
        residue_found = False
        for residue in chain:
            if (str(residue.id[1]) + residue.id[2].strip()) == mut_data.position_modeller[0]:
                residue_found = True
                break
        if not residue_found or not (structure_tools.convert_aa(residue.resname) == expecte_aa):
            logger.error('residue_found? %s' % residue_found)
            logger.error(residue.resname)
            logger.error(residue.id)
            logger.error(mut_data.position_modeller)
            logger.error(expecte_aa)
            raise errors.FoldXAAMismatchError('Expected and actual FoldX amino acids do not match!')




