# -*- coding: utf-8 -*-

import os
import subprocess
import cPickle as pickle

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser

import domain_alignment
import errors
import sql_db
import analyze_structure
import pdb_template
import call_foldx
import helper_functions as hf


###############################################################################
# Useful functions
def _score_match(pair_match, matrix_match):
    if pair_match not in matrix_match:
        return matrix_match[(tuple(reversed(pair_match)))]
    else:
        return matrix_match[pair_match]


def score_pairwise(seq1, seq2, matrix, gap_s, gap_e):
    """ Function needed to get the BLOSUM score """
    score = 0
    gap = False
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_s
            else:
                score += _score_match(pair, matrix)
        else:
            if '-' not in pair:
                gap = False
                score += _score_match(pair, matrix)
            else:
                score += gap_e
    return score


###############################################################################
# H     Alpha helix
# G     3-10 helix
# I     PI-helix
# E     Extended conformation
# B|b   Isolated bridge
# T     Turn
# C|-   Coil (none of the above) ()
secondary_structure_to_int = {
    '-': 0, 'C': 0, 'B': 1, 'b': 1, 'E': 2, 'G': 3, 'H': 4, 'I': 5, 'S': 6, 'T': 7}


def format_mutation_features(feature_df, core_or_interface):
    if core_or_interface == False or core_or_interface == 0 or core_or_interface == 'core':
        foldx_column_name = 'stability_energy'
        foldx_feature_names_wt = call_foldx.names_stability_wt
        foldx_feature_names_mut = call_foldx.names_stability_mut
    elif core_or_interface == True or core_or_interface == 1 or core_or_interface == 'interface':
        foldx_column_name = 'analyse_complex_energy'
        foldx_feature_names_wt = call_foldx.names_stability_complex_wt
        foldx_feature_names_mut = call_foldx.names_stability_complex_mut

    # Drop rows that have missing FoldX information
    # (should not happen when callced from inside the pipeline because we have only one column)
    feature_df = feature_df.dropna(subset=[foldx_column_name + '_wt', foldx_column_name + '_mut'])

    # FoldX output
    for column_index, column_name in enumerate(foldx_feature_names_wt):
        feature_df[column_name] = feature_df[foldx_column_name + '_wt'].apply(lambda x: float(x.split(',')[column_index]))
    del feature_df[foldx_column_name + '_wt']

    for column_index, column_name in enumerate(foldx_feature_names_mut):
        feature_df[column_name] = feature_df[foldx_column_name + '_mut'].apply(lambda x: float(x.split(',')[column_index]))
    del feature_df[foldx_column_name + '_mut']

    # PhysicoChemical properties
    names_phys_chem = ['pcv_salt_equal', 'pcv_salt_opposite', 'pcv_hbond', 'pcv_vdw']
    for column_index, column_name in enumerate(names_phys_chem):
        feature_df[column_name + '_wt'] = feature_df['physchem_wt'].apply(lambda x: int(x.split(',')[column_index]))
        feature_df[column_name + '_self_wt'] = feature_df['physchem_wt_ownchain'].apply(lambda x: int(x.split(',')[column_index]))
        feature_df[column_name + '_mut'] = feature_df['physchem_mut'].apply(lambda x: int(x.split(',')[column_index]))
        feature_df[column_name + '_self_mut'] = feature_df['physchem_mut_ownchain'].apply(lambda x: int(x.split(',')[column_index]))
    del feature_df['physchem_wt']
    del feature_df['physchem_wt_ownchain']
    del feature_df['physchem_mut']
    del feature_df['physchem_mut_ownchain']

    for col in feature_df.columns:
        if 'secondary_structure' in col:
            feature_df[col] = feature_df[col].apply(lambda x: secondary_structure_to_int[x])
    return feature_df


def get_mutation_features(d, mut, row_idx=0):
    """
    """
    feature_dict = {key: value for (key, value) in mut.__dict__.iteritems() if not key.startswith('_')}

    feature_dict.update({
        # Header columns
        # 'uniprot_id': mut.uniprot_id,
        # 'mutation': mut.mutation,
        't_date_modified': d.template.t_date_modified,
        'm_date_modified': d.template.model.m_date_modified,
        # 'mut_date_modified': mut.mut_date_modified,
        #
        'norm_dope': d.template.model.norm_dope,
    })

    if hasattr(d, 'uniprot_domain_id'):
        feature_dict.update({
            # Header columns
            'uniprot_domain_id': d.uniprot_domain_id,
            'cath_id': d.template.cath_id,
            'pfam_name': d.pdbfam_name,
            'clan_name': d.pfam_clan,
            #
            'alignment_identity': d.template.alignment_identity,
            'alignment_coverage': d.template.alignment_coverage,
            'alignment_score': d.template.alignment_score,
        })

    elif hasattr(d, 'uniprot_domain_pair_id'):
        feature_dict.update({
            # Header columns
            'uniprot_domain_pair_id': d.uniprot_domain_pair_id,
            'cath_id_1': d.template.cath_id_1,
            'cath_id_2': d.template.cath_id_2,
            'pfam_name_1': d.uniprot_domain_1.pdbfam_name,
            'pfam_name_2': d.uniprot_domain_2.pdbfam_name,
            'clan_name_1': d.uniprot_domain_1.pfam_clan,
            'clan_name_2': d.uniprot_domain_2.pfam_clan,
            #
            'alignment_identity': d.template.identical_1 + d.template.identical_2,
            'identical_1': d.template.identical_1,
            'identical_2': d.template.identical_2,
            'interface_area_hydrophobic': d.template.model.interface_area_hydrophobic,
            'interface_area_hydrophilic': d.template.model.interface_area_hydrophilic,
            'interface_area_total': d.template.model.interface_area_total,
        })

    feature_df = pd.DataFrame(feature_dict, index=[row_idx])

    return feature_df


def convert_features_to_differences(df, keep_mut=False):
    """
    Return a dataframe with all the `_mut` features replaced with `_change`
    features which coresspond to the difference between `_mut` and `_wt`.
    """
    column_list = []
    for column_name, column in df.iteritems():
        if ('_mut' in column_name and
            column_name.replace('_mut','_wt') in df.columns and
            df[column_name].dtype != object):
                if keep_mut:
                    column_list.append(column)
                new_column = column - df[column_name.replace('_mut','_wt')]
                if 'secondary_structure' in column_name:
                    new_column = new_column.apply(lambda x: 1 if x else 0)
                new_column.name = column_name.replace('_mut','_change')
                column_list.append(new_column)
        else:
            column_list.append(column)
#    new_df = pd.DataFrame(column_list).T
    new_df = pd.concat(column_list, axis=1)
    return new_df


def encode_list_as_text(list_of_lists):
    return ','.join([':'.join([str(x) for x in xx]) for xx in zip(*list_of_lists)])


def decode_text_as_list(list_string):
    str2num = lambda x: float(x) if '.' in x else int(x)
    return zip(*[[str2num(x) for x in sublist.split(':')] for sublist in list_string.split(',')])



###############################################################################
class MutationData(object):

    def __init__(self):
        self.uniprot_id_1 = None
        self.mutation = None
        self.pfam_name = None
        self.domain_start = None
        self.domain_end = None
        self.alignment = None
        self.alignment_id = None
        self.chains_pdb = None
        self.uniprot_sequences = None
        self.domain_sequences = None
        self.HETflag = None
        self.chains_modeller = None
        self.uniprot_domain_id = None
        self.path_to_provean_supset = None
        self.save_path = None
        self.pdbFile_wt = None
        self.position_domain = None
        self.mutation_domain = None
        self.structure = None
        self.position_modeller = None
        self.mutation_modeller = None


class GetMutation(object):

    def __init__(
            self, global_temp_path, temp_path, unique, pdb_path,
            db, logger, n_cores, bin_path, provean_temp_path, foldX_WATER, build_model_runs, matrix, gap_s, gap_e):
        self.global_temp_path = global_temp_path
        self.temp_path = temp_path
        self.unique = unique
        self.unique_temp_folder = temp_path + unique + '/'
        self.pdb_path = pdb_path
        self.db = db # sql database
        self.logger = logger
        self.bin_path = bin_path
        self.foldX_WATER = foldX_WATER
        self.build_model_runs = build_model_runs
        self.matrix = matrix
        self.gap_s = gap_s
        self.gap_e = gap_e
        self.n_cores = n_cores
        self.provean_temp_path = provean_temp_path

        self.clf_domain = pickle.load(open(bin_path + 'machine_learning_protherm_clf.pickle', 'rb'))
        self.clf_domain_features = pd.read_pickle(bin_path + 'machine_learning_protherm_features.pickle')
#        self.clf_interface = pickle.load(open(bin_path + 'clf_interface.pickle', 'rb'))
#        self.clf_interface_features = pd.read_pickle(bin_path + 'clf_interface_features.pickle')


    def get_mutation_data(self, d, uniprot_id_1, mutation):
        """
        """
        path_to_provean_supset = None
        #######################################################################
        # Core
        if isinstance(d, sql_db.UniprotDomain):
            d_1 = d
            self.logger.debug("Analyzing core mutation for uniprot: %s" % uniprot_id_1)
            pdbfam_name = d.pdbfam_name
            if d.template.model.model_domain_def != None:
                self.logger.debug('Using model domain definitions')
                domain_start, domain_end = hf.decode_domain_def(d.template.model.model_domain_def)
            else:
                self.logger.debug('Using template domain definitions')
                domain_start, domain_end = hf.decode_domain_def(d.template.domain_def)
            alignment, __ = self.db.get_alignment(d.template.model, d.path_to_data)
#            alignment_id = d.template.model.alignment_id
            chains_pdb = [d.template.domain.pdb_chain, ]
            uniprot_sequences = [self.db.get_uniprot_sequence(uniprot_id_1), ]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],]
            chains_modeller = [d.template.model.chain, ]
            uniprot_domain_id = d.uniprot_domain_id

        #######################################################################
        # Interface
        elif isinstance(d, sql_db.UniprotDomainPair):
            if uniprot_id_1 == d.uniprot_domain_1.uniprot_id:
                uniprot_id_2 = d.uniprot_domain_2.uniprot_id
                d_1, d_2 = d.uniprot_domain_1, d.uniprot_domain_2
                #
                if (d.template.model.model_domain_def_1 != None and
                    d.template.model.model_domain_def_2 != None):
                        domain_start, domain_end = hf.decode_domain_def(d.template.model.model_domain_def_1)
                        domain_2_start, domain_2_end = hf.decode_domain_def(d.template.model.model_domain_def_2)
                else:
                        domain_start, domain_end = hf.decode_domain_def(d.uniprot_domain_1.template.domain_def)
                        domain_2_start, domain_2_end = hf.decode_domain_def(d.uniprot_domain_2.template.domain_def)

                alignment, __ = self.db.get_alignment(d.template.model, d.path_to_data)
#                alignment_id = d.alignment_id_1
                chains_pdb = [d.template.domain_1.pdb_chain, d.template.domain_2.pdb_chain]
                interacting_aa = [int(uniprot_num) for uniprot_num in d.template.model.interacting_aa_1.split(',') if uniprot_num]
                chains_modeller = [d.template.model.chain_1, d.template.model.chain_2]

            elif uniprot_id_1 == d.uniprot_domain_2.uniprot_id:
                self.logger.debug('Mutated uniprot is uniprot 2. Rearranging...')
                uniprot_id_2 = d.uniprot_domain_1.uniprot_id
                d_1, d_2 = d.uniprot_domain_2, d.uniprot_domain_1

                if (d.template.model.model_domain_def_1 != None and
                    d.template.model.model_domain_def_2 != None):
                        domain_start, domain_end = hf.decode_domain_def(d.template.model.model_domain_def_2)
                        domain_2_start, domain_2_end = hf.decode_domain_def(d.template.model.model_domain_def_1)
                else:
                        domain_start, domain_end = hf.decode_domain_def(d.uniprot_domain_2.template.domain_def)
                        domain_2_start, domain_2_end = hf.decode_domain_def(d.uniprot_domain_1.template.domain_def)

                __, alignment = self.db.get_alignment(d.template.model, d.path_to_data)
#                alignment_id = d.alignment_id_2
                chains_pdb = [d.template.domain_2.pdb_chain, d.template.domain_1.pdb_chain]
                interacting_aa = [int(uniprot_num) for uniprot_num in d.template.model.interacting_aa_2.split(',') if uniprot_num]
                chains_modeller = [d.template.model.chain_2, d.template.model.chain_1]

            self.logger.debug('Analysing interface mutation between uniprots %s and %s' % (uniprot_id_1, uniprot_id_2,))
            uniprot_sequences = [self.db.get_uniprot_sequence(d_1.uniprot_id),
                                 self.db.get_uniprot_sequence(d_2.uniprot_id)]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],
                                uniprot_sequences[1][domain_2_start-1:domain_2_end]]

            if int(mutation[1:-1]) not in interacting_aa:
                raise errors.MutationOutsideInterfaceError('mutated residue not involved in the interaction')

        #######################################################################
        # Common
        if int(mutation[1:-1]) < domain_start or int(mutation[1:-1]) > domain_end:
            raise errors.MutationOutsideDomainError('Mutation falls outside domain')

        save_path = self.temp_path + d.path_to_data
        pdbFile_wt = d.template.model.model_filename
        pdbfam_name = d_1.pdbfam_name
        uniprot_domain_id = d_1.uniprot_domain_id

        # Provean
        if (d_1.uniprot_sequence.provean and
            d_1.uniprot_sequence.provean.provean_supset_filename):
                path_to_provean_supset = (
                    self.temp_path + hf.get_uniprot_base_path(d_1) +
                    d_1.uniprot_sequence.provean.provean_supset_filename )
                if not os.path.isfile(path_to_provean_supset):
                    error_message = (
                        'Provean supporting set sequence does not exist even though it should!\n{}'
                        .format(path_to_provean_supset)
                    )
                    self.logger.error(error_message)
                    self.logger.error('d_1: {}'.format(d_1))
                    self.logger.error('listdir: {}'.format(os.listdir(os.path.dirname(path_to_provean_supset))))
                    raise Exception(error_message)


        # Convert mutation from uniprot domain numbering to pdb domain numbering
        position_domain = int(mutation[1:-1]) - domain_start + 1 # +1 to have the numbering staring with 1
        mutation_domain = mutation[0] + str(position_domain) + mutation[-1]
#        mutation_structure = map_mutation_to_structure(alignment, alignment_id, mutation_domain, uniprot_id_1)
#        chain_mutation_structure = chains_pdb[0] + '_' + mutation_structure

        self.logger.debug('Modeller pdb file: {}'.format(save_path + pdbFile_wt))
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        structure = parser.get_structure('ID', save_path + pdbFile_wt)
        position_modeller = pdb_template.convert_position_to_resid(structure[0][chains_modeller[0]], [position_domain])
        mutation_modeller = mutation[0] + position_modeller[0] + mutation[-1]

        # Save the results
        mut_data = MutationData()

        mut_data.uniprot_id_1 = uniprot_id_1
        mut_data.mutation = mutation
        mut_data.pdbfam_name = pdbfam_name
        mut_data.domain_start = domain_start
        mut_data.domain_end = domain_end
        mut_data.alignment = alignment
#        mut_data.alignment_id = alignment_id
        mut_data.chains_pdb = chains_pdb
        mut_data.uniprot_sequences = uniprot_sequences
        mut_data.domain_sequences = domain_sequences
        mut_data.chains_modeller = chains_modeller
        mut_data.uniprot_domain_id = uniprot_domain_id
        mut_data.path_to_provean_supset = path_to_provean_supset
        mut_data.save_path = save_path
        mut_data.pdbFile_wt = pdbFile_wt
        mut_data.position_domain = position_domain
        mut_data.mutation_domain = mutation_domain
#        mut_data.structure = structure
        mut_data.position_modeller = position_modeller
        mut_data.mutation_modeller = mutation_modeller

        for key, value in mut_data.__dict__.iteritems():
            self.logger.debug(key + ':')
            self.logger.debug(value)

        # Run some sanity checks
        mutated_aa_uniprot = mut_data.uniprot_sequences[0][int(mutation[1:-1])-1]
        if mutated_aa_uniprot != mutation[0]:
            self.logger.error('Uniprot sequence: {}'.format(uniprot_sequences[0]))
            self.logger.error('Uniprot AA: {};\t Mutation AA: {}'.format(mutated_aa_uniprot, mutation[0]))
            raise Exception('Mutated amino acid was not found inside the specified uniprot!')

        mutated_aa_domain = str(mut_data.domain_sequences[0].seq)[int(mut_data.mutation_domain[1:-1])-1]
        if mutated_aa_domain != mut_data.mutation_domain[0]:
            self.logger.error('Domain sequence: {}'.format(str(mut_data.domain_sequences[0].seq)))
            self.logger.error('Domain AA: {};\t Mutation AA: {}'.format(mutated_aa_domain, mut_data.mutation_domain[0]))
            raise Exception('Mutated amino acid was not found inside the specified domain!')

        return mut_data


    def evaluate_mutation(self, d, mut_data, uniprot_mutation):
        """
        """
        if (mut_data.path_to_provean_supset and
            uniprot_mutation.provean_score in {None, 1.0}): # TODO: change back to None
            self.logger.debug('Calculating the provean score for the mutation...')
            try:
                provean_mutation, provean_score = self.get_provean_score(
                    mut_data.uniprot_domain_id, mut_data.mutation_domain,
                    mut_data.domain_sequences[0], mut_data.path_to_provean_supset)
            except errors.ProveanError as e:
                self.logger.error(str(type(e)) + ': ' + e.__str__())
                provean_mutation, provean_score = None, None
            self.logger.debug('provean mutation:')
            self.logger.debug(provean_mutation)
            self.logger.debug('provean score:')
            self.logger.debug(provean_score)
            uniprot_mutation.provean_score = provean_score

        if not uniprot_mutation.stability_energy_wt:
            self.logger.debug('Evaluating the structural impact of the mutation...')
            uniprot_mutation = self.evaluate_structural_impact(d, mut_data, uniprot_mutation)

        if (uniprot_mutation.provean_score and
            uniprot_mutation.stability_energy_wt and
            uniprot_mutation.ddg in {None, 1.0}): # TODO: Change back to None
                self.logger.debug('Predicting the thermodynamic effect of the mutation...')
                uniprot_mutation = self.predict_thermodynamic_effect(d, uniprot_mutation)

        return uniprot_mutation


    def get_provean_score(self, uniprot_domain_id, domain_mutation, domain_sequence, path_to_provean_supset):
        """
        """
#        uniprot_sequence = self.db.get_uniprot_sequence(d.uniprot_id).seq.tostring()
#        domain_def = sql_db.decode_domain_def(t.domain_def)
#        uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
        if isinstance(domain_sequence, str):
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

        path_to_provean_supset_local = self.unique_temp_folder + 'sequence_conservation/' + path_to_provean_supset.split('/')[-1]
        subprocess.check_call('cp -f ' + path_to_provean_supset + ' ' + path_to_provean_supset_local, shell=True)
        subprocess.check_call('cp -f ' + path_to_provean_supset + '.fasta ' + path_to_provean_supset_local + '.fasta', shell=True)

        result, error_message, return_code = domain_alignment.check_provean_supporting_set(
            self, domain_mutation, domain_sequence, str(uniprot_domain_id), path_to_provean_supset_local,
            save_supporting_set=False, check_mem_usage=False)

        self.logger.debug(result)
        while (return_code != 0 and
                ('IDs are not matched' in error_message or 'OID not found' in error_message)):
            self.logger.error(error_message)
            line_to_remove = error_message.split(':')[1].split(',')[0].strip()
            self.logger.error('Removing line with id: {} from the supporting set...'.format(line_to_remove))
            with open(path_to_provean_supset_local) as ifh, \
                    open(path_to_provean_supset_local + '.mod', 'w') as ofh:
                for line in ifh:
                    if line_to_remove not in line:
                        ofh.write(line)
            subprocess.check_call('mv -f ' + path_to_provean_supset_local + '.mod ' + path_to_provean_supset_local, shell=True)
            result, error_message, return_code = domain_alignment.check_provean_supporting_set(
                self, domain_mutation, domain_sequence, str(uniprot_domain_id), path_to_provean_supset_local,
                save_supporting_set=False, check_mem_usage=False)
            self.logger.debug(result)

        if return_code != 0:
            self.logger.error(error_message)
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
        provean_mutation, provean_score = None, None
        for line in result.split('\n'):
            if 'VARIATION\tSCORE' in line:
                variations_started = True
                continue
            if variations_started:
                provean_mutation, provean_score = line.split()
                break
        return provean_mutation, provean_score


    def evaluate_structural_impact(self, d, mut_data, uniprot_mutation):

        uniprot_id_1 = mut_data.uniprot_id_1
        mutation = mut_data.mutation

        #######################################################################
        # Copy the model pdb to the foldx folder
        system_command = (
            'cp -u ' + mut_data.save_path + mut_data.pdbFile_wt + ' ' +
            self.unique_temp_folder + 'FoldX/' + mut_data.pdbFile_wt.split('/')[-1])
        childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell=True)
        result, e = childProcess.communicate()
        if childProcess.returncode != 0:
            self.logger.error('cp result: {}'.format(result))
            self.logger.error('cp error: {}'.format(e))

        #######################################################################
        ## 2nd: use the 'Repair' feature of FoldX to optimise the structure
        fX = call_foldx.FoldX(self.unique_temp_folder,
                   mut_data.save_path + mut_data.pdbFile_wt,
                   mut_data.chains_modeller[0],
                   self.build_model_runs,
                   self.foldX_WATER,
                   self.logger)
        repairedPDB_wt = fX('RepairPDB')

        #######################################################################
        ## 3rd: introduce the mutation using FoldX
        if len(mut_data.domain_sequences) == 1:
            self.logger.debug(mut_data.domain_sequences[0].seq)
        else:
            self.logger.debug(mut_data.domain_sequences[0].seq)
            self.logger.debug(mut_data.domain_sequences[1].seq)
#        mutations_foldX = [prepareMutationFoldX(mut_data.domain_sequences[0], mut_data.mutation_domain),]
#        self.logger.debug("mutations_foldX:")
#        self.logger.debug(mutations_foldX)

        # compile a list of mutations
        mutCodes = [mutation[0] + mut_data.chains_modeller[0] + mut_data.position_modeller[0] + mutation[-1], ]
        self.logger.debug('Mutcodes for foldx:')
        self.logger.debug(mutCodes)

        # Introduce the mutation using foldX
        fX_wt = call_foldx.FoldX(self.unique_temp_folder,
                      repairedPDB_wt,
                      mut_data.chains_modeller[0],
                      self.build_model_runs,
                      self.foldX_WATER,
                      self.logger)
        repairedPDB_wt_list, repairedPDB_mut_list = fX_wt('BuildModel', mutCodes)

        wt_chain_sequences = pdb_template.get_structure_sequences(repairedPDB_wt_list[0])
        mut_chain_sequences = pdb_template.get_structure_sequences(repairedPDB_mut_list[0])
        self.logger.debug('repairedPDB_wt_list: %s' % str(repairedPDB_wt_list))
        self.logger.debug('wt_chain_sequences: %s' % str(wt_chain_sequences))
        self.logger.debug('repairedPDB_mut_list: %s' % str(repairedPDB_mut_list))
        self.logger.debug('mut_chain_sequences: %s' % str(mut_chain_sequences))
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
            fX_wt_list.append(call_foldx.FoldX(self.unique_temp_folder,
                                    wPDB,
                                    mut_data.chains_modeller[0],
                                    self.build_model_runs,
                                    self.foldX_WATER,
                                    self.logger))

        fX_mut_list = list()
        for mPDB in repairedPDB_mut_list:
            fX_mut_list.append(call_foldx.FoldX(self.unique_temp_folder,
                                     mPDB,
                                     mut_data.chains_modeller[0],
                                     self.build_model_runs,
                                     self.foldX_WATER,
                                     self.logger))

        #######################################################################
        ## 5th: calculate the energy for the wildtype
        #if isinstance(d, sql_db.UniprotDomain):
        stability_values_wt = encode_list_as_text([foldx('Stability') for foldx in fX_wt_list])
        stability_values_mut = encode_list_as_text([foldx('Stability') for foldx in fX_mut_list])

        if isinstance(d, sql_db.UniprotDomainPair):
            complex_stability_values_wt = encode_list_as_text([ foldx('AnalyseComplex') for foldx in fX_wt_list ])
            complex_stability_values_mut = encode_list_as_text([ foldx('AnalyseComplex') for foldx in fX_mut_list ])

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
            return [encode_list_as_text(opposite_chain_contact_vector_all),
                    encode_list_as_text(same_chain_contact_vector_all)]

        physi_chem = analyze_structure.PhysiChem(5.0, 4.0, self.unique_temp_folder, self.logger)
        opposite_chain_contact_vector_all_wt, same_chain_contact_vector_all_wt = get_contact_vectors(
            physi_chem, repairedPDB_wt_list, mut_data.chains_modeller[0], mut_data.mutation_domain)
        opposite_chain_contact_vector_all_mut, same_chain_contact_vector_all_mut = get_contact_vectors(
            physi_chem, repairedPDB_wt_list, mut_data.chains_modeller[0], mut_data.mutation_domain)


        #######################################################################
        # Calculate secondary structure, sasa, and interchain distance
        def obtain_additional_mutation_properties(repaired_pdb_list, is_domain_pair=False):
            analyze_structure_instance = analyze_structure.AnalyzeStructure(
                self.unique_temp_folder + 'FoldX/',
                self.unique_temp_folder + 'analyze_structure/',
                repaired_pdb_list[0].split('/')[-1], # dssp file wildtype
                mut_data.chains_modeller, None, self.logger)

            (seasa_by_chain_together, seasa_by_chain_separately,
            seasa_by_residue_together, seasa_by_residue_separately) = analyze_structure_instance.get_seasa()
            seasa_info = seasa_by_residue_separately[
                (seasa_by_residue_separately['pdb_chain']==mut_data.chains_modeller[0]) &
                (seasa_by_residue_separately['res_num']==mut_data.position_modeller[0])].iloc[0]
            if (pdb_template.convert_aa(seasa_info['res_name']) != mut_data.mutation_domain[0] and
                pdb_template.convert_aa(seasa_info['res_name']) != mut_data.mutation_domain[-1]):
                    self.logger.error('Wrong amino acid for msms mutant!')
                    self.logger.error(pdb_template.convert_aa(seasa_info['res_name']))
                    self.logger.error(mut_data.mutation_domain)
                    self.logger.error(seasa_info)
                    self.logger.error(seasa_by_residue_separately)
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
                    self.logger.error('Wrong amino acid for stride output!')
                    self.logger.error(secondary_structure_df.amino_acid)
                    self.logger.error(mut_data.mutation_domain)
                    raise Exception('surface area calculated for the wrong atom!')
            secondary_structure = secondary_structure_df.ss_code

            contact_distance = None
            if is_domain_pair:
                try:
                    contact_distance = analyze_structure_instance.get_interchain_distances(mut_data.chains_modeller[0], mut_data.mutation_modeller)
                    contact_distance = contact_distance[mut_data.chains_modeller[0]][mut_data.chains_modeller[1]]
                    self.logger.debug(
                        'The shortest interchain distance between chain {} and chain {} is {}'
                        .format(mut_data.chains_modeller[0], mut_data.chains_modeller[1], contact_distance))
                    if not contact_distance:
                        raise ValueError
                except (IndexError, KeyError, ValueError) as e:
                    self.logger.error('Could not calculate the shortest contact distance between two chains!')
                    self.logger.error(str(e))
                    self.logger.error(contact_distance)
                    raise e
            return solvent_accessibility, secondary_structure, contact_distance

        solvent_accessibility_wt, secondary_structure_wt, contact_distance_wt = \
            obtain_additional_mutation_properties(repairedPDB_wt_list, isinstance(d, sql_db.UniprotDomainPair))
        solvent_accessibility_mut, secondary_structure_mut, contact_distance_mut = \
            obtain_additional_mutation_properties(repairedPDB_mut_list, isinstance(d, sql_db.UniprotDomainPair))

        self.logger.debug('solvent_accessibility (wt/mut): ({}/{})'.format(solvent_accessibility_wt, solvent_accessibility_mut))
        self.logger.debug('secondary_structure (wt/mut): ({}/{})'.format(secondary_structure_wt, secondary_structure_mut))
        self.logger.debug('contact_distance (wt/mut): ({}/{})'.format(contact_distance_wt, contact_distance_mut))

        #######################################################################
        ## 11th: get the BLOSUM (or what ever matrix is given) score
        # self.mutationss is a list containing the mutations of the form A_H56T
        matrix_score = 0
        for mut in [mut_data.mutation_domain]:
            fromAA = mut[0]
            toAA   = mut[-1]
            matrix_score += score_pairwise(fromAA, toAA, self.matrix, self.gap_s, self.gap_e)

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

        uniprot_mutation.matrix_score = matrix_score

        uniprot_mutation.secondary_structure_wt = secondary_structure_wt
        uniprot_mutation.solvent_accessibility_wt = solvent_accessibility_wt
        uniprot_mutation.secondary_structure_mut = secondary_structure_mut
        uniprot_mutation.solvent_accessibility_mut = solvent_accessibility_mut

        if isinstance(d, sql_db.UniprotDomainPair):
            uniprot_mutation.analyse_complex_energy_wt = complex_stability_values_wt
            uniprot_mutation.analyse_complex_energy_mut = complex_stability_values_mut
            uniprot_mutation.contact_distance_wt = contact_distance_wt
            uniprot_mutation.contact_distance_mut = contact_distance_mut
        #######################################################################
        # Save alignments and modeller models to output database for storage
        # Move template files to the output folder as a tar archives
        # self.make_tarfile(self.HOME + self.outputPath + save_path + '_' + mutation + '.tar.bz2', save_path[:-1])

        #######################################################################
        self.logger.info('Finished processing template:')
        self.logger.info(mut_data.save_path.split('/')[-2])

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
        if not residue_found or not (pdb_template.convert_aa(residue.resname) == expecte_aa):
            self.logger.error('residue_found? %s' % residue_found)
            self.logger.error(residue.resname)
            self.logger.error(residue.id)
            self.logger.error(mut_data.position_modeller)
            self.logger.error(expecte_aa)
            raise Exception('Expected and actual FoldX amino acids do not match!')


    def predict_thermodynamic_effect(self, d, uniprot_mutation):
        """
        """
        if isinstance(d, sql_db.UniprotDomain):
            clf = self.clf_domain
            clf_features = self.clf_domain_features
        elif isinstance(d, sql_db.UniprotDomainPair):
            # TODO: FIX THIS
            #clf = self.clf_interface
            #clf_features = self.clf_interface_features
            uniprot_mutation.ddg = 1.0
            return uniprot_mutation.ddg

        feature_name_conversion = {
            'normDOPE': 'norm_dope',
            'seq_id_avg': 'alignment_identity'}
        clf_features = [feature_name_conversion.get(x, x) for x in clf_features]

        feature_df = get_mutation_features(d, uniprot_mutation)
        feature_df = format_mutation_features(feature_df, isinstance(d, sql_db.UniprotDomainPair))
        feature_df = convert_features_to_differences(feature_df, True) # keep mut, remove it in next step
        feature_df = feature_df[clf_features]

        uniprot_mutation.ddg = clf.predict(feature_df)[0]
        self.logger.debug('Predicted ddG: {}'.format(uniprot_mutation.ddg))

        return uniprot_mutation


###############################################################################
# Methods when the input is a raw crystal structure
# They are deprecated and do not work anymore...

    def get_pdb_and_mutations(self, mutations):
        # mutations is of the form I_V70A (2VIR_AB_C_TC131I, 1PGA_A__VA29P)
        pdbCode, chains_pdb1, chains_pdb2, mutations = mutations.split('_')
        pdb_type, pdb_resolution = self.pdb_resolution_database(pdbCode)

        # Equivalent to class_get_uniprot_template_core_and_interface.py output
        template = {'uniprot_id_1': '',
                    'uniprot_id_2': '',
                    'pfam_name_1': '',
                    'pfam_name_2': '',
                    'domain_def_1': '',
                    'domain_def_2': '',
                    'pdb_id': pdbCode,
                    'pdb_type': pdb_type,
                    'pdb_resolution': pdb_resolution,
                    'pdb_chain_1': chains_pdb1,
                    'pdb_chain_2': chains_pdb2,
                    'pdb_domain_def_1': tuple(),
                    'pdb_domain_def_2': tuple(),
                    'uniprot_domain_sequence_1': '',
                    'uniprot_domain_sequence_2': '',
                    'alignment_1': None,
                    'alignment_2': None,
                    'alignment_scores': (100, 100, 100,)}

        # Unique template identifier for the template
        saveFolder = (template['pdb_id'] + '_' + '-'.join(template['chains_pdb']))
        save_path =  self.unique_temp_folder + '/' + saveFolder + '/'
        subprocess.check_call('mkdir -p ' + save_path, shell=True)

        chains_get_pdb = [ item for item in chains_pdb1 ]
        chains_get_pdb.extend( [ item for item in chains_pdb2 ] )

        normDOPE_wt, pdbFile_wt, SWITCH_CHAIN, chains_get_pdb = self._getCrystalStructure(pdbCode, chains_get_pdb, save_path)

        template_mutations = {'mutation_positions': [], 'mutations': [], 'mutations_pdb': [], 'mutations_modeller': []}
        for mutation in mutations.split(','):
    #        mutations_pdb = [mutation[1] + '_' + mutation[0] + mutation[2:], ] # So its [Chain_AAwt10AAmut, ]
            # The section below is required for making foldx_mutations
            sequence, chainNumbering = self._get_pdb_sequence(pdbCode, mutation[1])
            mutation_pos = chainNumbering.index(int(mutation[2:-1]) + 1)  # +1 to have it start with one
            mutation_pdb = mutation[1] + '_' + mutation[0] + str(mutation_pos) + mutation[-1]

            template_mutations['mutation_positions'].append(int(mutation[2:-1]))
            template_mutations['mutations'].append(mutation)
            template_mutations['mutations_pdb'].append(mutation_pdb)
            template_mutations['mutations_modeller'].append(mutation_pdb)

        # Equivalent to class_multi.py additions
        template['is_in_core'] = True
        template['normDOPE_wt'] = normDOPE_wt
        template['saveFolder'] = saveFolder
        template['save_path'] = save_path
        template['pdbFile_wt'] = pdbFile_wt
        template['domain_sequences'] = [sequence,]
        template['chains_modeller'] = chains_get_pdb
        template['switch_chain'] = False

        return template, template_mutations



    def _get_pdb_sequence(self, pdbCode, chain):
        """
        Return the pdb file sequence (not SEQRES)
        """
        pdb = pdb_template.PDBTemplate(self.pdb_path, pdbCode, chain, [], self.unique_temp_folder, self.unique_temp_folder, self.logger)

        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()
        __, chain_numbering = pdb.get_chain_sequence_and_numbering(chain)
        if chain_numbering == []:
            raise errors.PDBError('Could not get the pdb numbering for ' + pdbCode + '_' + chain)

        return next(SeqIO.parse(self.unique_temp_folder + '/' + pdbCode + chain + '.seq.txt', 'fasta')), chain_numbering



    def _getCrystalStructure(self, pdbCode, chains, savePDB, FULL_PATH=False):
        pdb = pdb_template.PDBTemplate(self.pdb_path, pdbCode, chains, [], savePDB, self.unique_temp_folder, self.logger)
        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()

        SWITCH_CHAIN = False
        if chains != chains_pdb_order:
            SWITCH_CHAIN = True
            chains = chains_pdb_order

        return '-', pdbCode.split('.')[0] + ''.join(chains) + '.pdb', SWITCH_CHAIN, chains



