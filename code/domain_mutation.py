# -*- coding: utf-8 -*-

import os
import subprocess
import cPickle as pickle
import datetime

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser

import errors
import sql_db
import analyze_structure
import pdb_template
from pdb_template import PDBTemplate
import call_foldx
import helper_functions as hf


###############################################################################
# Useful functions

def score_pairwise(seq1, seq2, matrix, gap_s, gap_e):
    """ Required to get the BLOSUM score
    """
    def score_match(pair_match, matrix_match):
        if pair_match not in matrix_match:
            return matrix_match[(tuple(reversed(pair_match)))]
        else:
            return matrix_match[pair_match]

    score = 0
    gap = False
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_s
            else:
                score += score_match(pair, matrix)
        else:
            if '-' not in pair:
                gap = False
                score += score_match(pair, matrix)
            else:
                score += gap_e
    return score


###############################################################################
def get_mutation_feature_vector(d, t, m, mut):
    """
    Return two dataframes for a supplied mutation (sql_db.UniprotDomainMutation
    or sql_db.UniprotDomainPairMutation object). The first dataframe contains the
    header information for that mutation. The second dataframe contains all
    the features to be used in machine learning.
    """
    secondary_structure_to_int = {
        '-': 0, 'B': 1, 'E': 2, 'G': 3, 'H': 4, 'I': 5, 'S': 6, 'T': 7}

    def append_to_df(name_list_value_list, dfs, existing_columns=None):
        name_list, value_list = name_list_value_list[:]
        if isinstance(value_list, basestring):
            value_list = decode_text_as_list(value_list)[0]
        name_list = list(name_list)
        value_list = list(value_list)
        for i in range(len(value_list)):
            if isinstance(value_list[i], datetime.datetime):
                value_list[i] = value_list[i].strftime('%Y-%m-%d %H-%M-%S-%f')
        if existing_columns is not None:
            column_idx = 0
            while column_idx < len(name_list):
                if name_list[column_idx] in existing_columns:
                    # print 'Removing duplicate column {}'.format(name_list[column_idx])
                    del name_list[column_idx]
                    del value_list[column_idx]
                    continue
                column_idx += 1
            existing_columns.update(name_list)
        dfs.append(pd.DataFrame([value_list], columns=name_list))

    header_dfs = []
    def append_to_header_df(name_list_value_list):
        append_to_df(name_list_value_list, header_dfs)

    feature_dfs = []
    feature_existing_columns = set()
    def append_to_feature_df(name_list_value_list):
        append_to_df(name_list_value_list, feature_dfs, feature_existing_columns)

    ### Header information ----------------------------------------------------
    header_common = [
        ['uniprot_id', mut.uniprot_id],
        ['mutation', mut.mutation],
        ['t_date_modified', t.t_date_modified],
        ['m_date_modified', m.m_date_modified],
        ['mut_date_modified', mut.mut_date_modified]]
    append_to_header_df(zip(*header_common))

    if isinstance(d, sql_db.UniprotDomain):
        header_domain = [
            ['cath_id_1', t.domain.cath_id],
            ['pfam_name_1', d.pfam_name],
            ['clan_name_1', d.clan_name],
            ['uniprot_domain_id', d.uniprot_domain_id],]
        append_to_header_df(zip(*header_domain))

    elif isinstance(d, sql_db.UniprotDomainPair):
        header_domain_contact = [
            ['cath_id_1', t.domain_1.cath_id],
            ['cath_id_2', t.domain_2.cath_id],
            ['pfam_name_1', d.uniprot_domain_1.pfam_name],
            ['pfam_name_2', d.uniprot_domain_2.pfam_name],
            ['clan_name_1', d.uniprot_domain_1.clan_name],
            ['clan_name_2', d.uniprot_domain_2.clan_name],
            ['uniprot_domain_pair_id', d.uniprot_domain_pair_id],]
        append_to_header_df(zip(*header_domain_contact))

    ### Feature information ---------------------------------------------------
    # FoldX output
    if isinstance(d, sql_db.UniprotDomain):
        append_to_feature_df([call_foldx.names_stability_wt, mut.Stability_energy_wt])
        append_to_feature_df([call_foldx.names_stability_mut, mut.Stability_energy_mut])

    elif isinstance(d, sql_db.UniprotDomainPair):
        append_to_feature_df([call_foldx.names_stability_complex_wt, mut.AnalyseComplex_energy_wt])
        append_to_feature_df([call_foldx.names_stability_complex_mut, mut.AnalyseComplex_energy_mut])

    # PhysicoChemical properties
    names_phys_chem = ['pcv_salt_equal', 'pcv_salt_opposite', 'pcv_hbond', 'pcv_vdW']
    append_to_feature_df([[name + '_wt' for name in names_phys_chem], mut.physChem_wt])
    append_to_feature_df([[name + '_self_wt' for name in names_phys_chem], mut.physChem_wt_ownChain])
    append_to_feature_df([[name + '_mut' for name in names_phys_chem], mut.physChem_mut])
    append_to_feature_df([[name + '_self_mut' for name in names_phys_chem], mut.physChem_mut_ownChain])

    # Other features
    other_features_common = [
        ['secondary_structure_wt', mut.secondary_structure_wt],
        ['solvent_accessibility_wt', mut.solvent_accessibility_wt],
        ['secondary_structure_mut', mut.secondary_structure_mut],
        ['solvent_accessibility_mut', mut.solvent_accessibility_mut],
        ['sift_score', mut.provean_score],
        ['normDOPE', m.norm_dope],
        ['matrix_score', mut.matrix_score],]
    append_to_feature_df(zip(*other_features_common))

    if isinstance(d, sql_db.UniprotDomain):
        other_features_domain = [
            ['seq_id_avg', t.alignment_identity],]
        append_to_feature_df(zip(*other_features_domain))

    if isinstance(d, sql_db.UniprotDomainPair):
        other_features_domain_pair = [
            ['seq_id_avg', t.alignment_identity_1 + t.alignment_identity_2],
            ['seq_id_chain1', t.alignment_identity_1],
            ['seq_id_chain2', t.alignment_identity_2],
            ['if_hydrophobic', m.interface_area_hydrophobic],
            ['if_hydrophilic', m.interface_area_hydrophilic],
            ['if_total', m.interface_area_total],
            ['contact_distance_wt', mut.contact_distance_wt],
            ['contact_distance_mut', mut.contact_distance_mut],]
        append_to_feature_df(zip(*other_features_domain_pair))

    ### Join df lists
    header_dfs_joined = header_dfs[0]
    for header_df in header_dfs[1:]:
        header_dfs_joined = header_dfs_joined.join(header_df)
    feature_dfs_joined = feature_dfs[0]
    for feature_df in feature_dfs[1:]:
        feature_dfs_joined = feature_dfs_joined.join(feature_df)
    for name, value in zip(feature_dfs_joined.columns, feature_dfs_joined.loc[0].values):
        if 'secondary_structure' in name:
            feature_dfs_joined.loc[0, name] = secondary_structure_to_int[feature_dfs_joined.loc[0, name]]
    return header_dfs_joined, feature_dfs_joined


def convert_features_to_differences(df, keep_mut=False):
    """
    Return a dataframe with all the `_mut` features replaced with `_change`
    features which coresspond to the difference between `_mut` and `_wt`.
    """
    column_list = []
    for column_name, column in df.iteritems():
        if '_mut' in column_name and column_name.replace('_mut','_wt') in df.columns:
            if keep_mut:
                column_list.append(column)
            new_column = column - df[column_name.replace('_mut','_wt')]
            if 'secondary_structure' in column_name:
                new_column = new_column.apply(lambda x: 1 if x else 0)
            new_column.name = column_name.replace('_mut','_change')
            column_list.append(new_column)
        else:
            column_list.append(column)
    new_df = pd.DataFrame(column_list).T
    return new_df


def encode_list_as_text(list_of_lists):
    return ','.join([':'.join([str(x) for x in xx]) for xx in zip(*list_of_lists)])


def decode_text_as_list(list_string):
    str2num = lambda x: float(x) if '.' in x else int(x)
    return zip(*[[str2num(x) for x in sublist.split(':')] for sublist in list_string.split(',')])


###############################################################################
class MutationData(object):

    def __init__(self):
        self.d = None
        self.t = None
        self.m = None
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
            db, log, n_cores, bin_path, foldX_WATER, build_model_runs, matrix, gap_s, gap_e):
        """
        """
        self.global_temp_path = global_temp_path
        self.temp_path = temp_path
        self.unique = unique
        self.unique_temp_folder = temp_path + unique + '/'
        self.pdb_path = pdb_path
        self.db = db # sql database
        self.log = log
        self.bin_path = bin_path
        self.foldX_WATER = foldX_WATER
        self.build_model_runs = build_model_runs
        self.matrix = matrix
        self.gap_s = gap_s
        self.gap_e = gap_e
        self.n_cores = n_cores
        with open(bin_path + 'clf_domain.pickle', 'rb') as ifh:
            self.clf_domain, self.clf_domain_features = pickle.load(ifh)
        with open(bin_path + 'clf_interface.pickle', 'rb') as ifh:
            self.clf_interface, self.clf_interface_features = pickle.load(ifh)


    def get_mutation_data(self, d, t, m, uniprot_id_1, mutation):
        """
        """
        path_to_provean_supset = None
        #######################################################################
        # Core
        if isinstance(d, sql_db.UniprotDomain):
            self.logger.debug("Analyzing core mutation for uniprot: %s" % uniprot_id_1)
            pfam_name = d.pfam_name
            domain_start, domain_end = sql_db.decode_domain(t.domain_def)
            alignment, __ = self.db.get_alignment(t, d.path_to_data)
            alignment_id = t.alignment_id
            chains_pdb = [t.domain.pdb_chain, ]
            uniprot_sequences = [self.db.get_uniprot_sequence(uniprot_id_1), ]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],]
            chains_modeller = [m.chain, ]
            uniprot_domain_id = d.uniprot_domain_id
            if t.provean_supset_filename:
                path_to_provean_supset = self.temp_path + d.path_to_data + t.provean_supset_filename

        #######################################################################
        # Interface
        elif isinstance(d, sql_db.UniprotDomainPair):
            if uniprot_id_1 == d.uniprot_domain_1.uniprot_id:
                uniprot_id_2 = d.uniprot_domain_2.uniprot_id
                d_1, d_2 = d.uniprot_domain_1, d.uniprot_domain_2
                domain_start, domain_end = sql_db.decode_domain(t.domain_def_1)
                domain_2_start, domain_2_end = sql_db.decode_domain(t.domain_def_2)
                alignment, __ = self.db.get_alignment(t, d.path_to_data)
                alignment_id = t.alignment_id_1
                chains_pdb = [t.domain_1.pdb_chain, t.domain_2.pdb_chain]
                interacting_aa = [int(uniprot_num) for uniprot_num in m.interacting_aa_1.split(',') if uniprot_num]
                chains_modeller = [m.chain_1, m.chain_2]

            elif uniprot_id_1 == d.uniprot_domain_2.uniprot_id:
                self.logger.debug('Mutated uniprot is uniprot 2. Rearranging...')
                uniprot_id_2 = d.uniprot_domain_1.uniprot_id
                d_1, d_2 = d.uniprot_domain_2, d.uniprot_domain_1
                domain_start, domain_end = sql_db.decode_domain(t.domain_def_2)
                domain_2_start, domain_2_end = sql_db.decode_domain(t.domain_def_1)
                __, alignment = self.db.get_alignment(t, d.path_to_data)
                alignment_id = t.alignment_id_2
                chains_pdb = [t.domain_2.pdb_chain, t.domain_1.pdb_chain]
                interacting_aa = [int(uniprot_num) for uniprot_num in m.interacting_aa_2.split(',') if uniprot_num]
                chains_modeller = [m.chain_2, m.chain_1]

            self.logger.debug('Analysing interface mutation between uniprots %s and %s' % (uniprot_id_1, uniprot_id_2,) )
            pfam_name = d_1.pfam_name
            uniprot_sequences = [self.db.get_uniprot_sequence(d_1.uniprot_id),
                                 self.db.get_uniprot_sequence(d_2.uniprot_id)]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],
                                uniprot_sequences[1][domain_2_start-1:domain_2_end]]

            # HETflag is a dictionary, doesn't matter the order. But key-value should match
            uniprot_domain_id = d_1.uniprot_domain_id

            if int(mutation[1:-1]) not in interacting_aa:
                raise errors.MutationOutsideInterface('mutated residue not involved in the interaction')
            if d_1.template != None and d_1.template.provean_supset_filename != None:
                # This will only work if the core mutations are evaluated at the same time as interface mutations
                if os.path.isfile(self.temp_path + d_1.path_to_data + d_1.template.provean_supset_filename):
                    path_to_provean_supset = self.temp_path + d_1.path_to_data + d_1.template.provean_supset_filename
                else:
                    self.logger.error('Provean supporting set does not exist even though it should!!!')

        #######################################################################
        # Common
        if int(mutation[1:-1]) < domain_start or int(mutation[1:-1]) > domain_end:
#            raise error.MutationOutsideDomain('Mutation %s outside of domain %s with boundaries %i:%i' \
#            % (mutation, pfam_name, domain_start, domain_end))
            raise errors.MutationOutsideDomain('Mutation falls outside domain')

        save_path = self.temp_path + d.path_to_data
        pdbFile_wt = m.model_filename

        # Convert mutation from uniprot domain numbering to pdb domain numbering
        position_domain = int(mutation[1:-1]) - domain_start + 1 # +1 to have the numbering staring with 1
        mutation_domain = mutation[0] + str(position_domain) + mutation[-1]
#        mutation_structure = map_mutation_to_structure(alignment, alignment_id, mutation_domain, uniprot_id_1)
#        chain_mutation_structure = chains_pdb[0] + '_' + mutation_structure

        self.logger.debug(save_path + pdbFile_wt)
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        structure = parser.get_structure('ID', save_path + pdbFile_wt)
        position_modeller = pdb_template.convert_position_to_resid(structure[0], chains_modeller[0], [position_domain])
        mutation_modeller = mutation[0] + position_modeller[0] + mutation[-1]

        # Save the results
        mut_data = MutationData()

        mut_data.d = d
        mut_data.t = t
        mut_data.m = m
        mut_data.uniprot_id_1 = uniprot_id_1
        mut_data.mutation = mutation
        mut_data.pfam_name = pfam_name
        mut_data.domain_start = domain_start
        mut_data.domain_end = domain_end
        mut_data.alignment = alignment
        mut_data.alignment_id = alignment_id
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

        return mut_data


    def evaluate_mutation(self, mut_data, uniprot_mutation):
        """
        """
        if (mut_data.path_to_provean_supset and not uniprot_mutation.provean_score):
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

        if not uniprot_mutation.Stability_energy_wt:
            self.logger.debug('Evaluating the structural impact of the mutation...')
            uniprot_mutation = self.evaluate_structural_impact(mut_data, uniprot_mutation)

        if (not uniprot_mutation.ddG and
                (uniprot_mutation.provean_score and
                uniprot_mutation.Stability_energy_wt)):
            self.logger.debug('Predicting the thermodynamic effect of the mutation...')
            uniprot_mutation = self.predict_thermodynamic_effect(mut_data.d, mut_data.t, mut_data.m, uniprot_mutation)

        return uniprot_mutation


    def get_provean_score(self, uniprot_domain_id, domain_mutation, domain_sequence, path_to_provean_supset):
        """
        """
#        uniprot_sequence = self.db.get_uniprot_sequence(d.uniprot_id).seq.tostring()
#        domain_def = sql_db.decode_domain(t.domain_def)
#        uniprot_sequence_domain = uniprot_sequence[domain_def[0]-1:domain_def[1]]
        if isinstance(domain_sequence, str):
            pass
        elif isinstance(domain_sequence, Seq):
            domain_sequence = str(domain_sequence)
        elif isinstance(domain_sequence, SeqRecord):
            domain_sequence = str(domain_sequence.seq)
        else:
            raise Exception('Wrong class type %s for domain_sequence' % str(type(domain_sequence)))
        domain_sequence_seqrec = SeqRecord(seq=Seq(domain_sequence), id=str(uniprot_domain_id))
        SeqIO.write(domain_sequence_seqrec, self.unique_temp_folder + 'sequence_conservation/sequence.fasta', 'fasta')

#        first_aa = uniprot_sequence_domain[0]
    #    provean_supset_filename = t.alignment_filename.replace('.aln', '.supset')
#        provean_supset_filename = t.provean_supset_filename

        subprocess.check_call('echo ' + domain_mutation + ' > ' + self.unique_temp_folder + 'sequence_conservation/decoy.var', shell=True)
        path_to_provean_supset_local = self.unique_temp_folder + 'sequence_conservation/' + path_to_provean_supset.split('/')[-1]
        subprocess.check_call('cp -f ' + path_to_provean_supset + ' ' + path_to_provean_supset_local, shell=True)
        system_command = (
            './provean' +
            ' -q ' + './sequence.fasta' +
            ' -v ' + './decoy.var' +
            ' -d ' + self.global_temp_path + 'blast/db/nr' +
            ' --tmp_dir ' + self.unique_temp_folder +
            ' --num_threads ' + '{}'.format(self.n_cores) +
            ' --psiblast ' + hf.get_which('psiblast') +
            ' --blastdbcmd ' + hf.get_which('blastdbcmd') +
            ' --cdhit ' + hf.get_which('cd-hit') +
            ' --supporting_set ' + path_to_provean_supset_local)
        child_process = hf.run_subprocess_locally(
            self.unique_temp_folder + 'sequence_conservation/',
            system_command)
        result, error_message = child_process.communicate()
        self.logger.debug(result)
        while (child_process.returncode != 0
                and 'IDs are not matched' in error_message):
            self.logger.error(error_message)
            line_to_remove = error_message.split(':')[1].split(',')[0]
            self.logger.error('Removing line with id: {} from the supporting set...'.format())
            with open(path_to_provean_supset_local, 'r') as ifh, \
                    open(path_to_provean_supset_local + '.mod', 'w') as ofh:
                for line in ifh:
                    if line_to_remove not in line:
                        ofh.write(line)
            subprocess.check_call('mv -f ' + path_to_provean_supset_local + '.mod ' + path_to_provean_supset_local, shell=True)
            child_process = hf.run_subprocess_locally(
                self.unique_temp_folder + 'sequence_conservation/',
                system_command)
            result, error_message = child_process.communicate()
            self.logger.debug(result)

        if child_process.returncode != 0:
            self.logger.error(error_message)
            self.logger.error(system_command)
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


    def evaluate_structural_impact(self, mut_data, uniprot_mutation):
        d = mut_data.d
#        t = mut_data.t
#        model = mut_data.structure[0]
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
                   self.log)
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
        assert(str(mut_data.domain_sequences[0].seq)[int(mut_data.mutation_domain[1:-1])-1] == mut_data.mutation_domain[0])
        mutCodes = [mutation[0] + mut_data.chains_modeller[0] + mut_data.position_modeller[0] + mutation[-1], ]
        self.logger.debug('Mutcodes for foldx:')
        self.logger.debug(mutCodes)

        # Introduce the mutation using foldX
        fX_wt = call_foldx.FoldX(self.unique_temp_folder,
                      repairedPDB_wt,
                      mut_data.chains_modeller[0],
                      self.build_model_runs,
                      self.foldX_WATER,
                      self.log)
        repairedPDB_wt_list, repairedPDB_mut_list = fX_wt('BuildModel', mutCodes)

        wt_chain_sequences = pdb_template.get_chain_sequences(repairedPDB_wt_list[0])
        mut_chain_sequences = pdb_template.get_chain_sequences(repairedPDB_mut_list[0])
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
                                    self.log))

        fX_mut_list = list()
        for mPDB in repairedPDB_mut_list:
            fX_mut_list.append(call_foldx.FoldX(self.unique_temp_folder,
                                     mPDB,
                                     mut_data.chains_modeller[0],
                                     self.build_model_runs,
                                     self.foldX_WATER,
                                     self.log))

        #######################################################################
        ## 5th: calculate the energy for the wildtype
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

        physi_chem = analyze_structure.PhysiChem(5.0, 4.0, self.unique_temp_folder, self.log)
        opposite_chain_contact_vector_all_wt, same_chain_contact_vector_all_wt = get_contact_vectors(
            physi_chem, repairedPDB_wt_list, mut_data.chains_modeller[0], mut_data.mutation_domain)
        opposite_chain_contact_vector_all_mut, same_chain_contact_vector_all_mut = get_contact_vectors(
            physi_chem, repairedPDB_wt_list, mut_data.chains_modeller[0], mut_data.mutation_domain)

        #######################################################################
        # Calculate secondary structure, sasa, and interchain distance
        analyze_structure_wt = analyze_structure.AnalyzeStructure(
            self.unique_temp_folder + 'FoldX/',
            self.unique_temp_folder + 'analyze_structure/',
            repairedPDB_wt_list[0].split('/')[-1], # dssp file wildtype
            mut_data.chains_modeller, None, self.log)
        (seasa_by_chain_together, seasa_by_chain_separately,
        seasa_by_residue_together, seasa_by_residue_separately) = analyze_structure_wt.get_seasa()
        seasa_info_wt = seasa_by_residue_separately[
            (seasa_by_residue_separately['pdb_chain']==mut_data.chains_modeller[0]) &
            (seasa_by_residue_separately['res_num']==mut_data.position_modeller[0])].iloc[0]
        if pdb_template.convert_aa(seasa_info_wt['res_name']) != mut_data.mutation_domain[0]:
            self.logger.error('Wrong amino acid for msms wild-type!')
            self.logger.error(seasa_info_wt)
            self.logger.error(seasa_by_residue_separately)
            raise Exception('surface area calculated for the wrong atom!')
        solvent_accessibility_wt = seasa_info_wt['rel_sasa']

        secondary_structure_wt, solvent_accessibility_dssp_wt = analyze_structure_wt.get_dssp()
        secondary_structure_wt = secondary_structure_wt[mut_data.chains_modeller[0]][int(mut_data.mutation_domain[1:-1])-1]

        if isinstance(d, sql_db.UniprotDomainPair):
            contact_distance_wt = analyze_structure_wt.get_interchain_distances(mut_data.chains_modeller[0], mut_data.mutation_modeller)
            try:
                contact_distance_wt = contact_distance_wt[mut_data.chains_modeller[0]][0]
            except IndexError:
                self.logger.error(contact_distance_wt)
                raise

        analyze_structure_mut = analyze_structure.AnalyzeStructure(
            self.unique_temp_folder + 'FoldX/',
            self.unique_temp_folder + 'analyze_structure/',
            repairedPDB_mut_list[0].split('/')[-1], # dssp file mut
            mut_data.chains_modeller, None, self.log)
        (seasa_by_chain_together, seasa_by_chain_separately,
        seasa_by_residue_together, seasa_by_residue_separately) = analyze_structure_mut.get_seasa()
        seasa_info_mut = seasa_by_residue_separately[
            (seasa_by_residue_separately['pdb_chain']==mut_data.chains_modeller[0]) &
            (seasa_by_residue_separately['res_num']==mut_data.position_modeller[0])].iloc[0]
        if pdb_template.convert_aa(seasa_info_mut['res_name']) != mut_data.mutation_domain[-1]:
            self.logger.error('Wrong amino acid for msms mutant!')
            self.logger.error(pdb_template.convert_aa(seasa_info_mut['res_name']))
            self.logger.error(mut_data.mutation_domain[-1])
            self.logger.error(seasa_info_mut)
            self.logger.error(seasa_by_residue_separately)
            raise Exception('surface area calculated for the wrong atom!')
        solvent_accessibility_mut = seasa_info_mut['rel_sasa']

        secondary_structure_mut, solvent_accessibility_dssp_mut = analyze_structure_mut.get_dssp()
        secondary_structure_mut = secondary_structure_mut[mut_data.chains_modeller[0]][int(mut_data.mutation_domain[1:-1])-1]

        if isinstance(d, sql_db.UniprotDomainPair):
            contact_distance_mut = analyze_structure_wt.get_interchain_distances(mut_data.chains_modeller[0], mut_data.mutation_modeller)
            try:
                contact_distance_mut = contact_distance_mut[mut_data.chains_modeller[0]][0]
            except IndexError:
                self.logger.error(contact_distance_mut)
                raise

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

        uniprot_mutation.Stability_energy_wt = stability_values_wt
        uniprot_mutation.Stability_energy_mut = stability_values_mut

        uniprot_mutation.physChem_wt = opposite_chain_contact_vector_all_wt
        uniprot_mutation.physChem_wt_ownChain = same_chain_contact_vector_all_wt
        uniprot_mutation.physChem_mut = opposite_chain_contact_vector_all_mut
        uniprot_mutation.physChem_mut_ownChain = same_chain_contact_vector_all_mut

        uniprot_mutation.matrix_score = matrix_score

        uniprot_mutation.secondary_structure_wt = secondary_structure_wt
        uniprot_mutation.solvent_accessibility_wt = solvent_accessibility_wt
        uniprot_mutation.secondary_structure_mut = secondary_structure_mut
        uniprot_mutation.solvent_accessibility_mut = solvent_accessibility_mut

        if isinstance(d, sql_db.UniprotDomainPair):
            uniprot_mutation.AnalyseComplex_energy_wt = complex_stability_values_wt
            uniprot_mutation.AnalyseComplex_energy_mut = complex_stability_values_mut
            uniprot_mutation.contact_distance_wt = contact_distance_wt
            uniprot_mutation.contact_distance_mut = contact_distance_mut
        #######################################################################
        # Save alignments and modeller models to output database for storage
        # Move template files to the output folder as a tar archives
#        self.make_tarfile(self.HOME + self.outputPath + save_path + '_' + mutation + '.tar.bz2',
#                          save_path[:-1])

        #######################################################################
        self.logger.info('Finished processing template:')
        self.logger.info(mut_data.save_path.split('/')[-2])

        return uniprot_mutation


    def __check_structure_match(self, pdb_filename, mut_data, expecte_aa):
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


    def predict_thermodynamic_effect(self, d, t, m, mut):
        """
        """
        if isinstance(d, sql_db.UniprotDomain):
            clf = self.clf_domain
            clf_features = self.clf_domain_features
        elif isinstance(d, sql_db.UniprotDomainPair):
            clf = self.clf_interface
            clf_features = self.clf_interface_features
        header_df, feature_df = get_mutation_feature_vector(d, t, m, mut)
        feature_df = convert_features_to_differences(feature_df, True) # keep mut, remove it in next step
        for column_name in feature_df.columns:
            if column_name not in clf_features:
                feature_df.drop(column_name, axis=1, inplace=True)
        mut.ddG = clf.predict(feature_df)[0]
        return mut


###############################################################################
# Methods when the input is a raw crystal structure

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
        domains = [['Null', 'Null'], ] # extract the full sequence

        pdb = PDBTemplate(self.pdb_path, pdbCode, chain, domains, self.unique_temp_folder, self.unique_temp_folder, self.log)

        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()
        chainNumberingDomain = pdb.get_chain_numbering(chain)
        if chainNumberingDomain == []:
            raise errors.PDBError('Could not get the pdb numbering for ' + pdbCode + '_' + chain)

        return next(SeqIO.parse(self.unique_temp_folder + '/' + pdbCode + chain + '.seq.txt', 'fasta')), chainNumberingDomain



    def _getCrystalStructure(self, pdbCode, chains, savePDB, FULL_PATH=False):
        domains = [['Null', 'Null'] for i in range(len(chains)) ]
        pdb = PDBTemplate(self.pdb_path, pdbCode, chains, domains, savePDB, self.unique_temp_folder, self.log)
        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()

        SWITCH_CHAIN = False
        if chains != chains_pdb_order:
            SWITCH_CHAIN = True
            chains = chains_pdb_order

        return '-', pdbCode.split('.')[0] + ''.join(chains) + '.pdb', SWITCH_CHAIN, chains



