# -*- coding: utf-8 -*-

import os
import subprocess
import tarfile

import errors
import sql_db
import analyze_structure

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser

import pdb_template
from pdb_template import PDBTemplate
from physicochemical_properties import pysiChem
from call_foldx import foldX
import helper_functions as hf

###############################################################################
# Useful functions

def map_mutation_to_structure(alignment, alignment_id, mutation_domain):
    """
    """

    if alignment[0].id == alignment_id:
        alignment_uniprot = str(alignment[1].seq)
        alignment_protein = str(alignment[0].seq)
    elif alignment[1].id == alignment_id:
        alignment_uniprot = str(alignment[0].seq)
        alignment_protein = str(alignment[1].seq)

    # now get the position
    uniprot_position = 0
    pdb_position = 0
    for idx in range(len(alignment_uniprot)):
        if uniprot_position >= int(mutation_domain[1:-1])-1 and not alignment_uniprot[idx] == '-':
            break
        if alignment_uniprot[idx] != '-':
            uniprot_position += 1
        if alignment_protein[idx] != '-':
            pdb_position += 1

    assert alignment_uniprot[idx] == mutation_domain[0]
    mutation_structure = mutation_domain[0] + str(pdb_position + 1) + mutation_domain[-1]
    return mutation_structure



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



def prepareMutationFoldX(sequence, mutation):
    """
    Create the sequence snippets for the mutation with FoldX
    also, record the position of the mutation in the snippet
    """

    if int(mutation[1:-1]) < 7:
        left = 0
        m_pos = int(mutation[1:-1]) # the position of the muation in the snippet
    else:
        left = int(mutation[1:-1])-7
        m_pos = 7
    if int(mutation[1:-1]) >= len(sequence) - 8:
        right = len(sequence)
    else:
        right = int(mutation[1:-1])+8

    wt_seq = sequence.seq[left:right]
    if sequence.seq.count(wt_seq) > 1:
        m_pos = None

    mut_seq = sequence.seq[left:int(mutation[1:-1])-1] + mutation[-1] + sequence.seq[int(mutation[1:-1]):right]

    mutation_foldX = (m_pos, str(wt_seq) + '\n' + str(mut_seq))
    return mutation_foldX



def setChainID(chains, mutation, HETflag, SWITCH_CHAIN, do_modelling):
    """
    modeller renames the chains, in order to keep track of the chain IDs
    they are renamed here to match the modeller output. Modeller labels the
    chains subsequently, i.e. A, B, C etc., thus, depending on wether HETATMs
    are present as seperate chain, the chain labels are A, B or A, C
    """

    mutChain = str(chains.index(mutation[0]))

    # get the chain names correct
    if do_modelling:
        if len(chains) == 1:
            # if no HETATMs are present, the chain ID is empty
            # otherwise the chain IDs are A and B for the atoms
            # and HETATMs respectively
            if HETflag[chains[0]] == False:
                chains_new = [' ']
                mutChain = ' '
            else:
                chains_new = ['A']
                mutChain = 'A'
        else:
            a = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            index = 0
            chains_new = list()
            for chainID in chains:
                if HETflag[chainID] == True:
                    # if True, then there are HETATMs associated to this chain
                    # in modeller they will be listed as two different chains.
                    # thus the chain names have to jump by two
                    chains_new.append(a[index])
                    index += 2
                else:
                    chains_new.append(a[index])
                    index += 1

            if mutChain == '0':
                mutChain = 'A'
            elif mutChain == '1':
                if chains_new[1] == 'B':
                    mutChain = 'B'
                elif chains_new[1] == 'C':
                    mutChain = 'C'
                else:
                    mutChain = None

    return mutChain + '_' + mutation[2:]


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:bz2") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


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
            db, log, n_cores, foldX_WATER, build_model_runs, matrix, gap_s, gap_e):
        """
        """
        self.global_temp_path = global_temp_path
        self.temp_path = temp_path
        self.unique = unique
        self.unique_temp_folder = temp_path + unique + '/'
        self.pdb_path = pdb_path
        self.db = db # sql database
        self.log = log
        self.foldX_WATER = foldX_WATER
        self.build_model_runs = build_model_runs
        self.matrix = matrix
        self.gap_s = gap_s
        self.gap_e = gap_e
        self.n_cores = n_cores


    def get_mutation_data(self, d, t, m, uniprot_id_1, mutation):
        """
        """
        path_to_provean_supset = None
        #######################################################################
        # Core
        if isinstance(d, sql_db.UniprotDomain):
            self.log.debug("Analyzing core mutation for uniprot: %s" % uniprot_id_1)
            pfam_name = d.pfam_name
            domain_start, domain_end = sql_db.decode_domain(t.domain_def)
            alignment, __ = self.db.get_alignment(t, d.path_to_data)
            alignment_id = t.alignment_id
            chains_pdb = [t.domain.pdb_chain, ]
            uniprot_sequences = [self.db.get_uniprot_sequence(uniprot_id_1), ]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],]
            HETflag = {chains_pdb[0]: m.het_flag}
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
                self.log.debug('Mutated uniprot is uniprot 2. Rearranging...')
                uniprot_id_2 = d.uniprot_domain_1.uniprot_id
                d_1, d_2 = d.uniprot_domain_2, d.uniprot_domain_1
                domain_start, domain_end = sql_db.decode_domain(t.domain_def_2)
                domain_2_start, domain_2_end = sql_db.decode_domain(t.domain_def_1)
                __, alignment = self.db.get_alignment(t, d.path_to_data)
                alignment_id = t.alignment_id_2
                chains_pdb = [t.domain_2.pdb_chain, t.domain_1.pdb_chain]
                interacting_aa = [int(uniprot_num) for uniprot_num in m.interacting_aa_2.split(',') if uniprot_num]
                chains_modeller = [m.chain_2, m.chain_1]

            self.log.debug('Analysing interface mutation between uniprots %s and %s' % (uniprot_id_1, uniprot_id_2,) )
            pfam_name = d_1.pfam_name
            uniprot_sequences = [self.db.get_uniprot_sequence(d_1.uniprot_id),
                                 self.db.get_uniprot_sequence(d_2.uniprot_id)]
            domain_sequences = [uniprot_sequences[0][domain_start-1:domain_end],
                                uniprot_sequences[1][domain_2_start-1:domain_2_end]]
            # Commented out because sometimes the core domains have not been calculated and
            # the lines below give an error.
#            domain_core_start, domain_core_end = sql_db.decode_domain(d_1.template[0].domain_def)
#            position_core_domain = int(mutation[1:-1]) - domain_core_start + 1 # +1 to have the numbering staring with 1

            # HETflag is a dictionary, doesn't matter the order. But key-value should match
            HETflag = {t.domain_1.pdb_chain: m.het_flag_1,
                       t.domain_2.pdb_chain: m.het_flag_2}
            uniprot_domain_id = d_1.uniprot_domain_id

            if int(mutation[1:-1]) not in interacting_aa:
                raise errors.MutationOutsideInterface('mutated residue not involved in the interaction')
            if d_1.template != None and d_1.template.provean_supset_filename != None:
                # This will only work if the core mutations are evaluated at the same time as interface mutations
                if os.path.isfile(self.temp_path + d_1.path_to_data + d_1.template.provean_supset_filename):
                    path_to_provean_supset = self.temp_path + d_1.path_to_data + d_1.template.provean_supset_filename
                else:
                    self.log.error('Provean supporting set does not exist even though it should!!!')

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

        self.log.debug(save_path + pdbFile_wt)
        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
        structure = parser.get_structure('ID', save_path + pdbFile_wt)
        position_modeller = pdb_template.convert_position_to_resid(structure[0], chains_modeller[0], [position_domain])
        mutation_modeller = mutation[0] + position_modeller[0] + mutation[-1]

        # Save the results
        mut_data = MutationData()

        mut_data.d = d
        mut_data.t = t
#        mut_data.m = m
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
        mut_data.HETflag = HETflag
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
            self.log.debug(key + ':')
            self.log.debug(value)

        return mut_data


    def evaluate_mutation(self, mut_data):

        parser = PDBParser(QUIET=True) # set QUIET to False to output warnings like incomplete chains etc.
#        model = mut_data.structure[0]
        d, t = mut_data.d, mut_data.t
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
            self.log.error('cp result: {}'.format(result))
            self.log.error('cp error: {}'.format(e))


        #######################################################################
        ## 0th: get Provean score for the mutation
        provean_mutation, provean_score = None, None
        if mut_data.path_to_provean_supset:
            try:
                provean_mutation, provean_score = \
                self.get_provean_score(
                    mut_data.uniprot_domain_id,
                    mut_data.mutation_domain,
                    mut_data.domain_sequences[0],
                    mut_data.path_to_provean_supset)
            except errors.ProveanError as e:
                self.log.error(str(type(e)))
        self.log.debug('provean mutation:')
        self.log.debug(provean_mutation)
        self.log.debug('provean score:')
        self.log.debug(provean_score)


        #######################################################################
        ## 1st: get a snippet of AAs around the mutation, as well as the mutation residue (no chain)
        # This is required input for FoldX
        if len(mut_data.domain_sequences) == 1:
            self.log.debug(mut_data.domain_sequences[0].seq)
        else:
            self.log.debug(mut_data.domain_sequences[0].seq)
            self.log.debug(mut_data.domain_sequences[1].seq)
#        if not switch_chain:
#            mutations_foldX = [prepareMutationFoldX(domain_sequences[0], mutation_domain),]
#        else:
#            mutations_foldX = [prepareMutationFoldX(domain_sequences[1], mutation_domain),]

        assert(str(mut_data.domain_sequences[0].seq)[int(mut_data.mutation_domain[1:-1])-1] == mut_data.mutation_domain[0])
        mutations_foldX = [prepareMutationFoldX(mut_data.domain_sequences[0], mut_data.mutation_domain),]
        self.log.debug("mutations_foldX:")
        self.log.debug(mutations_foldX)


        #######################################################################
        ## 2nd: use the 'Repair' feature of FoldX to optimise the structure
        fX = foldX(self.unique_temp_folder,
                   mut_data.save_path + mut_data.pdbFile_wt,
                   mut_data.chains_modeller[0],
                   self.build_model_runs,
                   self.foldX_WATER,
                   self.log)
        repairedPDB_wt = fX.run('RepairPDB')


        #######################################################################
        ## 3rd: introduce the mutation using FoldX

        # compile a list of mutations
        mutCodes = list()
        for mut in mutations_foldX:
            if mut[0]:
                mutCodes.append(mut[1])
            else:
                mutCodes.append(mutation[0] + mut_data.chains_modeller[0] + mut_data.position_modeller[0] + mutation[-1])

        mutCodes = [mutation[0] + mut_data.chains_modeller[0] + mut_data.position_modeller[0] + mutation[-1], ]
        self.log.debug('Mutcodes for foldx:')
        self.log.debug(mutCodes)

        fX_wt = foldX(self.unique_temp_folder,
                      repairedPDB_wt,
                      mut_data.chains_modeller[0],
                      self.build_model_runs,
                      self.foldX_WATER,
                      self.log)

        # do the mutation with foldX
        repairedPDB_wt_list, repairedPDB_mut_list = fX_wt.run('BuildModel', mutCodes)

        wt_chain_sequences = pdb_template.get_chain_sequences(repairedPDB_wt_list[0])
        mut_chain_sequences = pdb_template.get_chain_sequences(repairedPDB_mut_list[0])
        self.log.debug('repairedPDB_wt_list: %s' % str(repairedPDB_wt_list))
        self.log.debug('wt_chain_sequences: %s' % str(wt_chain_sequences))
        self.log.debug('repairedPDB_mut_list: %s' % str(repairedPDB_mut_list))
        self.log.debug('mut_chain_sequences: %s' % str(mut_chain_sequences))

        def check_structure_match(pdb_filename, expecte_aa):
            structure = parser.get_structure('ID', pdb_filename)
            model = structure[0]
            chain = model[mut_data.chains_modeller[0]]
            residue_found = False
            for residue in chain:
                if (str(residue.id[1]) + residue.id[2].strip()) == mut_data.position_modeller[0]:
                    residue_found = True
                    break
            if not residue_found or not (pdb_template.convert_aa(residue.resname) == expecte_aa):
                self.log.error('residue_found? %s' % residue_found)
                self.log.error(residue.resname)
                self.log.error(residue.id)
                self.log.error(mut_data.position_modeller)
                self.log.error(expecte_aa)
                raise Exception('Expected and actual FoldX amino acids do not match!')

        check_structure_match(repairedPDB_wt_list[0], mutation[0])
        check_structure_match(repairedPDB_mut_list[0], mutation[-1])

        # Copy the foldX wildtype pdb file (use the first model if there are multiple)
        model_filename_wt = uniprot_id_1 + '_' + mutation + '/' +  repairedPDB_wt_list[0].split('/')[-1]
        model_filename_mut = uniprot_id_1 + '_' + mutation + '/MUT_' +  repairedPDB_mut_list[0].split('/')[-1]
        subprocess.check_call('mkdir -p ' + mut_data.save_path + uniprot_id_1 + '_' + mutation + '/', shell=True)
        subprocess.check_call('cp ' + repairedPDB_wt_list[0] + ' ' + mut_data.save_path + model_filename_wt, shell=True)
        subprocess.check_call('cp ' + repairedPDB_mut_list[0] + ' ' + mut_data.save_path + model_filename_mut, shell=True)


        #######################################################################
        ## 4th: set up the classes for the wildtype and the mutant structures
        fX_wt_list = list()
        for wPDB in repairedPDB_wt_list:
            fX_wt_list.append(foldX(self.unique_temp_folder,
                                    wPDB,
                                    mut_data.chains_modeller[0],
                                    self.build_model_runs,
                                    self.foldX_WATER,
                                    self.log))

        fX_mut_list = list()
        for mPDB in repairedPDB_mut_list:
            fX_mut_list.append(foldX(self.unique_temp_folder,
                                     mPDB,
                                     mut_data.chains_modeller[0],
                                     self.build_model_runs,
                                     self.foldX_WATER,
                                     self.log))

        #######################################################################
        ## 5th: calculate the energy for the wildtype
        FoldX_StabilityEnergy_wt = list()
        FoldX_AnalyseComplex_wt = list()
        for fX_wt in fX_wt_list:
            # stability of the interaction
            if len([ item for item in mut_data.chains_modeller if item != '' ]) == 1:
                AnalyseComplex = [ '0' for i in range(25) ]
            else:
                AnalyseComplex = fX_wt.run('AnalyseComplex')
                AnalyseComplex = [ item for item in AnalyseComplex ] # convert to list

            # stability
            StabilityEnergy = fX_wt.run('Stability')
            StabilityEnergy = [ item for item in StabilityEnergy ] # convert to list

            # append the results
            FoldX_StabilityEnergy_wt.append(StabilityEnergy)
            FoldX_AnalyseComplex_wt.append(AnalyseComplex)


        #######################################################################
        ## 6th: calculate the energy for the mutant
        FoldX_StabilityEnergy_mut = list()
        FoldX_AnalyseComplex_mut = list()
        for fX_mut in fX_mut_list:
            # stability of the interaction
            if len([ item for item in mut_data.chains_modeller if item != '' ]) == 1:
                AnalyseComplex = [ '0' for i in range(25) ]
            else:
                AnalyseComplex = fX_mut.run('AnalyseComplex')
                AnalyseComplex = [ item for item in AnalyseComplex ] # convert to list

            # stability
            StabilityEnergy = fX_mut.run('Stability')
            StabilityEnergy = [ item for item in StabilityEnergy ] # convert to list

            # append the results
            FoldX_StabilityEnergy_mut.append(StabilityEnergy)
            FoldX_AnalyseComplex_mut.append(AnalyseComplex)


        #######################################################################
        ## 7th: combine the energy results
        if isinstance(d, sql_db.UniprotDomainPair):
            # AnalyseComplex
            if len(FoldX_AnalyseComplex_wt) == 1:
                AnalyseComplex_energy_wt = ','.join(FoldX_AnalyseComplex_wt[0])
                AnalyseComplex_energy_mut = ','.join(FoldX_AnalyseComplex_mut[0])
            else:
                # zip the output for the wildtype
                AnalyseComplex_energy_wt = [ [] for i in range(len(FoldX_AnalyseComplex_wt[0])) ]
                for item in FoldX_AnalyseComplex_wt:
                    i = 0
                    for element in item:
                        AnalyseComplex_energy_wt[i].append(element)
                        i += 1
                AnalyseComplex_energy_wt = ':'.join([','.join(item) for item in AnalyseComplex_energy_wt])

                # zip the output for the mutant
                AnalyseComplex_energy_mut = [ [] for idx in range(len(FoldX_AnalyseComplex_mut[0])) ]
                for item in FoldX_AnalyseComplex_mut:
                    i = 0
                    for element in item:
                        AnalyseComplex_energy_mut[i].append(element)
                        i += 1
                AnalyseComplex_energy_mut = ':'.join([','.join(item) for item in AnalyseComplex_energy_mut])

        # Stability
        if len(FoldX_StabilityEnergy_wt) == 1:
            Stability_energy_wt = ','.join(FoldX_StabilityEnergy_wt[0])
            Stability_energy_mut = ','.join(FoldX_StabilityEnergy_mut[0])
        else:
            # zip the output for the wildtype
            Stability_energy_wt = [ [] for idx in range(len(FoldX_StabilityEnergy_wt[0])) ]
            for item in FoldX_StabilityEnergy_wt:
                i = 0
                for element in item:
                    Stability_energy_wt[i].append(element)
                    i += 1
            Stability_energy_wt = ':'.join([','.join(item) for item in Stability_energy_wt])

            # zip the output for the mutant
            Stability_energy_mut = [ [] for idx in range(len(FoldX_StabilityEnergy_mut[0])) ]
            for item in FoldX_StabilityEnergy_mut:
                i = 0
                for element in item:
                    Stability_energy_mut[i].append(element)
                    i += 1
            Stability_energy_mut = ':'.join([','.join(item) for item in Stability_energy_mut])


#        #######################################################################
#        ## 8th: calculate the interface size
#        ## don't do it if no complex is modelled
#        if len([ item for item in chains if item != '' ]) == 1:
#            interface_size = ['0', '0', '0']
#        else:
#            try:
#                pops = pdb_template.GetSASA('', self.unique_temp_folder + 'pops/')
#                # here the interface between chain 0 and all the other chains is
#                # calculated. If more than two chains are present, you can change
#                # this behaviour here by adding the addintional chain IDs to the list
#
#    #                chains_complex = [ chain for chain in chains[0] ]
#                interface_size = pops(repairedPDB_wt_list[0], ['A', ])
#            except:
#                print "pops did not work!"
#                self.log.error('pops did not work!')
#                interface_size = ['0', '0', '0']


        #######################################################################
        ## 9th: calculate the pysico-chemical properties
        ## copy the pdb files to the result folder
        get_atomicContactVector = pysiChem(5.0, 4.0, self.unique_temp_folder, self.log)
        # mutations is of the form I_V70A
        # pysiChem need is like I70 (chain that is mutated and the position)
        # self.mutationss is a list which can contain more than one mutations
        mut_physChem = [ mut_data.chains_modeller[0] + mut_data.mutation_domain[1:-1] ]
        physChem_mut = list()
        physChem_wt = list()
        physChem_mut_ownChain = list()
        physChem_wt_ownChain = list()
#        os.chdir(self.unique_temp_folder) # from os

        for item in repairedPDB_wt_list:
            # calculate the contact vector
            res_wt = [0, 0, 0, 0]
            res_wt_ownChain = [0, 0, 0, 0]
            for mut in mut_physChem:
                chains_complex = [ chain for chain in mut_data.chains_modeller[0] ]
#                atomicContactVector, atomicContactVector_ownChain, mutpos = get_atomicContactVector(item, mut[0], chains_complex, mut_snippet_pos, mut_snippet)
                atomicContactVector, atomicContactVector_ownChain = \
                    get_atomicContactVector(item, mut[0], chains_complex,
                                            int(mut_data.mutation_domain[1:-1]),
                                            mut_data.mutation_domain[0])
                for index in range(0,4):
                    res_wt[index] += atomicContactVector[index]
                for index in range(0,4):
                    res_wt_ownChain[index] += atomicContactVector_ownChain[index]
            physChem_wt.append(res_wt)
            physChem_wt_ownChain.append(res_wt_ownChain)

        for item in repairedPDB_mut_list:
            # calculate the contact vector
            res_mut = [0, 0, 0, 0]
            res_mut_ownChain = [0, 0, 0, 0]
            for mut in mut_physChem:
                chains_complex = [ chain for chain in mut_data.chains_modeller[0] ]
#                atomicContactVector, atomicContactVector_ownChain, mutpos = get_atomicContactVector(item, mut[0], chains_complex, mut_snippet_pos, mut_snippet)
                atomicContactVector, atomicContactVector_ownChain = \
                    get_atomicContactVector(item, mut[0], chains_complex,
                                            int(mut_data.mutation_domain[1:-1]),
                                            mut_data.mutation_domain[-1])
                # what ever happens here:
                # atomicContactVector seems not to be a normal list..
                # thus I convert it to one in this ugly way...
                for index in range(0,4):
                    res_mut[index] += atomicContactVector[index]
                for index in range(0,4):
                    res_mut_ownChain[index] += atomicContactVector_ownChain[index]
            physChem_mut.append(res_mut)
            physChem_mut_ownChain.append(res_mut_ownChain)


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
            self.log.error('Wrong amino acid for msms wild-type!')
            self.log.error(seasa_info_wt)
            self.log.error(seasa_by_residue_separately)
            raise Exception('surface area calculated for the wrong atom!')
        solvent_accessibility_wt = seasa_info_wt['rel_sasa']

        secondary_structure_wt, solvent_accessibility_dssp_wt = analyze_structure_wt.get_dssp()
        secondary_structure_wt = secondary_structure_wt[mut_data.chains_modeller[0]][int(mut_data.mutation_domain[1:-1])-1]

        if isinstance(d, sql_db.UniprotDomainPair):
            contact_distance_wt = analyze_structure_wt.get_interchain_distances(mut_data.chains_modeller[0], mut_data.mutation_modeller)
            try:
                contact_distance_wt = contact_distance_wt[mut_data.chains_modeller[0]][0]
            except IndexError:
                self.log.error(contact_distance_wt)
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
            self.log.error('Wrong amino acid for msms mutant!')
            self.log.error(pdb_template.convert_aa(seasa_info_mut['res_name']))
            self.log.error(mut_data.mutation_domain[-1])
            self.log.error(seasa_info_mut)
            self.log.error(seasa_by_residue_separately)
            raise Exception('surface area calculated for the wrong atom!')
        solvent_accessibility_mut = seasa_info_mut['rel_sasa']

        secondary_structure_mut, solvent_accessibility_dssp_mut = analyze_structure_mut.get_dssp()
        secondary_structure_mut = secondary_structure_mut[mut_data.chains_modeller[0]][int(mut_data.mutation_domain[1:-1])-1]

        if isinstance(d, sql_db.UniprotDomainPair):
            contact_distance_mut = analyze_structure_wt.get_interchain_distances(mut_data.chains_modeller[0], mut_data.mutation_modeller)
            try:
                contact_distance_mut = contact_distance_mut[mut_data.chains_modeller[0]][0]
            except IndexError:
                self.log.error(contact_distance_mut)
                raise
#            except:
#                print "DSSP DID NOT WORK!"
#                self.log.error('DSSP did not work!')
#                solvent_accessibility_wt, secondary_structure_wt = '-1', '-1'
#                solvent_accessibility_mut, secondary_structure_mut = '-1', '-1'
#                mutation_errors += 'dssp did not work; '


        #######################################################################
        ## 10th: cobine the results
        # convert the elements to string and join them
        # they are lists in a list of the form: physChem_mut = [[0, 0, 0, 0], [0, 0, 0, 0]]
        # a little bit confusing but in the following are two list comprehensions
        # it is [item, item] with each item = [element, element, element, element]
        # first convert every element to str
        physChem_mut          = [ [ str(element) for element in item ] for item in physChem_mut ]
        physChem_mut_ownChain = [ [ str(element) for element in item ] for item in physChem_mut_ownChain ]
        physChem_wt           = [ [ str(element) for element in item ] for item in physChem_wt ]
        physChem_wt_ownChain  = [ [ str(element) for element in item ] for item in physChem_wt_ownChain ]
        # then join every 'element' via ':'
        physChem_mut          = ':'.join([','.join(item) for item in physChem_mut])
        physChem_mut_ownChain = ':'.join([','.join(item) for item in physChem_mut_ownChain])
        physChem_wt           = ':'.join([','.join(item) for item in physChem_wt])
        physChem_wt_ownChain  = ':'.join([','.join(item) for item in physChem_wt_ownChain])


        #######################################################################
        ## 11th: get the BLOSUM (or what ever matrix is given) score
        # self.mutationss is a list containing the mutations of the form A_H56T
        matrix_score = 0
        for mut in [mut_data.mutation_domain]:
            fromAA = mut[0]
            toAA   = mut[-1]
            matrix_score += score_pairwise(fromAA, toAA, self.matrix, self.gap_s, self.gap_e)


        #######################################################################
        #

        if isinstance(d, sql_db.UniprotDomain):
            uniprot_mutation = sql_db.UniprotDomainMutation()
            uniprot_mutation.uniprot_domain_id = t.uniprot_domain_id

        elif isinstance(d, sql_db.UniprotDomainPair):
            uniprot_mutation = sql_db.UniprotDomainPairMutation()
            uniprot_mutation.uniprot_domain_pair_id = t.uniprot_domain_pair_id
            uniprot_mutation.AnalyseComplex_energy_wt = AnalyseComplex_energy_wt
            uniprot_mutation.AnalyseComplex_energy_mut = AnalyseComplex_energy_mut
            uniprot_mutation.contact_distance_wt = contact_distance_wt
            uniprot_mutation.contact_distance_mut = contact_distance_mut

        uniprot_mutation.uniprot_id = uniprot_id_1
        uniprot_mutation.mutation = mutation

        uniprot_mutation.chain_modeller = mut_data.chains_modeller[0]
        uniprot_mutation.mutation_modeller = mut_data.mutation_modeller

        uniprot_mutation.model_filename_wt = model_filename_wt
        uniprot_mutation.model_filename_mut = model_filename_mut

        uniprot_mutation.Stability_energy_wt = Stability_energy_wt
        uniprot_mutation.Stability_energy_mut = Stability_energy_mut

        uniprot_mutation.physChem_wt = physChem_wt
        uniprot_mutation.physChem_wt_ownChain = physChem_wt_ownChain
        uniprot_mutation.physChem_mut = physChem_mut
        uniprot_mutation.physChem_mut_ownChain = physChem_mut_ownChain

        uniprot_mutation.matrix_score = matrix_score

        uniprot_mutation.secondary_structure_wt = secondary_structure_wt
        uniprot_mutation.solvent_accessibility_wt = solvent_accessibility_wt
        uniprot_mutation.secondary_structure_mut = secondary_structure_mut
        uniprot_mutation.solvent_accessibility_mut = solvent_accessibility_mut

        uniprot_mutation.provean_score = provean_score
        #######################################################################
        # Save alignments and modeller models to output database for storage
        # Move template files to the output folder as a tar archives
#        self.make_tarfile(self.HOME + self.outputPath + save_path + '_' + mutation + '.tar.bz2',
#                          save_path[:-1])

        #######################################################################
        self.log.info('Finished processing template:')
        self.log.info(mut_data.save_path.split('/')[-2])

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
            ' --supporting_set ' + path_to_provean_supset)
        child_process = hf.run_subprocess_locally(
            self.unique_temp_folder + 'sequence_conservation/',
            system_command)
        result, error_message = child_process.communicate()
        self.log.debug(result)
        if child_process.returncode != 0:
            self.log.error(error_message)
            self.log.error(system_command)
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



