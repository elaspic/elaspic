# -*- coding: utf-8 -*-
import os.path as op
import logging 
import json

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import conf, errors, structure_tools, structure_analysis, sequence, model, predictor
from .pipeline import Pipeline, execute_and_remember

logger = logging.getLogger(__name__)
configs = conf.Configs()

sql_db = None
domain_alignment = None
domain_model = None
domain_mutation = None



#%%
class LocalPipeline(Pipeline):

    def __init__(self, structure_file, sequence_file=None, mutations=None, configurations=None):
        """
        TODO: Add an option to store provean results based on sequence hash.
        """
        super().__init__(configurations)
        
        # Input parameters
        self.pdb_id = op.splitext(op.basename(structure_file))[0]
        self.pdb_file = structure_file
        self.sequence_file = sequence_file if sequence_file else ''

        logger.info('pdb_file: {}'.format(self.pdb_file))
        logger.info('pwd: {}'.format(self.PWD))
        
        ### Load PDB structure and extract required sequences and chains.
        #fix_pdb(self.pdb_file, self.pdb_file)
        self.sp = structure_tools.StructureParser(self.pdb_file)
        self.sp.extract()
        self.sp.save_structure(configs['unique_temp_dir'])
        self.sp.save_sequences(configs['unique_temp_dir'])
        
        if not mutations:
            mutations = []
        elif isinstance(mutations, str):
            mutations = mutations.split(',')
        else:
            mutations = mutations

        if self.sequence_file is None:
            # Read template sequences from the PDB
            self.seqrecords = tuple([
                SeqRecord(
                    id='{}{}_{}'.format(self.pdb_id, chain_id, i), 
                    seq=Seq(self.sp.chain_sequence_dict[chain_id]))
                for (i, chain_id) in enumerate(self.sp.chain_ids)
            ])
            ### Parse mutations
            self.mutations = dict()
            for mutation_in in mutations:
                mutation_chain, mutation_resnum = mutation_in.split('_')
                # This is the index of the chain that is being mutated
                mutation_idx = self.sp.chain_ids.index(mutation_chain)
                # Converting mutation from PDB to index coordinates
                mutation_id = (' ', int(mutation_resnum[1:-1]), ' ',)
                chain_aa_residues = (
                    structure_tools.get_aa_residues(self.sp.structure[0][mutation_chain])
                )
                mutation_pos = chain_aa_residues.index(mutation_id) + 1
                mutation = mutation_resnum[0] + str(mutation_pos) + mutation_resnum[-1]
                # Validating
                mutation_expected_aa = (
                    structure_tools.AAA_DICT
                    [self.sp.structure[0][mutation_chain][mutation_id].resname]
                )
                if mutation_resnum[0] != mutation_expected_aa:
                    raise errors.MutationMismatchError()
                self.mutations[(mutation_idx, mutation,)] = mutation_in
        else:
            # Read template sequences from the sequence file
            self.seqrecords = list(SeqIO.parse(self.sequence_file, 'fasta'))
            for i, seqrec in enumerate(self.seqrecords):
                seqrec.id = '{}_{}'.format(seqrec.id, str(i))
            ### Parse mutations
            self.mutations = []
            for mutation_in in mutations:
                mutation_idx, mutation = mutation_in.split('_')
                mutation_idx = int(mutation_idx)
                # Validation
                mutation_expected_aa = (
                    str(self.seqrecords[mutation_idx].seq)[int(mutation[1:-1])]
                )
                if mutation[0] != mutation_expected_aa:
                    raise errors.MutationMismatchError()
                self.mutations[(mutation_idx, mutation,)] = mutation_in


    #==============================================================================
    # Run methods
    #==============================================================================

    def run(self):
        self.run_all_sequences()
        self.run_all_models()
        self.run_all_mutations()
    
    
    def run_all_sequences(self):
        sequence_results = []
        sequence_results_file = op.join(configs['unique_temp_dir'], 'sequence.json')
        if op.isfile(sequence_results_file):
            logger.debug('Results file for model already exists: {}'.format(sequence_results_file))
            return
        for chain_id in self.sp.chain_ids:
            if chain_id == self.sp.hetatm_chain_id:
                continue
            idx = self._get_chain_idx(chain_id)
            sequence = self.get_sequence(idx)
            sequence_result = sequence.result
            sequence_result['idx'] = idx
            sequence_results.append(sequence_result)
        with open(sequence_results_file, 'w') as ofh:
            json.dump(sequence_results, ofh)
    
    
    def run_all_models(self):
        model_results = []
        model_results_file = op.join(configs['unique_temp_dir'], 'model.json')
        if op.isfile(model_results_file):
            logger.debug('Results file for model already exists: {}'.format(model_results_file))
            return
        for chain_id in self.sp.chain_ids:
            if chain_id == self.sp.hetatm_chain_id:
                continue
            idx = self._get_chain_idx(chain_id)
            model = self.get_model(idx)
            model_result = model.result
            model_result['idx'] = self.sp.chain_ids.index(chain_id)
            model_results.append(model_result)
        for idxs in self.sp.interacting_chain_idxs:
            model = self.get_model(idxs)
            model_result = model.result
            model_result['idxs'] = tuple(idxs)
            model_results.append(model_result)
        with open(model_results_file, 'w') as ofh:
            json.dump(model_results, ofh)
    
    
    def run_all_mutations(self):
        for (mutation_idx, mutation), mutation_in in self.mutations.items():
            mutation_results = []            
            mutation_results_file = op.join(
                configs['unique_temp_dir'], 'mutation_{}.json'.format(mutation_in)
            )
            if op.isfile(mutation_results_file):
                logger.debug(
                    'Results file for mutation {} already exists: {}'
                    .format(mutation_in, mutation_results_file)
                )
                continue
            mutation_result = self.get_mutation_score(mutation_idx, mutation_idx, mutation)
            mutation_result['idx'] = mutation_idx
            mutation_results.append(mutation_result)
            for idxs in self.sp.interacting_chain_idxs:
                if mutation_idx in idxs:
                    mutation_result = self.get_mutation_score(idxs, mutation_idx, mutation)
                    mutation_result['idx'] = mutation_idx
                    mutation_result['idxs'] = tuple(idxs)
                    mutation_results.append(mutation_result)
            with open(mutation_results_file, 'w') as ofh:
                json.dump(mutation_results, ofh)


    #==============================================================================
    # Get methods
    #==============================================================================

    def get_sequence(self, idx):
        """
        """
        logger.debug('-' * 80)
        logger.debug('get_sequence({})'.format(idx))
        return PrepareSequence(self.seqrecords, idx, None)


    def get_model(self, idxs):
        """
        Use modeller to make a homology model for each uniprot domain that
        has a template in pdbfam
        """
        logger.debug('-' * 80)
        logger.debug('get_model({})'.format(idxs))
        idxs = self._sort_chain_idxs(idxs)
        return PrepareModel(self.seqrecords, self.sp, idxs)


    def get_mutation_score(self, idxs, mutation_idx, mutation):
        logger.debug('-' * 80)
        logger.debug('get_mutation_score({}, {}, {})'.format(idxs, mutation_idx, mutation))
        idxs = self._sort_chain_idxs(idxs)
        sequence = self.get_sequence(mutation_idx)
        model = self.get_model(idxs)
        return PrepareMutation(sequence, model, idxs.index(mutation_idx), mutation)
        
        
    #==============================================================================
    # Helper functions
    #==============================================================================
    
    def _get_chain_idx(self, chain_id):
        """chain_id -> chain_idx
        """
        chain_idx = [
            i for (i, chain) 
            in enumerate(self.sp.structure[0].child_list) 
            if chain.id == chain_id
        ]
        if len(chain_idx) == 0:
            raise errors.PDBChainError(
                'Chain {} was not found in PDB {}!'.format(chain_id, self.sp.pdb_file))
        elif len(chain_idx) > 1:
            raise errors.PDBChainError(
                'Chain {} was found more than once in PDB {}!'.format(chain_id, self.sp.pdb_file))
        return chain_idx[0]
        
    
    def _sort_chain_idxs(self, idxs):
        """Sort positions and return as as tuple.
        """
        if not hasattr(idxs, '__getitem__'):
            idxs = (idxs,)
        else:
            idxs = tuple(sorted(idxs))
        return idxs
            


#%%
@execute_and_remember
class PrepareSequence:
    """
    Raises
    -------
    errors.ProveanError
    errors.ProveanResourceError
    """
    def __init__(self, seqrecords, position, provean_supset_file):
        self.seqrecord = seqrecords[position]
        self.provean_supset_file = provean_supset_file
        self.sequence_file = None
        self.sequence = None
        
    def __bool__(self):
        return True
        
    def __enter__(self):
        sequence_file = op.join(configs['sequence_dir'], self.seqrecord.id + '.fasta')
        with open(sequence_file, 'w') as ofh:
            SeqIO.write(self.seqrecord, ofh, 'fasta')
        self.sequence_file = sequence_file
        
    def run(self):
        self.sequence = sequence.Sequence(self.sequence_file, self.provean_supset_file)
        
    def __exit__(self, exc_type, exc_value, traceback):
        return False
        
    @property
    def result(self):
        return self.sequence



#%%
@execute_and_remember
class PrepareModel:
    """
    Returns
    -------
    model_dict : dict
        Contains all important values calculated by modeller.
        
    Raises
    -------
    errors.ModellerError
    errors.PDBChainError
    errors.PDBEmptySequenceError
    errors.PDBNotFoundError
    """

    def __init__(self, seqrecords, sp, positions):
        self.seqrecords = [seqrecords[pos] for pos in positions]
        self.positions = positions
        self.sp = sp
        self.sequence_file = None
        self.structure_file = None
        self.model = None
        
    def __bool__(self):
        return True
        
    def __enter__(self):
        # Target sequence file
        self.sequence_file = op.join(
            configs['model_dir'], '_'.join(seqrec.id for seqrec in self.seqrecords) + '.fasta'
        )
        with open(self.sequence_file, 'w') as ofh:
            SeqIO.write(self.seqrecords, ofh, 'fasta')
        assert op.isfile(self.sequence_file)

        # Template structure file
        chain_string = ''.join(self.sp.structure[0].child_list[pos].id for pos in self.positions)
        self.structure_file = (
            op.join(configs['unique_temp_dir'], self.sp.pdb_id + chain_string + '.pdb')
        )
        assert op.isfile(self.structure_file)

    def run(self):
        self.model = model.Model(self.sequence_file, self.structure_file)
        
    def __exit__(self, exc_type, exc_value, traceback):
        return False
        
    @property
    def result(self):
        return self.model



#%%
@execute_and_remember
class PrepareMutation:
    """
    Raises
    ------
    errors.PDBError
    errors.FoldxError
    errors.ResourceError
    errors.FoldXAAMismatchError
    errors.MutationOutsideDomainError
    errors.MutationOutsideInterfaceError

    """

    def __init__(self, sequence, model, mutation_idx, mutation):
        self.sequence = sequence
        self.model = model
        self.mutation_idx = mutation_idx
        self.mutation = mutation
    
    def __bool__(self):
        return True
    
    def __enter__(self):
        pass
        
    def run(self):
        features = dict()
        
        ### Sequence features
        results = self.sequence.mutate(self.mutation)
        features['provean_score'] = results['provean_score']
        features['matrix_score'] = results['matrix_score']        

        ### Structure features
        results = self.model.mutate(self.mutation_idx, self.mutation)
        
        features['norm_dope'] = self.model.modeller_results['norm_dope']
        (features['alignment_identity'], 
         features['alignment_coverage'], 
         features['alignment_score']) = (
             self.model.modeller_results['alignment_stats'][self.mutation_idx]
        )
        
        features['model_filename_wt'] = results['model_filename_wt']
        features['model_filename_mut'] = results['model_filename_mut']

        features['stability_energy_wt'] = results['stability_energy_wt']
        features['stability_energy_mut'] = results['stability_energy_mut']

        features['physchem_wt'] = (
            '{},{},{},{}'.format(*results['physchem_wt'])
        )
        features['physchem_wt_ownchain'] = (
            '{},{},{},{}'.format(*results['physchem_ownchain_wt'])
        )
        features['physchem_mut'] = (
            '{},{},{},{}'.format(*results['physchem_mut'])
        )
        features['physchem_mut_ownchain'] = (
            '{},{},{},{}'.format(*results['physchem_ownchain_mut'])
        )
        
        features['secondary_structure_wt'] = results['secondary_structure_wt']
        features['solvent_accessibility_wt'] = results['solvent_accessibility_wt']
        features['secondary_structure_mut'] = results['secondary_structure_mut']
        features['solvent_accessibility_mut'] = results['solvent_accessibility_mut']
        
        if len(self.model.sequence_seqrecords) > 1:
            features['interface_area_hydrophobic'] = self.model.interface_area_hydrophobic
            features['interface_area_hydrophilic'] = self.model.interface_area_hydrophilic
            features['interface_area_total'] = self.model.interface_area_total

            features['analyse_complex_energy_wt'] = results['analyse_complex_energy_wt']
            features['analyse_complex_energy_mut'] = results['analyse_complex_energy_mut']
            features['contact_distance_wt'] = results['contact_distance_wt']
            features['contact_distance_mut'] = results['contact_distance_mut']

                
        logger.debug('feature_dict: {}'.format(features))
        feature_df = pd.DataFrame(features, index=[0])
    
        pred = predictor.Predictor()
        features['ddg'] = pred.score(feature_df, len(self.model.sequence_seqrecords) > 1)
        logger.debug('Predicted ddG: {}'.format(features['ddg']))
        
        self.mutation_features = features
        

    def __exit__(self, exc_type, exc_value, traceback):
        return False


    @property
    def result(self):
        return self.mutation_features



#%%
class Unused:
    
    def get_structure_sequence(self):
        """
        """
        structure_sequence = []
        # After running `sp.extract()`, all but the last chain are guaranteed to be without HATATMS
        for chain in self.sp.structure.child_list[0].child_list[:-1]:
            chain_sequence, chain_numbering_extended = self.sp.get_chain_sequence_and_numbering(chain.id)
            structure_sequence += [chain_sequence]
        
        # If the last chain starts with a HETATM, it should be entirely HETATMs.
        last_chain = self.sp.structure.child_list[0].child_list[-1]
        if last_chain.child_list[0].resname in structure_tools.AAA_DICT:
            chain_sequence, chain_numbering_extended = self.sp.get_chain_sequence_and_numbering(last_chain.id)
            structure_sequence += [chain_sequence]
        else:
            structure_sequence += ['.' * len(last_chain.child_list)]
            
        structure_sequence = '/'.join(structure_sequence)
        return structure_sequence
    
        
    
    def get_mutation_data(self, pdb_chain, pdb_mutation, provean_results, modeller_results):
        """
        Create a MutData class that holds all the information we need
        about a given mutation.
        """    
        protein_id = "{}{}".format(self.pdb_id, pdb_chain)
        logger.info(provean_results)
        provean_supset_filename = provean_results['provean_supset_filename']
        model_file = modeller_results['model_file']
        pdb_chains = [chain.id for chain in self.sp.structure[0].child_list]

        mutation_chain_pos = 1
        for aa, aa_pos in zip(*self.sp.get_chain_sequence_and_numbering(pdb_chain)):
            if aa_pos != pdb_mutation[1:-1]:
                mutation_chain_pos += 1
            else:
                break
        mutation_chain = pdb_mutation[0] + str(mutation_chain_pos) + pdb_mutation[-1]
        
        mutation_modeller = self.mutation_to_modeller(pdb_chain, pdb_mutation)
        mutation_modeller_pos = mutation_modeller[1:-1]

        ### Assign values to appropriate variables
        mut_data = domain_mutation.MutationData()
        
        # Mutation info
        mut_data.uniprot_id_1 = protein_id 
        mut_data.mutation = mutation_modeller
        mut_data.chains_modeller = pdb_chains
        mut_data.domain_sequences = [
            self.sp.get_chain_sequence_and_numbering(chain_id)[0]
            for chain_id in pdb_chains
        ]
        mut_data.path_to_provean_supset = op.join(
            self.unique_temp_folder, 'sequence_conservation', provean_supset_filename
        )    
        mut_data.save_path = self.unique_temp_folder
        mut_data.pdbFile_wt = op.join(self.unique_temp_folder, 'modeller', model_file)
        mut_data.position_domain = mutation_chain_pos
        mut_data.mutation_domain = mutation_chain
        mut_data.position_modeller = mutation_modeller_pos
        mut_data.mutation_modeller = mutation_modeller
        
        # Provean results
        mut_data.provean_mutation = None
        mut_data.provean    
        
        
    
    def _mutation_to_modeller(self, pdb_chain, pdb_mutation):
        """
        Modeller numbers all residues from 1 to n_residues, for all chains. 
        Use this function to convert mutations from PDB coordinates to Modeller coordinates.
        """
        mutation_id = pdb_mutation[:-1]
        mutation_idx = 1
        
        for chain in self.sp.structure.child_list[0].child_list:
            if chain.id != pdb_chain:
                mutation_idx += len(chain.child_list)
            else:
                break
            
        for residue in chain.child_list:
            residue_id = (structure_tools.AAA_DICT[residue.resname] + str(residue.id[1]))
            if residue_id != mutation_id:
                mutation_idx += 1
            else:
                break
    
        if residue_id != mutation_id:
            return None
        else:
            mutation_modeller = "{}{}{}".format(pdb_mutation[0], mutation_idx, pdb_mutation[-1])
            return mutation_modeller
    
    
    
    def _get_interactions(self, pdb_chain, pdb_mutation):
        interactions = structure_analysis.get_interactions(
            self.sp.structure.child_list[0], pdb_chain, self.R_CUTOFF)

        mutation_contacts = []
        for chain in interactions:
            for aa_pos, aa in interactions[chain]:
                if aa + aa_pos == pdb_mutation[:-1]:
                    mutation_contacts.append(chain)

        return mutation_contacts

