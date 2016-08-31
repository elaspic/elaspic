""".

TODO: The model object has two serialization steps:
    1. Inside the modeller class to save modeller results.
    2. In the local_pipeline to save all results.
"""
import os.path as op
import logging
import json

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import (
    CACHE_DIR, conf, helper, errors, structure_tools, elaspic_sequence,
    elaspic_model, elaspic_predictor
)
from .pipeline import Pipeline, execute_and_remember

logger = logging.getLogger(__name__)

sql_db = None
domain_alignment = None
domain_model = None
domain_mutation = None


class StandalonePipeline(Pipeline):
    """Run elaspic locally without a database.

    Parameters
    ----------
    structure_file : str
        Full path to the structure *pdb* file.
    sequence_file : str
        Full path to the sequence *fasta* file.
    mutations : str
        Comma-separated list of mutations.
    configurations : str
        Full path to the configurations *ini* file.
    mutation_format : int, default None
        1. {pdb_chain}_{pdb_mutation},...
        2. {pdb_chain}_{sequence_mutation},...
        3. {sequence_pos}_{sequence_mutation}...

        If `sequence_file` is None, this does not matter (always {pdb_chain}_{pdb_mutation}).

    .. todo:: Add an option to store provean results based on sequence hash.
    """

    def __init__(
            self, structure_file, sequence_file=None, mutations=None, configurations=None,
            mutation_format=None, run_type='5'):
        super().__init__(configurations)

        # Input parameters
        self.pdb_id = op.splitext(op.basename(structure_file))[0]
        self.pdb_file = structure_file
        self.run_type = self._validate_run_type(run_type)

        logger.info('pdb_file: {}'.format(self.pdb_file))
        logger.info('pwd: {}'.format(self.PWD))

        # Load PDB structure and extract required sequences and chains.
        # fix_pdb(self.pdb_file, self.pdb_file)
        self.sp = structure_tools.StructureParser(self.pdb_file)
        self.sp.extract()
        self.sp.save_structure(conf.CONFIGS['unique_temp_dir'])
        self.sp.save_sequences(conf.CONFIGS['unique_temp_dir'])

        if sequence_file in ['', 'None', None]:
            self.sequence_file = ''
        else:
            self.sequence_file = sequence_file
        logger.debug('self.sequence_file: {}'.format(self.sequence_file))

        # Use the PDB chain to index mutations both with and without the index file
        if not self.sequence_file:
            # Read template sequences from the PDB
            self.seqrecords = tuple([
                SeqRecord(
                    id='{}{}_{}'.format(self.pdb_id, chain_id, i),
                    seq=Seq(self.sp.chain_sequence_dict[chain_id]))
                for (i, chain_id) in enumerate(self.sp.chain_ids)
            ])
        else:
            # Read template sequences from the sequence file
            self.seqrecords = tuple(SeqIO.parse(self.sequence_file, 'fasta'))
            for i, seqrec in enumerate(self.seqrecords):
                seqrec.id = helper.slugify('{}_{}'.format(seqrec.id, str(i)))

        self.mutations = self._split_mutations(mutations)
        if 'mutation' in self.run_type:
            self.mutations = self.parse_mutations(self.mutations, mutation_format)
        logger.debug('mutations: {}'.format(self.mutations))

        if len(self.sp.chain_ids) != len(self.seqrecords):
            logger.warning(
                'The number of chain ids ({}) does not match the number of sequences ({})!'
                .format(len(self.sp.chain_ids), len(self.seqrecords))
            )

    def parse_mutations(self, mutations, mutation_format):
        # Parse mutations
        # There are many ways mutations can be specified here...
        # try one at at a time until something succeeds
        parsed_mutations = dict()

        if self.sequence_file:
            possible_mutation_formats = ['3', '2', '1']
        else:
            possible_mutation_formats = ['1', '2', '3']

        if mutation_format is not None:
            parsed_mutations = self._parse_mutations(mutations, mutation_format)
        else:
            for mutation_format in possible_mutation_formats:
                try:
                    parsed_mutations = self._parse_mutations(mutations, mutation_format)
                    break
                except (IndexError, ValueError, errors.MutationMismatchError) as e:
                    error_message = (
                        "Error parsing mutations '{}' using mutation_format '{}':\n{} {}"
                        .format(mutations, mutation_format, type(e), e)
                    )
                    logger.error(error_message)
                    continue
            if not parsed_mutations:
                raise errors.MutationMismatchError()
        logger.debug('parsed mutations: {}'.format(self.mutations))
        return parsed_mutations

    def _parse_mutations(self, mutations, mutation_format='1'):
        """Parse mutations.

        Parse mutations provided using the {pdb_chain}_{mutation} naming scheme (default)
        or the {mutation_pos}_{mutation} naming scheme (fallback for when `pdb_chain`
        is not found in the structure).
        """
        logger.debug("Parsing mutations using mutation_format: '{}'".format(mutation_format))
        mutations_out = dict()
        for mutation_in in mutations:
            mutation_chain, mutation_residue = mutation_in.split('_')

            if mutation_format in ['1', '2']:
                # This is the index of the chain that is being mutated
                mutation_idx = self.sp.chain_ids.index(mutation_chain)
            elif mutation_format in ['3']:
                # Mutation pos starts at 1
                mutation_idx = int(mutation_chain) - 1

            if mutation_format in ['1']:
                # Converting mutation from PDB to sequence coordinates
                resnum = int(''.join(n for n in mutation_residue[1:-1] if n.isdigit()))
                resnum_suffix = (
                    ''.join(n for n in mutation_residue[1:-1] if not n.isdigit()).upper()
                )
                mutation_id = (' ', resnum, resnum_suffix or ' ')
                chain_aa_residues = (
                    structure_tools.get_aa_residues(self.sp.structure[0][mutation_chain])
                )
                mutation_pos = chain_aa_residues.index(mutation_id) + 1
                mutation = mutation_residue[0] + str(mutation_pos) + mutation_residue[-1]
                # Validation
                mutation_expected_aa = (
                    structure_tools.AAA_DICT
                    [self.sp.structure[0][mutation_chain][mutation_id].resname]
                )
                if mutation_residue[0] != mutation_expected_aa:
                    raise errors.MutationMismatchError()
            elif mutation_format in ['2', '3']:
                # Mutation is already in sequence coordinates
                mutation = mutation_residue
                # Validation
                mutation_expected_aa = (
                    str(self.seqrecords[mutation_idx].seq)[int(mutation[1:-1]) - 1]
                )
                if mutation[0] != mutation_expected_aa:
                    raise errors.MutationMismatchError()
                mutations_out[(mutation_idx, mutation,)] = mutation_in

            mutations_out[(mutation_idx, mutation,)] = mutation_in
        return mutations_out

    # === Run methods ===

    def run(self):
        if 'sequence' in self.run_type:
            self.run_all_sequences()
        if 'model' in self.run_type:
            self.run_all_models()
        if 'mutation' in self.run_type:
            self.run_all_mutations()

    def run_all_sequences(self):
        sequence_results = []
        sequence_results_file = op.join(conf.CONFIGS['unique_temp_dir'], 'sequence.json')
        if op.isfile(sequence_results_file):
            logger.debug('Results file for sequence already exists: {}'
                         .format(sequence_results_file))
            return
        for chain_id, _ in zip(self.sp.chain_ids, self.seqrecords):
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
        model_results_file = op.join(conf.CONFIGS['unique_temp_dir'], 'model.json')
        if op.isfile(model_results_file):
            logger.debug('Results file for model already exists: {}'.format(model_results_file))
            return
        for chain_id, _ in zip(self.sp.chain_ids, self.seqrecords):
            if chain_id == self.sp.hetatm_chain_id:
                continue
            idx = self._get_chain_idx(chain_id)
            model = self.get_model(idx)
            model_result = model.result
            model_result['idx'] = self.sp.chain_ids.index(chain_id)
            model_results.append(model_result)
        for idxs in self.sp.interacting_chain_idxs:
            if not all(i in range(len(self.seqrecords)) for i in idxs):
                warning = (
                    "Skipping idxs: '{}' because we lack the corresponding seqrecord!"
                    .format(idxs)
                )
                logger.warning(warning)
                continue
            model = self.get_model(idxs)
            if model is None:
                continue
            model_result = model.result
            model_result['idxs'] = tuple(idxs)
            model_results.append(model_result)
        with open(model_results_file, 'w') as ofh:
            json.dump(model_results, ofh)

    def run_all_mutations(self):
        handled_errors = (
            errors.ChainsNotInteractingError,
            errors.MutationOutsideDomainError,
            errors.MutationOutsideInterfaceError,
        )
        for (mutation_idx, mutation), mutation_in in self.mutations.items():
            mutation_results = []
            mutation_results_file = op.join(
                conf.CONFIGS['unique_temp_dir'], 'mutation_{}.json'.format(mutation_in)
            )
            if op.isfile(mutation_results_file):
                logger.debug(
                    'Results file for mutation {} already exists: {}'
                    .format(mutation_in, mutation_results_file)
                )
                continue
            try:
                mutation_result = self.get_mutation_score(mutation_idx, mutation_idx, mutation)
            except handled_errors as e:
                logger.error(e)
                continue
            mutation_result['idx'] = mutation_idx
            mutation_results.append(mutation_result)
            for idxs in self.sp.interacting_chain_idxs:
                if not all(i in range(len(self.seqrecords)) for i in idxs):
                    warning = (
                        "Skipping idxs: '{}' because we lack the corresponding seqrecord!"
                        .format(idxs)
                    )
                    logger.warning(warning)
                    continue
                if mutation_idx in idxs:
                    try:
                        mutation_result = self.get_mutation_score(idxs, mutation_idx, mutation)
                    except handled_errors as e:
                        logger.error(e)
                        continue
                    mutation_result['idx'] = mutation_idx
                    mutation_result['idxs'] = tuple(idxs)
                    mutation_results.append(mutation_result)
            with open(mutation_results_file, 'w') as ofh:
                json.dump(mutation_results, ofh)

    # === Get methods ===

    def get_sequence(self, idx):
        """
        """
        logger.debug('-' * 80)
        logger.debug('get_sequence({})'.format(idx))
        return PrepareSequence(self.seqrecords, idx, None)

    def get_model(self, idxs):
        """Make a homology model for each uniprot domain with a template in pdbfam."""
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

    # === Helper functions ===

    def _get_chain_idx(self, chain_id):
        """chain_id -> chain_idx."""
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
        """Sort positions and return as as tuple."""
        if not hasattr(idxs, '__getitem__'):
            idxs = (idxs,)
        else:
            idxs = tuple(sorted(idxs))
        return idxs


@execute_and_remember
class PrepareSequence:
    """.

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
        sequence_file = op.join(
            conf.CONFIGS['sequence_dir'],
            helper.slugify(self.seqrecord.id + '.fasta'))
        with open(sequence_file, 'w') as ofh:
            SeqIO.write(self.seqrecord, ofh, 'fasta')
        self.sequence_file = sequence_file

    def run(self):
        self.sequence = elaspic_sequence.Sequence(self.sequence_file, self.provean_supset_file)

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    @property
    def result(self):
        return self.sequence


@execute_and_remember
class PrepareModel:
    """.

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
            conf.CONFIGS['model_dir'],
            helper.slugify('_'.join(seqrec.id for seqrec in self.seqrecords) + '.fasta')
        )
        with open(self.sequence_file, 'w') as ofh:
            SeqIO.write(self.seqrecords, ofh, 'fasta')
        assert op.isfile(self.sequence_file)

        # Template structure file
        chain_string = ''.join(self.sp.structure[0].child_list[pos].id for pos in self.positions)
        self.structure_file = op.join(
            conf.CONFIGS['unique_temp_dir'],
            helper.slugify(self.sp.pdb_id + chain_string + '.pdb')
        )
        assert op.isfile(self.structure_file)

    def run(self):
        self.model = elaspic_model.Model(self.sequence_file, self.structure_file)

    def __exit__(self, exc_type, exc_value, traceback):
        handled_errors = (
            errors.ChainsNotInteractingError,
        )
        if exc_type in handled_errors:
            logger.debug("Caught the following error: '{}'".format(exc_type))
            return True
        return False

    @property
    def result(self):
        return self.model


@execute_and_remember
class PrepareMutation:
    """.

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
        if not self.sequence or not self.model:
            raise errors.ChainsNotInteractingError

        features = dict()
        features['mutation'] = self.mutation

        # Sequence features
        results = self.sequence.mutate(self.mutation)
        features['provean_score'] = results['provean_score']
        features['matrix_score'] = results['matrix_score']

        # Structure features
        results = self.model.mutate(self.mutation_idx, self.mutation)

        features['norm_dope'] = self.model.modeller_results['norm_dope']
        (features['alignment_identity'],
         features['alignment_coverage'],
         features['alignment_score']) = (
             self.model.modeller_results['alignment_stats'][self.mutation_idx]
        )
        assert features['alignment_identity'] > 0.01 and features['alignment_identity'] <= 1
        assert features['alignment_coverage'] > 0.01 and features['alignment_coverage'] <= 1

        features['model_file_wt'] = results['model_file_wt']
        features['model_file_mut'] = results['model_file_mut']

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

        # new additions
        features['mutation_errors'] = results['mutation_errors']
        features['chain_modeller'] = results['chain_modeller']
        features['mutation_modeller'] = results['mutation_modeller']

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

        if len(self.model.sequence_seqrecords) == 1:
            pred = elaspic_predictor.CorePredictor()
        else:
            pred = elaspic_predictor.InterfacePredictor()
        pred.load(CACHE_DIR)
        features['ddg'] = pred.score(feature_df)[0]
        logger.debug('Predicted ddG: {}'.format(features['ddg']))

        self.mutation_features = features

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    @property
    def result(self):
        return self.mutation_features
