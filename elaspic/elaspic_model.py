import os
import os.path as op
import logging
import shutil
import json
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBIO
from kmtools.system_tools import switch_paths
from . import (
    conf, errors, structure_tools, structure_analysis,
    call_modeller, call_tcoffee, call_foldx
)

logger = logging.getLogger(__name__)


class Model:
    """Structural homology model.

    Parameters
    ----------
    sequence_file
        fasta file containing the sequence of the protein that should be mutated.
    structure_file
        pdb file containing the structure to be used as a template for homology modelling.
    modeller_results_file
        Precalculated data from a previous modeller run.
    """

    def __init__(self, sequence_file, structure_file, modeller_results_file=None):
        logger.debug('Initialising a Model instance with parameters:')
        logger.debug('sequence_file: {}:'.format(sequence_file))
        logger.debug('structure_file: {}:'.format(structure_file))

        # Target sequences
        self.sequence_file = sequence_file
        self.sequence_seqrecords = list(SeqIO.parse(self.sequence_file, 'fasta'))
        self.sequence_id = op.splitext(op.basename(self.sequence_file))[0].replace(':', '.')
        self._validate_sequence_seqrecords()
        logger.debug('sequence_seqrecords: {}'.format(self.sequence_seqrecords))

        # Template structures
        self.structure_file = structure_file
        self.structure = structure_tools.get_pdb_structure(self.structure_file)
        self.structure_id = self.structure.id.replace(':', '.')
        self.structure_seqrecords = [
            SeqRecord(
                id='{}{}'.format(self.structure_id, chain.id),
                seq=Seq(
                    structure_tools
                    .get_chain_sequence_and_numbering(chain, include_hetatms=True)[0])
            ) for chain in self.structure[0].child_list
        ]
        self.chain_ids = [chain.id for chain in self.structure.child_list[0].child_list]
        logger.debug('structure_seqrecords: {}'.format(self.structure_seqrecords))

        # Homology modelling
        if self.sequence_id == self.structure_id:
            self.sequence_id += '_sequence'
        self.model_id = '{}-{}'.format(self.sequence_id, self.structure_id)

        # Check for precalculated data
        self.modeller_results_file = op.join(conf.CONFIGS['model_dir'], self.model_id + '.json')
        if (modeller_results_file is not None and
                modeller_results_file != self.modeller_results_file):
            logger.debug(
                'Copying precalculated modeller results file from {} to {}...'
                .format(modeller_results_file, self.modeller_results_file)
            )
            shutil.copy(modeller_results_file, self.modeller_results_file)
        if op.isfile(self.modeller_results_file):
            logger.debug(
                'Loading precalculated modeller results from file: {}'
                .format(self.modeller_results_file)
            )
            with open(self.modeller_results_file) as ifh:
                self.modeller_results = json.load(ifh)
        else:
            logger.debug('Creating sequence alignments and building a homology model')
            self._create_alignments_and_model()
            # Save model into a json file for faster future use
            with open(self.modeller_results_file, 'w') as ofh:
                json.dump(self.modeller_results, ofh)

        # Get interacting amino acids and interface area
        self.modeller_structure = (
            structure_tools.get_pdb_structure(
                op.join(conf.CONFIGS['unique_temp_dir'], self.modeller_results['model_file']))
        )
        self.modeller_chain_ids = [
            chain.id for chain in self.modeller_structure[0]
        ]
        self._analyse_core()
        if len(self.sequence_seqrecords) > 1:
            self._analyse_interface()

        self.mutations = {}
        self.errors = []

    @property
    def core_or_interface(self):
        if len(self.sequence_seqrecords) == 1:
            return 'core'
        else:
            return 'interface'

    def _validate_sequence_seqrecords(self):
        if len(self.sequence_seqrecords) > 2:
            message = (
                "ELASPIC is designed to predict the effect of mutations on the folding "
                "of a single domain or the interaction between two domains. It cannot predict "
                "the effect of mutations on the interaction between more than two domains. "
            )
            logger.warning(message)

    def _align_with_tcoffee(self, sequence_seqrec, structure_seqrec):
        alignment_fasta_file = op.join(
            conf.CONFIGS['tcoffee_dir'],
            '{}-{}.fasta'.format(sequence_seqrec.id, structure_seqrec.id)
        )
        with open(alignment_fasta_file, 'w') as ofh:
            SeqIO.write([sequence_seqrec, structure_seqrec], ofh, 'fasta')
        tc = call_tcoffee.TCoffee(
            alignment_fasta_file, pdb_file=self.structure_file, mode='3dcoffee')
        alignment_output_file = tc.align()
        return alignment_output_file

    def _create_pir_alignment(self):
        pir_alignment_file = op.join(conf.CONFIGS['model_dir'], self.model_id + '.pir')
        with open(pir_alignment_file, 'w') as ofh:
            write_to_pir_alignment(
                ofh, 'sequence', self.sequence_id,
                '/'.join(str(seqrec.seq) for seqrec in self.sequence_seqrecords_aligned)
            )
            write_to_pir_alignment(
                ofh, 'structure', self.structure_id,
                '/'.join(str(seqrec.seq) for seqrec in self.structure_seqrecords_aligned)
            )
        return pir_alignment_file

    def _create_alignments_and_model(self):
        # Align sequence to structure.
        alignment_files = []
        domain_def_offsets = []
        model_domain_defs = []
        alignment_stats = []
        self.sequence_seqrecords_aligned, self.structure_seqrecords_aligned = [], []
        for sequence_seqrec, structure_seqrec in zip(
                self.sequence_seqrecords, self.structure_seqrecords):
            if str(sequence_seqrec.seq) != str(structure_seqrec.seq):
                # Sequence and structure are different, so perform alignment
                alignment_output_file = self._align_with_tcoffee(sequence_seqrec, structure_seqrec)
                alignment = AlignIO.read(alignment_output_file, 'fasta')
                assert len(alignment) == 2
                # Check to make sure that the sequence does not have very large overhangs
                # over the structure.
                # TODO: Do something similar for very large gaps
                # (long region of sequence without structure)
                # Right now Modeller will try to model those regions as loops (which end up looking
                # very unnatural.
                domain_def_offset = get_alignment_overhangs(alignment)
                if any(domain_def_offset):
                    logger.debug(
                        'Shortening uniprot domain sequence because the alignment had large '
                        'overhangs... (domain_def_offset: {})'.format(domain_def_offset)
                    )
                    cut_from_start = domain_def_offset[0] if domain_def_offset[0] else None
                    cut_from_end = -domain_def_offset[1] if domain_def_offset[1] else None
                    sequence_seqrec.seq = (
                        Seq(str(sequence_seqrec.seq)[cut_from_start:cut_from_end])
                    )
                    alignment_output_file = (
                        self._align_with_tcoffee(sequence_seqrec, structure_seqrec)
                    )
                    alignment = AlignIO.read(alignment_output_file, 'fasta')
                    assert len(alignment) == 2
                # Analyse the quality of the alignment
                alignment_identity, alignment_coverage, __, __ = analyze_alignment(alignment)
                alignment_score = score_alignment(alignment_identity, alignment_coverage)
                # Save results
                alignment_stats.append(
                    (alignment_identity, alignment_coverage, alignment_score)
                )
                alignment_files.append(alignment_output_file)
                self.sequence_seqrecords_aligned.append(alignment[0])
                self.structure_seqrecords_aligned.append(alignment[1])
            else:
                # Sequence and structure are the same; no need for alignment. Save dummy results.
                alignment_stats.append((1.0, 1.0, 1.0,))
                alignment_output_file = (
                    op.join(
                        conf.CONFIGS['model_dir'],
                        '{}-{}.aln'.format(sequence_seqrec.id, structure_seqrec.id))
                )
                with open(alignment_output_file, 'w') as ofh:
                    SeqIO.write([sequence_seqrec, structure_seqrec], ofh, 'clustal')
                alignment_files.append(alignment_output_file)
                domain_def_offset = (0, 0,)
                self.sequence_seqrecords_aligned.append(sequence_seqrec)
                self.structure_seqrecords_aligned.append(structure_seqrec)
            # either way
            domain_def_offsets.append(domain_def_offset)
            model_domain_def = (
                (domain_def_offset[0] + 1,
                 domain_def_offset[0] + 1 + len(sequence_seqrec), )
            )
            model_domain_defs.append(model_domain_def)

        # Add the HETATM chain if necesasry.
        assert len(self.sequence_seqrecords_aligned) == len(self.structure_seqrecords_aligned)
        # TODO: This looks wrong...
        if len(self.structure_seqrecords) == len(self.structure_seqrecords_aligned) + 1:
            self.sequence_seqrecords_aligned.append(self.structure_seqrecords[-1])
            self.structure_seqrecords_aligned.append(self.structure_seqrecords[-1])

        # Write *.pir alignment.
        self.pir_alignment_file = self._create_pir_alignment()
        logger.debug('Created pir alignment: {}'.format(self.pir_alignment_file))

        # Run modeller.
        self.modeller_results = run_modeller(
            self.pir_alignment_file, self.sequence_id, self.structure_id,
            new_chains=''.join(self.chain_ids)
        )

        # Save additional alignment info
        self.modeller_results['alignment_files'] = [
            op.relpath(f, conf.CONFIGS['unique_temp_dir'])
            for f in alignment_files]
        assert len(domain_def_offsets) <= 2
        self.modeller_results['domain_def_offsets'] = domain_def_offsets
        assert len(model_domain_defs) <= 2
        self.modeller_results['model_domain_defs'] = model_domain_defs
        self.modeller_results['alignment_stats'] = alignment_stats

    def _analyse_core(self):
        # Run the homology model through msms and get dataframes with all the
        # per atom and per residue SASA values
        analyze_structure = structure_analysis.AnalyzeStructure(
            op.join(conf.CONFIGS['unique_temp_dir'], self.modeller_results['model_file']),
            conf.CONFIGS['modeller_dir']
        )
        __, seasa_by_chain_separately, __, seasa_by_residue_separately = (
            analyze_structure.get_seasa()
        )

        # Get SASA only for amino acids in the chain of interest
        def _filter_df(df, chain_id, resname, resnum):
            df2 = df[
                (df['pdb_chain'] == chain_id) &
                (df['res_name'] == resname) &
                (df['res_num'] == resnum)
            ]
            return df2.iloc[0]['rel_sasa']

        self.relative_sasa_scores = {}
        for chain_id in self.modeller_chain_ids:
            self.relative_sasa_scores[chain_id] = []
            chain = self.modeller_structure[0][chain_id]
            for residue in chain:
                if residue.resname in structure_tools.AAA_DICT:
                    relative_sasa_score = _filter_df(
                        seasa_by_residue_separately,
                        chain_id=chain_id,
                        resname=residue.resname,
                        resnum=(str(residue.id[1]) + residue.id[2].strip())
                    )
                    self.relative_sasa_scores[chain_id].append(relative_sasa_score)
            number_of_aa = len(structure_tools.get_chain_sequence_and_numbering(chain)[0])
            if (number_of_aa != len(self.relative_sasa_scores[chain_id])):
                logger.error(
                    'Chain has {} non-hetatm AA, but we have SASA score for only {} AA.'
                    .format(number_of_aa, len(self.relative_sasa_scores[chain_id]))
                )
                raise errors.MSMSError()

    def _analyse_interface(self):
        # Get a dictionary of interacting residues
        interacting_residues = (
            structure_tools.get_interacting_residues(self.modeller_structure[0], r_cutoff=6.0)
        )
        _interacting_residues_complement = dict()
        for key, values in interacting_residues.items():
            for value in values:
                _interacting_residues_complement.setdefault(value, set()).add(key)
        interacting_residues.update(_interacting_residues_complement)

        # Get interacting residues (and interacting resnum) for chain 1 and 2
        def _get_a2b_contacts(a_idx, b_idx):
            # 2 if you want to get AA indexes (starting from 0)
            # 3 if you want to get AA residue numbering
            a2b_contacts = set()
            for key in interacting_residues:
                if key[0] == a_idx:
                    for value in interacting_residues[key]:
                        if value[0] == b_idx:
                            a2b_contacts.add(tuple(key[2:]))
            return a2b_contacts

        a2b_contacts = _get_a2b_contacts(0, 1)
        b2a_contacts = _get_a2b_contacts(1, 0)

        if not a2b_contacts or not b2a_contacts:
            logger.error('Chains are not interacting!')
            logger.error("interacting_residues: {}".format(interacting_residues))
            logger.error('a2b_contacts: {}'.format(a2b_contacts))
            logger.error('b2a_contacts: {}'.format(b2a_contacts))
            raise errors.ChainsNotInteractingError()

        def _validate_a2b_contacts(a2b_contacts, chain_idx):
            logger.debug('Validating chain {} interacting AA...'.format(chain_idx))
            interface_aa_a = ''.join([
                i[2] for i in a2b_contacts
            ])
            try:
                interface_aa_b = ''.join([
                    str(self.sequence_seqrecords[chain_idx].seq)[i[0]]
                    for i in a2b_contacts
                ])
            except IndexError as e:
                logger.error('{}: {}'.format(type(e), e))
                interface_aa_b = None

            logger.debug('interface_aa_a: {}'.format(interface_aa_a))
            logger.debug('interface_aa_b: {}'.format(interface_aa_b))
            logger.debug('a2b_contacts: {}'.format(a2b_contacts))
            logger.debug(
                'self.sequence_seqrecords[chain_idx].seq: {}'
                .format(self.sequence_seqrecords[chain_idx].seq)
            )
            logger.debug(
                "domain_def_offsets: {}".format(self.modeller_results['domain_def_offsets']))
            if interface_aa_a != interface_aa_b:
                raise errors.InterfaceMismatchError()

        _validate_a2b_contacts(a2b_contacts, 0)
        _validate_a2b_contacts(b2a_contacts, 1)

        # Using residue indexes
        # self.interacting_residues_x uses the POSITION of the residue
        # (e.g. [1,2,3] means the first three residues are interacting)
        self.interacting_aa_1 = sorted(i[0] + 1 for i in a2b_contacts)
        self.interacting_aa_2 = sorted(i[0] + 1 for i in b2a_contacts)

        # Interface area
        analyze_structure = structure_analysis.AnalyzeStructure(
            op.join(conf.CONFIGS['unique_temp_dir'], self.modeller_results['model_file']),
            conf.CONFIGS['modeller_dir']
        )
        (self.interface_area_hydrophobic,
         self.interface_area_hydrophilic,
         self.interface_area_total) = (
             analyze_structure.get_interface_area(self.modeller_chain_ids[:2])
        )

    def mutate(self, sequence_idx, mutation):
        """Introduce mutation into model.

        Parameters
        ----------
        sequence_idx : int
            Integer describing whether the mutation is on the first domain (`0`)
            or on the second domain (`1`).

        Raises
        ------
        MutationOutsideDomainError
        MutationOutsideInterfaceError
        """
        if (sequence_idx, mutation) in self.mutations:
            return self.mutations[(sequence_idx, mutation)]

        protein_id = self.sequence_seqrecords[sequence_idx].id
        chain_id = self.modeller_structure.child_list[0].child_list[sequence_idx].id

        mutation_errors = ''

        # Domain definitions, in case not the entire sequence was modelled
        domain_def_offset = self.modeller_results['domain_def_offsets'][sequence_idx]
        domain_def = (
            domain_def_offset[0],
            len(self.sequence_seqrecords[sequence_idx].seq) - domain_def_offset[1]
        )
        mutation_pos = int(mutation[1:-1])
        if mutation_pos > (domain_def[1] - domain_def[0] + 1):
            raise errors.MutationOutsideDomainError()

        position_modeller = (
            structure_tools.convert_position_to_resid(
                self.modeller_structure[0][chain_id],
                [mutation_pos])[0]
        )
        mutation_modeller = (mutation[0] + str(position_modeller) + mutation[-1])
        logger.debug('mutation: {}'.format(mutation))
        logger.debug('position_modeller: {}'.format(position_modeller))
        logger.debug('mutation_modeller: {}'.format(mutation_modeller))

        if len(self.sequence_seqrecords) == 1:
            partner_chain_idx = None
            partner_protein_id = ''
            partner_chain_id = None
        else:
            # TODO: slight hack getting partner_chain_idx
            partner_chain_idx = [
                i for i in range(len(self.sequence_seqrecords)) if i != sequence_idx
            ][0]
            partner_protein_id = self.sequence_seqrecords[partner_chain_idx].id
            partner_chain_id = (
                self.modeller_structure.child_list[0].child_list[partner_chain_idx].id
            )
            logger.debug('sequence_idx: {}'.format(sequence_idx))
            logger.debug('partner_chain_idx: {}'.format(partner_chain_idx))
            if sequence_idx == 0:
                logger.debug('interacting_aa_1: {}'.format(self.interacting_aa_1))
                if int(mutation[1:-1]) not in self.interacting_aa_1:
                    raise errors.MutationOutsideInterfaceError()
            elif sequence_idx == 1:
                logger.debug('interacting_aa_2: {}'.format(self.interacting_aa_2))
                if int(mutation[1:-1]) not in self.interacting_aa_2:
                    raise errors.MutationOutsideInterfaceError()
            else:
                logger.warning(
                    "Can't make sure that a mutation is inside an interface if there are only "
                    "two chains!"
                )

        mutation_id = '{}-{}-{}'.format(protein_id, partner_protein_id, mutation)

        if mutation_errors:
            results = dict(
                protein_id=protein_id,
                sequence_idx=sequence_idx,
                chain_modeller=chain_id,
                partner_chain_id=partner_chain_id,
                mutation_id=mutation_id,
                mutation_domain=mutation,
                mutation_errors=mutation_errors,
            )
            self.mutations[(sequence_idx, mutation)] = results
            return results

        # ...
        logger.debug('Running mutation with mutation_id: {}'.format(mutation_id))
        logger.debug('chain_id: {}'.format(chain_id))
        logger.debug('partner_chain_id: {}'.format(partner_chain_id))

        #######################################################################
        # Create a folder for all mutation data.
        mutation_dir = op.join(conf.CONFIGS['model_dir'], 'mutations', mutation_id)
        os.makedirs(mutation_dir, exist_ok=True)
        os.makedirs(mutation_dir, exist_ok=True)
        shutil.copy(op.join(conf.CONFIGS['data_dir'], 'rotabase.txt'), mutation_dir)

        #######################################################################
        # Copy the homology model to the mutation folder
        model_file = op.join(mutation_dir, op.basename(self.modeller_results['model_file']))
        shutil.copy(
            op.join(conf.CONFIGS['unique_temp_dir'], self.modeller_results['model_file']),
            model_file)

        #######################################################################
        # 2nd: use the 'Repair' feature of FoldX to optimise the structure
        fX = call_foldx.FoldX(model_file, chain_id, mutation_dir)
        repairedPDB_wt = fX('RepairPDB')

        #######################################################################
        # 3rd: introduce the mutation using FoldX
        mutCodes = [mutation_modeller[0] + chain_id + mutation_modeller[1:], ]
        logger.debug('Mutcodes for foldx: {}'.format(mutCodes))

        # Introduce the mutation using foldX
        fX_wt = call_foldx.FoldX(repairedPDB_wt, chain_id, mutation_dir)
        repairedPDB_wt_list, repairedPDB_mut_list = fX_wt('BuildModel', mutCodes)

        logger.debug('repairedPDB_wt_list: %s' % str(repairedPDB_wt_list))
        logger.debug('repairedPDB_mut_list: %s' % str(repairedPDB_mut_list))

        wt_chain_sequences = structure_tools.get_structure_sequences(repairedPDB_wt_list[0])
        mut_chain_sequences = structure_tools.get_structure_sequences(repairedPDB_mut_list[0])

        logger.debug('wt_chain_sequences: %s' % str(wt_chain_sequences))
        logger.debug('mut_chain_sequences: %s' % str(mut_chain_sequences))

        # Copy the foldX wildtype and mutant pdb files (use the first model if there are multiple)
        model_file_wt = op.join(mutation_dir, mutation_id + '-wt.pdb')
        model_file_mut = op.join(mutation_dir, mutation_id + '-mut.pdb')
        shutil.copy(repairedPDB_wt_list[0], model_file_wt)
        shutil.copy(repairedPDB_mut_list[0], model_file_mut)

        #######################################################################
        # 4th: set up the classes for the wildtype and the mutant structures
        fX_wt_list = list()
        for wPDB in repairedPDB_wt_list:
            fX_wt_list.append(call_foldx.FoldX(wPDB, chain_id, mutation_dir))

        fX_mut_list = list()
        for mPDB in repairedPDB_mut_list:
            fX_mut_list.append(call_foldx.FoldX(mPDB, chain_id, mutation_dir))

        #######################################################################
        # 5th: Calculate energies
        assert len(fX_wt_list) == 1
        stability_values_wt = ','.join(
            '{}'.format(f) for f in fX_wt_list[0]('Stability')
        )
        assert len(fX_wt_list) == 1
        stability_values_mut = ','.join(
            '{}'.format(f) for f in fX_mut_list[0]('Stability')
        )

        if len(self.sequence_seqrecords) == 1:
            complex_stability_values_wt = None
            complex_stability_values_mut = None
        else:
            assert len(fX_wt_list) == 1
            complex_stability_values_wt = ','.join(
                '{}'.format(f) for f in fX_wt_list[0]('AnalyseComplex')
            )
            assert len(fX_mut_list) == 1
            complex_stability_values_mut = ','.join(
                '{}'.format(f) for f in fX_mut_list[0]('AnalyseComplex')
            )

        #######################################################################
        # 6: Calculate all other relevant properties
        # (This also verifies that mutations match mutated residues in pdb structures).
        analyze_structure_wt = structure_analysis.AnalyzeStructure(
            repairedPDB_wt_list[0], mutation_dir,
        )
        analyze_structure_results_wt = analyze_structure_wt(
            chain_id, mutation_modeller, partner_chain_id)

        analyze_structure_mut = structure_analysis.AnalyzeStructure(
            repairedPDB_mut_list[0], mutation_dir,
        )
        analyze_structure_results_mut = analyze_structure_mut(
            chain_id, mutation_modeller, partner_chain_id)

        logger.debug('analyze_structure_results_wt: {}'.format(analyze_structure_results_wt))
        logger.debug('analyze_structure_results_mut: {}'.format(analyze_structure_results_mut))

        #######################################################################
        # 5th: calculate the energy for the wildtype
        results = dict(
            protein_id=protein_id,
            sequence_idx=sequence_idx,
            chain_modeller=chain_id,
            partner_chain_id=partner_chain_id,
            mutation_id=mutation_id,
            mutation_domain=mutation,
            mutation_errors=mutation_errors,
            #
            mutation_dir=mutation_dir,
            mutation_modeller=mutation_modeller,
            mutation_foldx=','.join(mutCodes),
            model_file_wt=model_file_wt,
            model_file_mut=model_file_mut,
            stability_energy_wt=stability_values_wt,
            stability_energy_mut=stability_values_mut,
            analyse_complex_energy_wt=complex_stability_values_wt,
            analyse_complex_energy_mut=complex_stability_values_mut,
        )
        for key, value in analyze_structure_results_wt.items():
            results[key + '_wt'] = value
        for key, value in analyze_structure_results_mut.items():
            results[key + '_mut'] = value

        # Another exit point
        self.mutations[(sequence_idx, mutation)] = results
        return results

    @property
    def result(self):
        result = dict(
            model_id=self.model_id,
            structure_file=op.relpath(self.structure_file, conf.CONFIGS['unique_temp_dir']),
            structure_id=self.structure_id,
            sequence_file=op.relpath(self.sequence_file, conf.CONFIGS['unique_temp_dir']),
            sequence_id=self.sequence_id,
            chain_ids=tuple(self.chain_ids),
            mutations=self.mutations,
            modeller_results_file=op.relpath(
                self.modeller_results_file, conf.CONFIGS['unique_temp_dir']),
            modeller_chain_ids=tuple(self.modeller_chain_ids),
            relative_sasa_scores=self.relative_sasa_scores,
            core_or_interface=self.core_or_interface,
        )
        # Dump modeller resutls
        for key, value in self.modeller_results.items():
            result[key] = value
        # For interfaces
        if len(self.sequence_seqrecords) > 1:
            result_interface = dict(
                interacting_aa_1=self.interacting_aa_1,
                interacting_aa_2=self.interacting_aa_2,
                interface_area_hydrophobic=self.interface_area_hydrophobic,
                interface_area_hydrophilic=self.interface_area_hydrophilic,
                interface_area_total=self.interface_area_total,
            )
            result.update(result_interface)
        return result


def perform_alignment(self, uniprot_seqrecord, pdb_seqrecord, mode, path_to_data):
    """
    """
    # Perform the alignment
    t_coffee_parameters = [
        uniprot_seqrecord,
        pdb_seqrecord,
        mode,
        logger,
    ]
    logger.debug(
        "Calling t_coffee with parameters:\n" +
        ', '.join(['{}'.format(x) for x in t_coffee_parameters]))
    tcoffee = call_tcoffee.tcoffee_alignment(*t_coffee_parameters)
    alignments = tcoffee.align()
    assert len(alignments) == 1
    alignment = alignments[0]

    # Save the alignment
    logger.debug(alignment)
    alignment_filename = alignment[0].id + '_' + alignment[1].id + '.aln'
    try:
        AlignIO.write(
            alignment, self.unique_temp_folder + 'tcoffee/' + alignment_filename, 'clustal')
    except IndexError as e:
        raise errors.EmptyPDBSequenceError('{}: {}'.format(type(e), e))
    temp_save_path = self.temp_archive_path + path_to_data
    subprocess.check_call("mkdir -p '{}'".format(temp_save_path), shell=True)
    subprocess.check_call("cp -f '{}' '{}'".format(
        self.unique_temp_folder + 'tcoffee/' + alignment_filename,
        temp_save_path + alignment_filename), shell=True)

    return alignment, alignment_filename


def analyze_alignment(alignment, pdb_contact_idxs=[]):
    """Return scores describing the qualit of the alignment.

    Returns
    -------
    identity : float <= 1
        Core identity.
    coverage : float <= 1
        Core coverage.
    if_identity : float <= 1
        Interface identity.
    if_coverage : float <= 1
        Interface coverage.
    """
    pdb_aa_idx = -1
    sequence_1_length = 0
    sequence_1_identity = 0
    sequence_1_coverage = 0
    interface_1_identity = 0
    interface_1_coverage = 0

    for aa_1, aa_2 in zip(*alignment):
        is_interface = False
        # Check if the amino acid falls in a gap
        if aa_1 == '-':
            continue
        sequence_1_length += 1
        # Check if the template is in a gap
        if aa_2 == '-':
            continue
        pdb_aa_idx += 1
        if pdb_aa_idx in pdb_contact_idxs:
            is_interface = True  # This is an interface amino acid
        sequence_1_coverage += 1  # Count as coverage
        if is_interface:
            interface_1_coverage += 1  # Count as coverage
        # Check if the template is identical
        if aa_1 != aa_2:
            continue
        sequence_1_identity += 1  # Count as identity
        if is_interface:
            interface_1_identity += 1  # Count as identity

    identity = sequence_1_identity / float(sequence_1_length)
    coverage = sequence_1_coverage / float(sequence_1_length)
    assert identity <= 1
    assert coverage <= 1

    if pdb_contact_idxs:
        if_identity = sequence_1_identity / float(len(pdb_contact_idxs))
        if_coverage = interface_1_coverage / float(len(pdb_contact_idxs))
        assert if_identity <= 1
        assert if_coverage <= 1
    else:
        if_identity = None
        if_coverage = None

    return identity, coverage, if_identity, if_coverage


def score_alignment(identity, coverage, alpha=0.95):
    """T-score from the interactome3d paper."""
    return alpha * (identity) * (coverage) + (1.0 - alpha) * (coverage)


def get_alignment_overhangs(alignment):
    """Remove gap overhangs from the alignments.

    There are cases where no template sequence is availible for a big chunk
    of the protein. Return the number of amino acids that should be removed
    from the start and end of the query sequence in order to match the template.
    """
    n_gaps_start = 0
    n_gaps_end = 0
    for aa_query, aa_template in zip(*alignment):
        if aa_query != '-' and aa_template == '-':
            n_gaps_start += 1
        else:
            break
    for aa_query, aa_template in reversed(list(zip(*alignment))):
        if aa_query != '-' and aa_template == '-':
            n_gaps_end += 1
        else:
            break
    return n_gaps_start, n_gaps_end


def write_to_pir_alignment(pir_alignment_filehandle, seq_type, seq_name, seq):
    """Write the `*.pir` alignment compatible with modeller.

    Parameters
    ----------
    seq_type : str
        One of: ['sequence', 'structure'], in that order.
    seq_name : str
        Name to appear in the alignment.
    seq : str
        Alignment sequence.
    """
    pir_alignment_filehandle.write('>P1;' + seq_name + '\n')
    pir_alignment_filehandle.write(seq_type + ':' + seq_name + ':.:.:.:.::::\n')
    pir_alignment_filehandle.write(seq + '*')
    pir_alignment_filehandle.write('\n\n')


def run_modeller(
        pir_alignment_file, target_id, template_id, new_chains='ABCDEFGHIJKLMNOPQRSTUVWXYZ'):
    """
    """
    logger.debug(
        "Calling modeller with parameters:\n" +
        'pir_alignment_file: {}\n'.format([pir_alignment_file]) +
        'target_id: {}\n'.format(target_id) +
        'template_id: {}\n'.format(template_id)
    )
    modeller = call_modeller.Modeller(
        [pir_alignment_file], target_id, template_id, conf.CONFIGS['unique_temp_dir'])

    with switch_paths(conf.CONFIGS['modeller_dir']):
        norm_dope, pdb_filename = modeller.run()

    raw_model_file = op.join(conf.CONFIGS['modeller_dir'], pdb_filename)

    # If there is only one chain in the pdb, label that chain 'A'
    io = PDBIO()
    structure = structure_tools.get_pdb_structure(raw_model_file)
    chains = structure[0].child_list
    logger.debug('Modeller chain ids: ' + ', '.join(chain.id for chain in chains))
    for i in range(len(chains)):
        chains[i].id = new_chains[i]
    logger.debug('Corrected chain ids: ' + ', '.join(chain.id for chain in chains))
    io.set_structure(structure)
    model_file = op.splitext(pir_alignment_file)[0] + '.pdb'
    io.save(model_file)

    results = {
        'model_file': op.relpath(model_file, conf.CONFIGS['unique_temp_dir']),
        'raw_model_file': op.relpath(raw_model_file, conf.CONFIGS['unique_temp_dir']),
        'norm_dope': norm_dope,
        'pir_alignment_file': op.relpath(pir_alignment_file, conf.CONFIGS['unique_temp_dir']),
    }
    return results
