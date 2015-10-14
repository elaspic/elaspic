# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()

import os.path as op
import logging 
from functools import wraps

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

    def __init__(self, structure_file, sequence_file=None, configurations=None):
        """
        TODO: Add an option to store provean results based on sequence hash.
        """
        super().__init__(configurations)
        
        # Input parameters
        self.pdb_id = op.splitext(op.basename(structure_file))[0]
        self.pdb_file = structure_file
        self.sequence_file = sequence_file

        logger.info('pdb_file: {}'.format(self.pdb_file))
        logger.info('pwd: {}'.format(self.PWD))
        
        ### Load PDB structure and extract required sequences and chains.
        #fix_pdb(self.pdb_file, self.pdb_file)
        self.sp = structure_tools.StructureParser(self.pdb_file)
        self.sp.extract()
        self.sp.save_structure(configs['unique_temp_dir'])
        self.sp.save_sequences(configs['unique_temp_dir'])
            
        if self.sequence_file is not None:
            # Read template sequences from the sequence file
            self.seqrecords = list(SeqIO.parse(self.sequence_file, 'fasta'))
            for i, seqrec in enumerate(self.seqrecords):
                seqrec.id = '{}_{}'.format(seqrec.id, str(i))
        else:
            # Read template sequences from the PDB
            self.seqrecords = [
                SeqRecord(
                    id='{}{}_{}'.format(self.pdb_id, chain_id, i), 
                    seq=Seq(self.sp.chain_sequence_dict[chain_id]))
                for (i, chain_id) in enumerate(self.sp.chain_ids)
            ]
            
    
    def precalculate_all(self):
        """
        """
        for chain_id in self.sp.chain_ids:
            if chain_id == self.sp.hetatm_chain_id:
                continue
            chain_pos = self._get_chain_pos(chain_id)
            Factory.get('sequence', self.seqrecords, chain_pos, None)
            Factory.get('sequence', self.seqrecords, self.sp, (chain_pos,))
        for chain_idxs in self.sp.interacting_chain_idxs:
            chain_idxs = self._get_position_tuple(chain_idxs)
            Factory.get('sequence', self.seqrecords, self.sp, chain_idxs)


    def mutate_all(self, mutation_idx, mutation):
        """       
        """
        self.get_mutation_score(mutation_idx, mutation_idx, mutation)
        for idxs in self.sp.interacting_chain_idxs:
            if mutation_idx in idxs:
                self.get_mutation_score(idxs, mutation_idx, mutation)
        

    def get_sequence(self, idx):
        return Factory.get('sequence', self.seqrecords, idx, None)


    def get_model(self, idxs):
        idxs = self._get_position_tuple(idxs)
        return Factory.get('model', self.seqrecords, self.sp, idxs)


    def get_mutation_score(self, idxs, mutation_idx, mutation):
        sequence = self.get_sequence(mutation_idx)
        model = self.get_model(idxs)
        return Factory.get('mutation', sequence, model, idxs.index(mutation_idx), mutation)


    def _get_chain_pos(self, chain_id):
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
        
    
    def _get_position_tuple(self, positions):
        """Sort positions and return as as tuple.
        """
        if not hasattr(positions, '__getitem__'):
            positions = (positions,)
        else:
            positions = tuple(sorted(positions))
            


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
        
    def __exit__(self):
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
        
        ...
    
    def __bool__(self):
        return True
    
    def __enter__(self):
        sequence.mutate(self.mutation)
        model.mutate(self.mutation_idx, self.mutation)
        self.feature_dict = dict(
            
        )        
        
    def run(self):
        self.sequences[sequence_to_mutate].mutate(mutation)
        
        self.models[sequences_to_model].mutate(sequences_to_model.index(sequence_to_mutate), mutation)
        
    def __exit__(self):
        ...
        






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


     

#%% PDB only
#class UserStructureMutation(Base):
#    __tablename__ = 'user_structure_mutation'
#    _indexes = [
#        ['pdb_id', 'pdb_chain', 'pdb_chain_2', 'pdb_mutation'],
#    ]
#    __table_args__ = get_table_args(__tablename__, _indexes, ['schema_version_tuple'])
#    
#    id = sa.Column(sa.Integer, nullable=False, primary_key=True, autoincrement=True)
#    
#    ### Input
#    unique = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')), 
#        primary_key=True)
#    pdb_file = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    pdb_id = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    pdb_chains = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    
#    
#    ### Mutation    
#    # Original PDB coordinates
#    pdb_mutation = sa.Column(sa.String(SHORT), nullable=False)
#    # Where 1 is the start of the domain (i.e. chain)
#    domain_mutation = sa.Column(sa.String(SHORT), nullable=False)
#    # Where 1 is the start of the protein (i.e. chain A residue 1)
#    modeller_mutation = sa.Column(sa.String(SHORT), nullable=False)
#    
#    
#    ### Homology Model
#    # Core / Interface
#    model_filename = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    norm_dope = sa.Column(sa.Float)
#    sasa_score = sa.Column(sa.Text)
#        
#    # Interface
#    pdb_chain = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    pdb_chain_2 = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')))
#    interface_area_hydrophobic = sa.Column(sa.Float)
#    interface_area_hydrophilic = sa.Column(sa.Float)
#    interface_area_total = sa.Column(sa.Float)
#    interacting_aa_1 = sa.Column(sa.Text)
#    interacting_aa_2 = sa.Column(sa.Text)
#    
#    
#    
#    ### Mutation
#    mutation = sa.Column(sa.String(SHORT), nullable=False)
#    mutation_errors = sa.Column(sa.Text)
#    model_filename_wt = sa.Column(sa.String(MEDIUM))
#    model_filename_mut = sa.Column(sa.String(MEDIUM))
#    chain_modeller = sa.Column(sa.String(SHORT))
#    mutation_modeller = sa.Column(sa.String(SHORT))
#    analyse_complex_energy_wt = sa.Column(sa.Text)
#    stability_energy_wt = sa.Column(sa.Text)
#    analyse_complex_energy_mut = sa.Column(sa.Text)
#    stability_energy_mut = sa.Column(sa.Text)
#    physchem_wt = sa.Column(sa.Text)
#    physchem_wt_ownchain = sa.Column(sa.Text)
#    physchem_mut = sa.Column(sa.Text)
#    physchem_mut_ownchain = sa.Column(sa.Text)
#    matrix_score = sa.Column(sa.Float)
#    secondary_structure_wt = sa.Column(sa.Text)
#    solvent_accessibility_wt = sa.Column(sa.Float)
#    secondary_structure_mut = sa.Column(sa.Text)
#    solvent_accessibility_mut = sa.Column(sa.Float)
#    contact_distance_wt = sa.Column(sa.Float)
#    contact_distance_mut = sa.Column(sa.Float)
#    provean_score = sa.Column(sa.Float)
#    ddg = sa.Column(sa.Float, index=False)
#    mut_date_modified = sa.Column(
#        sa.DateTime, default=datetime.datetime.utcnow,
#        onupdate=datetime.datetime.utcnow, nullable=False)
#        
#    
#    ### Relationships
#    model = sa.orm.relationship(
#        UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
#        backref=sa.orm.backref('mutations', cascade='expunge')) # many to one
#
