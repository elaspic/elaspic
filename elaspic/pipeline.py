# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import object

import os
import os.path as op
import json
import signal
import six
import logging 
from functools import wraps

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import conf, errors, helper, structure_tools, structure_analysis, sequence, model, predictor

logger = logging.getLogger(__name__)
configs = conf.Configs()

sql_db = None
domain_alignment = None
domain_model = None
domain_mutation = None

MAX_DISTANCE_BETWEEN_INTERACTING_CHAINS = 6 # Angstrom
ELASPIC_LOGO = """

8888888888 888             d8888  .d8888b.  8888888b. 8888888 .d8888b.  
888        888            d88888 d88P  Y88b 888   Y88b  888  d88P  Y88b 
888        888           d88P888 Y88b.      888    888  888  888    888 
8888888    888          d88P 888  "Y888b.   888   d88P  888  888        
888        888         d88P  888     "Y88b. 8888888P"   888  888        
888        888        d88P   888       "888 888         888  888    888 
888        888       d8888888888 Y88b  d88P 888         888  Y88b  d88P 
8888888888 88888888 d88P     888  "Y8888P"  888       8888888 "Y8888P"  

"""

#%% 
class Pipeline(object):

    def __init__(self, configurations):
        """
        It should be possible to initialize one pipeline and call it in parallel using different
        mutations as input
        """       
        # Read the configuration file and set the variables
        if isinstance(configurations, six.string_types):
            conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            configs.update(**configurations)
        
        # TODO: remove so error message does not appear in a production release.  
        self._validate_temp_path(conf.configs)

        # Initialize a logger
        for line in ELASPIC_LOGO.split('\n'):
            logger.info(line)

        self.PWD = os.getcwd()
        
        self.sequences = {}
        self.models = {}
        self.predictions = {}



    def _validate_temp_path(self, configs):
        """
        Make sure that we are using a job specific temporary folder if we are on a cluster.
        """
        hostname = helper.get_hostname()
        no_job_specific_folder = configs['temp_dir'].startswith(configs['global_temp_dir'])
        on_node_with_manditory_job_specific_folder = (
            any([(x.lower() in hostname) for x in ['node', 'behemoth', 'grendel', 'beagle']])
        )
        if no_job_specific_folder and on_node_with_manditory_job_specific_folder:
            raise Exception('You should be using a temp folder that it specific to the particular job!')
    

    def run(self):
        raise NotImplementedError()
        
        
    def run_sequence(self, sequence):
        raise NotImplementedError()


    def _get_model(self):
        raise NotImplementedError()
    
    
    def run_model(self, model):
        raise NotImplementedError()


    def _get_mutation(self, model, mutation_data):
        raise NotImplementedError


    def run_mutation(self, mutation):
        raise NotImplementedError()


## Header
#pdb_chain = 'A'
#mutation = 'M1P'
#protein_id = pdb_id + pdb_chain
#chain_sequence, chain_numbering = sp.sp.get_chain_sequence_and_numbering(pdb_chain)
#sequence_file = op.join(working_dir, protein_id + '.fasta')
#structure_file = op.join(working_dir, protein_id + '.pdb')
#
#seqrecord = SeqRecord(id=protein_id, seq=Seq(chain_sequence))
#with open(op.join(working_dir, protein_id + '.fasta'), 'w') as ofh:
#    SeqIO.write(seqrecord, ofh, 'fasta')
#
#my_sequence = sequence.Sequence(sequence_file)
#assert op.isfile(op.join(working_dir, 'sequence', protein_id + '_provean_supset'))
#assert op.isfile(op.join(working_dir, 'sequence', protein_id + '_provean_supset.fasta'))
#
#my_sequence.mutate(mutation)
#
#
#
##%% Model
#
## Header
#pdb_chain = 'A'
#mutation = 'M1P'
#protein_id = pdb_id + pdb_chain
#chain_sequence, chain_numbering = sp.sp.get_chain_sequence_and_numbering(pdb_chain)
#sequence_file = op.join(working_dir, protein_id + '.fasta')
#structure_file = op.join(working_dir, protein_id + '.pdb')
#
#my_model = model.Model(sequence_file, structure_file)
#
#my_model.mutate(1, mutation)
#
#
#
##%%
#from importlib import reload
#reload(conf)
#reload(pipeline)
#reload(structure_tools)
#reload(structure_analysis)
#reload(sequence)
#reload(model)
#reload(predictor)
#
#sequence_data = my_sequence.__dict__
#structure_data = my_structure.__dict__
#
#my_predictor = predictor.Predictor()



#%%
class DatabasePipeline(Pipeline):
    
    def __init__(self, protein_id, configurations=None):
        """
        """
        super().__init__(configurations)




    def make_sequence(self, sequence_idx):
        ...
        
        
    def make_model(self):
        ...
        
        
    def mutate(self):
        ...






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
            self.seqrecords = list(SeqIO.parse(self.sequence_file, 'fasta'))
            for i, seqrec in enumerate(self.seqrecords):
                seqrec.id = '{}_{}'.format(seqrec.id, str(i))
        else:
            self.seqrecords = [
                SeqRecord(
                    id='{}{}_{}'.format(self.pdb_id, chain_id, i), 
                    seq=Seq(self.sp.chain_sequence_dict[chain_id]))
                for (i, chain_id) in enumerate(self.sp.chain_ids)
            ]
        

    def run(self, chain, mutation):
        """Run all steps of the pipeline.
        
        Made this a separate method because it makes it easier to test
        :py:func:`compute_provean`, :py:func:`compute_model`, and :py:func:`compute_mutation`.
        """
        ...
    
    
    def precalculate_all(self):
        """
        """
        for chain_id in self.sp.chain_ids:
            if chain_id == self.sp.hetatm_chain_id:
                continue
            chain_pos = self._get_chain_pos(chain_id)
            self.make_sequence(chain_pos)
            self.make_model(chain_pos)
        for chain_idxs in self.sp.interacting_chain_idxs:
            self.make_model(chain_idxs)


    def _get_chain_pos(self, chain_id):
        chain_idx = [i for (i, chain) in enumerate(self.sp.structure[0].child_list) if chain.id == chain_id]
        if len(chain_idx) == 0:
            raise errors.PDBChainError(
                'Chain {} was not found in PDB {}!'.format(chain_id, self.sp.pdb_file))
        elif len(chain_idx) > 1:
            raise errors.PDBChainError(
                'Chain {} was found more than once in PDB {}!'.format(chain_id, self.sp.pdb_file))
        return chain_idx[0]
                

    def make_sequence(self, position):
        """
        Raises
        -------
        errors.ProveanError
        errors.ProveanResourceError
        """
        if position in self.sequences:
            return self.sequences[position]
        
        seqrecord = self.seqrecords[position]
        sequence_file = op.join(
            configs['sequence_dir'],
            seqrecord.id + '.fasta'
        )
        with open(sequence_file, 'w') as ofh:
            SeqIO.write(seqrecord, ofh, 'fasta')
        self.sequences[position] = sequence.Sequence(sequence_file)
        
        
    def make_model(self, sequences_to_model):
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
        if not hasattr(sequences_to_model, '__getitem__'):
            positions = (sequences_to_model,)
        else:
            positions = tuple(sorted(sequences_to_model))
            
        if positions in self.models:
            return self.models[positions]
        
        # Target sequence file
        seqrecords = [self.seqrecords[pos] for pos in positions]
        sequence_file = op.join(
            configs['model_dir'],        
           '_'.join(seqrec.id for seqrec in seqrecords) + '.fasta'
        )
        with open(sequence_file, 'w') as ofh:
            SeqIO.write(seqrecords, ofh, 'fasta')

        # Template structure file
        chain_string = ''.join(self.sp.structure[0].child_list[pos].id for pos in positions)
        structure_file = op.join(configs['unique_temp_dir'], self.pdb_id + chain_string + '.pdb')
        assert op.isfile(structure_file)
        
        self.models[positions] = model.Model(sequence_file, structure_file)



    def mutate_all(self, sequence_to_mutate, mutation):
        """
        Wrapper around `mutate`.            
        """
        for chain_idxs in self.sp.interacting_chain_idxs:
            if sequence_to_mutate in chain_idxs:
                self.mutate(chain_idxs, sequence_to_mutate, mutation)
        
        

    def mutate(self, sequences_to_model, sequence_to_mutate, mutation):
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
        if not hasattr(sequences_to_model, '__getitem__'):
            sequences_to_model = (sequences_to_model,)
        else:
            sequences_to_model = tuple(sorted(sequences_to_model))
        assert sequence_to_mutate in sequences_to_model
        
        if sequence_to_mutate not in self.sequences:
            self.make_sequence(sequence_to_mutate)
            assert sequence_to_mutate in self.sequences
            
        if sequences_to_model not in self.models:
            self.make_model(sequences_to_model)
            assert sequences_to_model in self.models
        
        self.sequences[sequence_to_mutate].mutate(mutation)
        
        self.models[sequences_to_model].mutate(sequences_to_model.index(sequence_to_mutate), mutation)



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
