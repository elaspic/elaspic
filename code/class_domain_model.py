# -*- coding: utf-8 -*-

import os
import subprocess
import tarfile
import json

from Bio import SeqIO, AlignIO

from class_pdbTemplate import pdbTemplate
from class_time import runTime as rt
from class_interfaceSize import interfaceSize, getDSSP
from class_pysicoChemicalProperties import pysiChem
from class_foldX import foldX
from HelperFunction_prepareModeller import prepareModeller
from modeller import ModellerError

import class_error as error
import class_modeller as mod
import class_sql as sql



class GetModel(object):
    """
    """
    
    def __init__(self, tmpPath, unique, pdbPath, savePDB, saveAlignments,
                  pool, semaphore, db, log, path_to_archive, get_template):
        """
        input:
        tmpPath             type 'str'
        unique              type 'str'
        pdbPath             type 'str'
        savePDB             type 'str'
        saveAlignments      type 'str'
        pool                type class '__main__.ActivePool'
        semaphore           type class 'multiprocessing.synchronize.Semaphore'
        """
        self.tmpPath = tmpPath
        self.unique = unique + '/'
        self.pdbPath = pdbPath
        self.savePDB = savePDB
        self.saveAlignments = saveAlignments
        self.pool = pool
        self.semaphore = semaphore
        self.db = db
        
        # get the logger from the parent and add a handler
        self.log = log
        self.path_to_archive = path_to_archive
        self.get_template = get_template
    
    
    def __call__(self, uniprot_domain, uniprot_template):
        
        # Folder for storing files for export to output
        save_path = self.tmpPath + self.unique + '/' + uniprot_domain.path_to_data + '/'
        
        if type(uniprot_domain) == sql.UniprotDomain:
            
            # Get the canonical uniprot sequence
            uniprot_sequence = self.db.get_uniprot_sequence(uniprot_domain.uniprot_id)
                
            # Cut it to domain boundaries
            uniprot_sequence_domain = self.make_SeqRec_object(
                                    uniprot_sequence, 
                                    sql.decode_domain(uniprot_template.domain_def),
                                    uniprot_domain.uniprot_id)
            
            # Load previously-calculated alignments
            alignment, __ = self.get_template.get_alignment(uniprot_template)
            
            # Put the data into lists so that it's independent of whether we
            # are looking at a core or interface mutation
            chains_pdb = [uniprot_template.domain.pdb_chain, ]
            pdb_domain_definitions = [sql.decode_domain(uniprot_template.domain.pdb_domain_def), ]
            sequences = [uniprot_sequence_domain, ]
            alignments = [alignment, ]
            
        
        if type(uniprot_domain) == sql.UniprotDomainPair:
            
            # Get sequences of the first and second protein          
            uniprot_sequence_1 = self.db.get_uniprot_sequence(uniprot_template.domain_1.uniprot_id)
            uniprot_sequence_2 = self.db.get_uniprot_sequence(uniprot_template.domain_2.uniprot_id)
                
            # Cut sequence to boundaries and set sequence ID            
            uniprot_sequence_1_domain = self.make_SeqRec_object(
                                    uniprot_sequence_1, 
                                    sql.decode_domain(uniprot_template.domain_def_1),
                                    uniprot_domain.uniprot_domain_1.uniprot_id)
                
            uniprot_sequence_2_domain = self.make_SeqRec_object(
                                    uniprot_sequence_2, 
                                    sql.decode_domain(uniprot_template.domain_def_2), 
                                    uniprot_domain.uniprot_domain_2.uniprot_id)
                
            # Load previously-calculated alignments
            alignment_1, alignment_2 = self.get_template.get_alignment(uniprot_template)
            
            chains_pdb = [uniprot_template.domain_1.pdb_chain, uniprot_template.domain_2.pdb_chain]
            pdb_domain_definitions = [sql.decode_domain(uniprot_template.domain_1.pdb_domain_def),
                                           sql.decode_domain(uniprot_template.domain_2.pdb_domain_def)]
            sequences = [uniprot_sequence_1_domain, uniprot_sequence_2_domain]
            alignments = [alignment_1, alignment_2]


        # Common -----------------
        # save_path is where the pdb sequences and the pdb with the required chains are saved
        sequences, alignments, chains_modeller, switch_chain, HETflag, HETATMsInChain_SEQnumbering = \
            self.prepareInput(uniprot_template.pdb_id, chains_pdb, pdb_domain_definitions, sequences, alignments, save_path)
        
        # copy the chain IDs... they are used afterwards for renaming and copy 
        # is needed so that they are not overwritten
        target_ids = []
        template_ids = []
        for i in range(0,len(alignments)):
            target_ids.append(alignments[i][0].id)
            template_ids.append(alignments[i][1].id)
        
        normDOPE_wt, pdbFile_wt = self._getModel(alignments, target_ids, template_ids, HETATMsInChain_SEQnumbering, chains_modeller, save_path)
        
        # Save another copy of the alignment in the export folder
        modeller_path = self.tmpPath + self.unique + '/modeller/'
        archive_save_path = self.path_to_archive + uniprot_domain.path_to_data + '/'
        subprocess.check_call('cp ' + modeller_path + pdbFile_wt + ' ' + save_path + pdbFile_wt, shell=True)
        subprocess.check_call('cp ' + modeller_path + pdbFile_wt + ' ' + archive_save_path + pdbFile_wt, shell=True)

        # Rename the wildtype pdb file
        pdbFile_wt_renamed = self.uniprot_id + '_' + self.mutations + '.pdb'
        system_command = 'mv ' + modeller_path + pdbFile_wt + ' ' + modeller_path + pdbFile_wt_renamed
        subprocess.check_call(system_command, shell=True)

        # Save the data to mutation row objects
        if type(uniprot_domain) == sql.UniprotDomain:
            uniprot_model = sql.UniprotDomainModel()
            uniprot_model.uniprot_domain_id = uniprot_template.uniprot_domain_id
        
            uniprot_model.model_filename = pdbFile_wt
            uniprot_model.chain = chains_modeller[0]
        
            uniprot_model.norm_dope = normDOPE_wt
            uniprot_model.het_flag = HETflag
            uniprot_model.switch_chain = switch_chain
            
            uniprot_model.interface_aa = ''

        if type(uniprot_domain) == sql.UniprotDomainPair:
            uniprot_model = sql.UniprotDomainPairModel()
            uniprot_model.uniprot_domain_pair_id = uniprot_template.uniprot_domain_pair_id
        
            uniprot_model.model_filename = pdbFile_wt
            uniprot_model.norm_dope = normDOPE_wt
            uniprot_model.het_flag = HETflag
            uniprot_model.switch_chain = switch_chain
            
            uniprot_model.chain_1 = chains_modeller[0]
            uniprot_model.chain_2 = chains_modeller[1]
            uniprot_model.interface_aa_1 = ''
            uniprot_model.interface_aa_2 = ''
        
        # Save the model data to the output folder
        with open(archive_save_path + 'model.json', 'wb') as fh:
            json.dump(sql.row2dict(uniprot_model), fh, indent=4, separators=(',', ': ')) # separators to avoid tailing whitespace        
        
        return uniprot_model
        

    def prepareInput(self, pdbCode, chains, domains, sequences, alignments, savePDB):

        pdb = pdbTemplate(self.pdbPath, pdbCode, chains, domains, savePDB)

        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract() # dict with chain ID as key and Nr. of HETATM as value

        SWITCH_CHAIN = False
        if chains != chains_pdb_order:
            SWITCH_CHAIN = True
            sequence_new_order = list()
            alignments_new_order = list()
            for item in chains_pdb_order:
                i = chains.index(item)
                sequence_new_order.append(sequences[i])
                alignments_new_order.append(alignments[i])
            sequences = sequence_new_order
            alignments = alignments_new_order
            
            chains = chains_pdb_order
        
        # to get the range of the domain make sure to call extract before getChainNumberingNOHETATMS!!
        chainNumberingDomain = dict()

        for i in range(0,len(chains)):
            chainNumberingDomain[chains[i]] = pdb.getChainNumbering(chains[i])

        HETATMsInChain_SEQnumbering = [ dict() for x in range(0,len(chains)) ]
        
        for i in range(0,len(chains)):
            HETATMsInChain_SEQnumbering[i][chains[i]] = list()
            for item in HETATMsInChain_PDBnumbering[chains[i]]:
                    try:
                        HETATMsInChain_SEQnumbering[i][chains[i]].append(chainNumberingDomain[chains[i]].index(item))
                    except ValueError:
                        HETATMsInChain_SEQnumbering[i][chains[i]].append(10001) # add one add the end that will be outside of the sequence

        return sequences, alignments, chains, SWITCH_CHAIN, HETflag, HETATMsInChain_SEQnumbering
        
        
    def _getModel(self, alignments, target_ids, template_ids,
                       HETATMsInChain_SEQnumbering,
                       chains,
                       savePDB
                       ):
        outFile = self.tmpPath + self.unique + '/outFile_wildtype'
        
        # generate the input for modeller from the above generated alignment
        prepareModeller(outFile, 
                        alignments, 
                        target_ids, 
                        template_ids, 
                        HETATMsInChain_SEQnumbering,
                        chains
                        )
        
        inFile = outFile

        modeller_target_id = '_'.join(target_ids)
        
        if len(template_ids) == 1:
            modeller_template_id = template_ids[0]
        else:
            modeller_template_id = template_ids[0] + template_ids[1][-1] # if more than two chains are used this has to be adjusted

        modeller_path = self.tmpPath + self.unique + '/modeller/'
        os.chdir(modeller_path) # from os

        modeller = mod.modeller([inFile], 
                                modeller_target_id, 
                                modeller_template_id, 
                                savePDB, 
                                self.tmpPath + self.unique + '/', 
                                self.modeller_runs,
                                loopRefinement=True)
        normDOPE, pdbFile = modeller.run()

        os.chdir(self.PWD) # go back to the working directory
        
        return normDOPE, pdbFile        