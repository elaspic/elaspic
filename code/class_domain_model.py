# -*- coding: utf-8 -*-

import os
import subprocess

from class_pdbTemplate import pdbTemplate
from HelperFunction_prepareModeller import prepareModeller
import class_error as error
import class_modeller as mod
import class_sql as sql
import class_analyze_structure


class GetModel(object):
    """
    """
    
    def __init__(self, tmpPath, unique, pdbPath, db, log, modeller_runs, buildModel_runs, PWD):
        """
        Produces a model of a single uniprot domain or a domain pair
        
        Parameters
        ----------
        tmpPath : 'str'
        unique : 'str'
        pdbPath : 'str'
        savePDB : 'str'
        
        Returns
        -------
        model : sql.object
            Object having all the information about the produced model. The model
            itself is stored in a file on disk.
        """
        self.tmpPath = tmpPath
        self.unique = unique + '/'
        self.pdbPath = pdbPath
        self.db = db
        
        # get the logger from the parent and add a handler
        self.log = log
        self.modeller_runs = modeller_runs
        self.buildModel_runs = buildModel_runs
        self.PWD = PWD

    
    def __call__(self, uniprot_domain, uniprot_template):
        
        # Folder for storing files for export to output
        save_path = self.tmpPath + str(uniprot_domain.path_to_data)
        
        if type(uniprot_domain) == sql.UniprotDomain:
            # Get the canonical uniprot sequence
            uniprot_sequence = self.db.get_uniprot_sequence(uniprot_domain.uniprot_id)
                
            # Cut it to domain boundaries
            uniprot_sequence_domain = uniprot_sequence[
                sql.decode_domain(uniprot_template.domain_def)[0]-1:
                sql.decode_domain(uniprot_template.domain_def)[1]]
            
            # Load previously-calculated alignments
            alignment, __ = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)
            
            # Put the data into lists so that it's independent of whether we
            # are looking at a core or interface mutation
            pdb_id = uniprot_template.domain.pdb_id
            chains_pdb = [uniprot_template.domain.pdb_chain, ]
            pdb_domain_definitions = [sql.decode_domain(uniprot_template.domain.pdb_domain_def), ]
            sequences = [uniprot_sequence_domain, ]
            alignments = [alignment, ]
            
        
        if type(uniprot_domain) == sql.UniprotDomainPair:
            # Get sequences of the first and second protein          
            uniprot_sequence_1 = self.db.get_uniprot_sequence(uniprot_domain.uniprot_domain_1.uniprot_id)
            uniprot_sequence_2 = self.db.get_uniprot_sequence(uniprot_domain.uniprot_domain_2.uniprot_id)
                
            # Cut sequence to boundaries and set sequence ID            
            uniprot_sequence_1_domain = uniprot_sequence_1[
                sql.decode_domain(uniprot_template.domain_def_1)[0]-1:
                sql.decode_domain(uniprot_template.domain_def_1)[1]]
                
            uniprot_sequence_2_domain = uniprot_sequence_2[
                sql.decode_domain(uniprot_template.domain_def_2)[0]-1:
                sql.decode_domain(uniprot_template.domain_def_2)[1]]
                
            # Load previously-calculated alignments
            alignment_1, alignment_2 = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)
            
            pdb_id = uniprot_template.domain_1.pdb_id
            chains_pdb = [uniprot_template.domain_1.pdb_chain, uniprot_template.domain_2.pdb_chain]
            pdb_domain_definitions = [sql.decode_domain(uniprot_template.domain_1.pdb_domain_def),
                                           sql.decode_domain(uniprot_template.domain_2.pdb_domain_def)]
            sequences = [uniprot_sequence_1_domain, uniprot_sequence_2_domain]
            alignments = [alignment_1, alignment_2]

        model_errors = ''
        
        for a, b in zip(alignments, chains_pdb):
            self.log.debug('Alignment: ')
            self.log.debug(a)
            self.log.debug('Chain: %s' % b)
        
        # copy the chain IDs... they are used afterwards for renaming and copy 
        # is needed so that they are not overwritten
        target_ids = []
        template_ids = []
        for i in range(0,len(alignments)):
            target_ids.append(alignments[i][0].id)
            template_ids.append(alignments[i][1].id)


        chains = [template_id[-1] for template_id in template_ids]
        # Common---------------------------------------------------------------
        # save_path is where the pdb sequences and the pdb with the required chains are saved
        self.log.debug('Prepare input...')
        sequences, alignments, chains_modeller, switch_chain, HETflag, HETATMsInChain_SEQnumbering = \
            self.prepareInput(pdb_id, chains, pdb_domain_definitions, sequences, alignments, save_path)
        
        normDOPE_wt, pdbFile_wt, knotted = self._getModel(alignments, target_ids, template_ids, HETATMsInChain_SEQnumbering, chains_modeller, save_path)
        if knotted:
            model_errors += 'knotted, '
        self.log.debug('Model pdb file: %s, knotted: %s' % (pdbFile_wt, knotted,))
                
        # Save another copy of the modelled structure in the tmp export folder
        modeller_path = self.tmpPath + self.unique + 'modeller/'
        subprocess.check_call('cp ' + modeller_path + pdbFile_wt + ' ' + save_path + pdbFile_wt, shell=True)
        
        if type(uniprot_domain) == sql.UniprotDomain:
            uniprot_model = sql.UniprotDomainModel()
            uniprot_model.uniprot_domain_id = uniprot_template.uniprot_domain_id
        
            uniprot_model.model_filename = pdbFile_wt
            uniprot_model.chain = chains_modeller[0]
            uniprot_model.het_flag = HETflag[chains_modeller[0]]
            uniprot_model.switch_chain = switch_chain
            
            uniprot_model.norm_dope = normDOPE_wt
            
            # Get SASA using pops
            analyze_structure = class_analyze_structure.AnalyzeStructure(modeller_path, pdbFile_wt, ['A'], self.log, domain_defs=None)
            try:
                self.log.debug('chains: %s, chains_modeller: %s' % (chains, chains_modeller,) )
                sasa_score = analyze_structure.get_sasa()['A']
                
                if len(sasa_score) != \
                (sql.decode_domain(uniprot_template.domain_def)[1] - sql.decode_domain(uniprot_template.domain_def)[0] + 1):
                    model_errors += 'pops score length mismatch, '
            except error.PopsError as e:
                self.log.error('PopsError:')
                self.log.error(e.error)
                model_errors += 'pops error, '
                sasa_score = ''
            
            uniprot_model.sasa_score = ','.join(['%.2f' % (x*100) for x in sasa_score])            
            uniprot_model.model_errors = model_errors
        
        
        if type(uniprot_domain) == sql.UniprotDomainPair:
            uniprot_model = sql.UniprotDomainPairModel()
            uniprot_model.uniprot_domain_pair_id = uniprot_template.uniprot_domain_pair_id
        
            uniprot_model.model_filename = pdbFile_wt
            uniprot_model.chain_1 = chains_modeller[0]
            uniprot_model.chain_2 = chains_modeller[1]
            uniprot_model.het_flag_1 = HETflag[chains_modeller[0]]
            uniprot_model.het_flag_2 = HETflag[chains_modeller[1]]
            uniprot_model.switch_chain = switch_chain            
            
            uniprot_model.norm_dope = normDOPE_wt
            
            # Get interacting amino acids and interface area          
            analyze_structure = class_analyze_structure.AnalyzeStructure(modeller_path, pdbFile_wt, [chains_modeller[0]], self.log, domain_defs=None)
            interacting_aa = analyze_structure.get_interacting_aa()
            interacting_aa_1 = interacting_aa[chains_modeller[0]]
            uniprot_model.interacting_aa_1 = ','.join(['%i' % x for x in interacting_aa_1])
            
            interacting_aa_2 = interacting_aa[chains_modeller[1]]
            uniprot_model.interacting_aa_2 = ','.join(['%i' % x for x in interacting_aa_2])
            
            interface_area = analyze_structure.get_interface_area()
            uniprot_model.interface_area_hydrophobic = interface_area[0]
            uniprot_model.interface_area_hydrophilic = interface_area[1]
            uniprot_model.interface_area_total = interface_area[2]
            
            
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
                       save_path
                       ):
        outFile = self.tmpPath + self.unique + 'modeller/outFile_wildtype'
        
        for a, b in zip(alignments, chains):
            self.log.debug('Alignment: ')
            self.log.debug(a)
            self.log.debug('Chain: %s' % b)
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

        modeller_path = self.tmpPath + self.unique + 'modeller/'
        os.chdir(modeller_path) # from os

        modeller = mod.modeller([inFile], 
                                modeller_target_id, 
                                modeller_template_id, 
                                save_path, # path_to_pdb_for_modeller
                                self.tmpPath + self.unique + '/', # path to folders with executables
                                self.modeller_runs,
                                loopRefinement=True)
        normDOPE, pdbFile, knotted = modeller.run()

        os.chdir(self.PWD) # go back to the working directory
        
        return normDOPE, pdbFile, knotted