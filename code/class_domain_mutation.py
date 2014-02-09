# -*- coding: utf-8 -*-

import os
import subprocess
import tarfile

from Bio import SeqIO

import class_error as error
import class_sql as sql
import class_analyze_structure 

from class_pdbTemplate import pdbTemplate
from class_pysicoChemicalProperties import pysiChem
from class_foldX import foldX


class GetMutation(object):
    """
    """
    
    def __init__(self, tmpPath, unique, pdbPath, db, log,
                  foldX_WATER, buildModel_runs, PWD, matrix, gap_s, gap_e):
        """
        input:
        tmpPath             type 'str'
        unique              type 'str'
        pdbPath             type 'str'
        savePDB             type 'str'
        """
        self.tmpPath = tmpPath
        self.unique = unique + '/'
        self.pdbPath = pdbPath
        self.db = db
        
        # get the logger from the parent and add a handler
        self.log = log
        
        self.foldX_WATER = foldX_WATER
        self.buildModel_runs = buildModel_runs
        self.PWD = PWD
        
        self.matrix = matrix
        self.gap_s = gap_s
        self.gap_e = gap_e


    def __call__(self, uniprot_domain, uniprot_template, uniprot_model, mutated_uniprot_id, mutation):
        """
        AS: renamed the __run() function with some restructuring.
        
        input:  save_path        path to PDB files for modeller
                self.tmpPath    path to store some tmp data (not tcoffee)
                self.unique     set environment for tcoffee (for parallelization)
                pdbFile_wt      name of the wildtype structure file (either raw or from modeller)
                mutations_pdb   list of mutations in `Chain`_`AAwt``AAnumber``AAmut` format, where AAnumber is 
                                the index of the AA in the pdb
                chains_modeller     list of chains in pdbFile_wt, with the first chain in the list being mutated
                mutations_modeller  same as above, only compensated for chain switching and other changes made by modeller
                sequences       list of sequences of each chain in the pdb
                switch_chain    whether or not the order of the sequences was reversed during modelling
        """
        
        if type(uniprot_domain) == sql.UniprotDomain:
            self.log.error("single domain")
            domain_def = uniprot_template.domain_def
            
#            pdb_id = uniprot_template.domain.pdb_id
            chains_pdb = [uniprot_template.domain.pdb_chain, ]
            chains_modeller = [uniprot_model.chain, ]
                
            uniprot_sequences = [self.db.get_uniprot_sequence(mutated_uniprot_id), ]
            domain_sequences = [uniprot_sequences[0][sql.decode_domain(domain_def)[0]-1:
                                                        sql.decode_domain(domain_def)[1]], ]
            alignment, __ = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)
            alignment_id = uniprot_template.alignment_id
            HETflag = {chains_pdb[0]: uniprot_model.het_flag}
            switch_chain = uniprot_model.switch_chain
            
        
        elif type(uniprot_domain) == sql.UniprotDomainPair:
            
            if mutated_uniprot_id == uniprot_domain.uniprot_domain_1.uniprot_id:
                self.log.debug('mutated_uniprot_id: %s, uniprot_id: %s' 
                    % (mutated_uniprot_id, uniprot_domain.uniprot_domain_1.uniprot_id,) )
                self.log.debug("domain pair right order")
                
                domain_def = uniprot_template.domain_def_1
                
#                pdb_id = uniprot_template.domain_1.pdb_id
                chains_pdb = [uniprot_template.domain_1.pdb_chain, uniprot_template.domain_2.pdb_chain]
                chains_modeller = [uniprot_model.chain_1, uniprot_model.chain_2]
                
                uniprot_sequences = [self.db.get_uniprot_sequence(mutated_uniprot_id), 
                                     self.db.get_uniprot_sequence(uniprot_domain.uniprot_domain_2.uniprot_id)]
                domain_sequences = [uniprot_sequences[0][sql.decode_domain(uniprot_template.domain_def_1)[0]-1:
                                                            sql.decode_domain(uniprot_template.domain_def_1)[1]], 
                                    uniprot_sequences[1][sql.decode_domain(uniprot_template.domain_def_2)[0]-1:
                                                            sql.decode_domain(uniprot_template.domain_def_2)[1]]]
                alignment_1, alignment_2 = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)
                alignment = alignment_1
                alignment_id = uniprot_template.alignment_id_1
                # HETflag is a dictionary, doesn't matter the order. But key-value should match
                HETflag = {uniprot_template.domain_1.pdb_chain: uniprot_model.het_flag_1,
                           uniprot_template.domain_2.pdb_chain: uniprot_model.het_flag_2}
                switch_chain = uniprot_model.switch_chain
            
            
            elif mutated_uniprot_id == uniprot_domain.uniprot_domain_2.uniprot_id:
                self.log.debug('mutated_uniprot_id: %s, uniprot_id: %s' 
                    % (mutated_uniprot_id, uniprot_domain.uniprot_domain_1.uniprot_id,) )
                self.log.debug("domain pair reversed order")
                
                alignment_2, alignment_1 = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)
                
                domain_def = uniprot_template.domain_def_2
                
#                pdb_id = uniprot_template.domain_2.pdb_id
                chains_pdb = [uniprot_template.domain_1.pdb_chain, uniprot_template.domain_2.pdb_chain]
                chains_modeller = [uniprot_model.chain_1, uniprot_model.chain_2]
                
                uniprot_sequences = [self.db.get_uniprot_sequence(uniprot_domain.uniprot_domain_2.uniprot_id),
                                     self.db.get_uniprot_sequence(uniprot_domain.uniprot_domain_1.uniprot_id)]
                domain_sequences = [uniprot_sequences[0][sql.decode_domain(uniprot_template.domain_def_2)[0]-1:
                                                            sql.decode_domain(uniprot_template.domain_def_2)[1]],
                                    uniprot_sequences[1][sql.decode_domain(uniprot_template.domain_def_1)[0]-1:
                                                            sql.decode_domain(uniprot_template.domain_def_1)[1]]]
                alignment_1, alignment_2 = self.db.get_alignment(uniprot_template, uniprot_domain.path_to_data)
                alignment = alignment_2
                alignment_id = uniprot_template.alignment_id_2
                # HETflag is a dictionary, doesn't matter the order. But key-value should match
                HETflag = {uniprot_template.domain_1.pdb_chain: uniprot_model.het_flag_1,
                           uniprot_template.domain_2.pdb_chain: uniprot_model.het_flag_2}
                switch_chain = not uniprot_model.switch_chain
        
        mutation_errors = ''
        
        save_path = self.tmpPath + uniprot_domain.path_to_data
        pdbFile_wt = uniprot_model.model_filename
        
        chains = chains_modeller
        sequences = domain_sequences
        #######################################################################
        
        # Convert mutation from uniprot domain numbering to pdb domain numbering
        mutation_uniprot_domain_no_chain = int(mutation[1:-1]) - sql.decode_domain(domain_def)[0] + 1 # +1 to have the numbering staring with 1                    
        
        mutation_pdb_no_chain = self.map_to_pdb_sequence(alignment, alignment_id, mutation_uniprot_domain_no_chain)
        if mutation_pdb_no_chain == 'in gap':
            raise error.MutationOutsideDomain()
        
        mutation_uniprot = chains_pdb[0] + '_' + mutation[0] + str(mutation_uniprot_domain_no_chain) + mutation[-1]
        mutations_uniprot = [mutation_uniprot,]
        
        mutation_pdb = chains_pdb[0] + '_' + mutation[0] + str(mutation_pdb_no_chain) + mutation[-1]
        mutations_pdb = [mutation_pdb,]
        
        # Rename chains in mutation by modeller output
        mutations_modeller = self.setChainID(chains_modeller, mutations_pdb, HETflag, switch_chain, True) 
        

        
        # Have to figure out a better way of doing this anyway
#        if not contacts_chain1[chains_pdb[1]]:
#            raise error.NotInteracting
            
        #######################################################################
        # Logging
        self.log.debug("save_path:")
        self.log.debug(save_path)          
        self.log.debug("pdbFile_wt:")
        self.log.debug(pdbFile_wt)
        self.log.debug("chains:")
        self.log.debug(chains)
        
        
        #######################################################################
        ## 1st: get a snippet of AAs around the mutation, as well as the mutation residue (no chain)
        # This is required input for FoldX
        if not switch_chain:
            mutations_foldX = self.prepareMutationFoldX(sequences[0], mutations_uniprot)
        else:
            mutations_foldX = self.prepareMutationFoldX(sequences[1], mutations_uniprot)

        self.log.debug("mutations_foldX:")
        self.log.debug(mutations_foldX)   

       
        #######################################################################
        ## 2nd: use the 'Repair' feature of FoldX to optimise the structure
        foldX_path = self.tmpPath + self.unique + 'FoldX/'
        try:
            os.chdir(foldX_path) # from os
            fX = foldX(self.tmpPath + self.unique, 
                               save_path + pdbFile_wt, 
                               chains[0], 
                               self.unique, 
                               self.buildModel_runs,
                               self.foldX_WATER
                               )
            repairedPDB_wt = fX.run('RepairPDB')
        except:
            raise
        finally:
            os.chdir(self.PWD)
        
        
        #######################################################################
        ## 3rd: introduce the mutation using FoldX
        
        # compile a list of mutations
        mutCodes = list()
        for mut in mutations_foldX:
            mutCodes.append(mut[1])
        
        fX_wt = foldX(self.tmpPath + self.unique, 
                      repairedPDB_wt, 
                      chains[0], 
                      self.unique, 
                      self.buildModel_runs,
                      self.foldX_WATER
                      )
        
        # do the mutation with foldX
        repairedPDB_wt_list, repairedPDB_mut_list = fX_wt.run('BuildModel', mutCodes)
        
        # Copy the foldX wildtype pdb file (use the first model if there are multiple)
        model_filename_wt = mutated_uniprot_id + '_' + mutation + '/' +  repairedPDB_wt_list[0].split('/')[-1]
        model_filename_mut = mutated_uniprot_id + '_' + mutation + '/MUT_' +  repairedPDB_mut_list[0].split('/')[-1]
        subprocess.check_call('mkdir -p ' + save_path + mutated_uniprot_id + '_' + mutation + '/', shell=True)
        subprocess.check_call('cp ' + repairedPDB_wt_list[0] + ' ' + save_path + model_filename_wt, shell=True)
        subprocess.check_call('cp ' + repairedPDB_mut_list[0] + ' ' + save_path + model_filename_mut, shell=True)
        
        
        #######################################################################
        ## 4th: set up the classes for the wildtype and the mutant structures
        fX_wt_list = list()       
        for wPDB in repairedPDB_wt_list:
            fX_wt_list.append(foldX(self.tmpPath + self.unique, 
                                    wPDB, 
                                    chains[0], 
                                    self.unique, 
                                    self.buildModel_runs,
                                    self.foldX_WATER
                                    ))        
        
        fX_mut_list = list()
        for mPDB in repairedPDB_mut_list:
            fX_mut_list.append(foldX(self.tmpPath + self.unique, 
                                     mPDB, 
                                     chains[0], 
                                     self.unique, 
                                     self.buildModel_runs,
                                     self.foldX_WATER
                                     ))        
        
        
        #######################################################################
        ## 5th: calculate the energy for the wildtype
        FoldX_StabilityEnergy_wt = list()
        FoldX_AnalyseComplex_wt  = list()
        for fX_wt in fX_wt_list:
            # stability of the interaction
            if len([ item for item in chains if item != '' ]) == 1:
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
            if len([ item for item in chains if item != '' ]) == 1:
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
        # AnalyseComplex
        if len(FoldX_AnalyseComplex_wt) == 1:
            AnalyseComplex_energy_wt = FoldX_AnalyseComplex_wt[0]
            AnalyseComplex_energy_mut = FoldX_AnalyseComplex_mut[0]
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
            AnalyseComplex_energy_mut = [ [] for i in range(len(FoldX_AnalyseComplex_mut[0])) ]
            for item in FoldX_AnalyseComplex_mut:
                i = 0
                for element in item:
                    AnalyseComplex_energy_mut[i].append(element)
                    i += 1
            AnalyseComplex_energy_mut = ':'.join([','.join(item) for item in AnalyseComplex_energy_mut])

        
        # Stability
        if len(FoldX_StabilityEnergy_wt) == 1:
            Stability_energy_wt = FoldX_StabilityEnergy_wt[0]
            Stability_energy_mut = FoldX_StabilityEnergy_mut[0]
        else:
            # zip the output for the wildtype
            Stability_energy_wt = [ [] for i in range(len(FoldX_StabilityEnergy_wt[0])) ]
            for item in FoldX_StabilityEnergy_wt:
                i = 0
                for element in item:
                    Stability_energy_wt[i].append(element)
                    i += 1
            Stability_energy_wt = ':'.join([','.join(item) for item in Stability_energy_wt])
            
            # zip the output for the mutant
            Stability_energy_mut = [ [] for i in range(len(FoldX_StabilityEnergy_mut[0])) ]
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
#                pops = analyze_structure.GetSASA('', self.tmpPath + self.unique + 'pops/')
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
        get_atomicContactVector = pysiChem(5.0, 4.0, self.tmpPath + self.unique)
        # mutations is of the form I_V70A
        # pysiChem need is like I70 (chain that is mutated and the position)
        # self.mutationss is a list which can contain more than one mutations
        mut_physChem = [ mut[0] + mut[3:-1] for mut in mutations_modeller ]
        physChem_mut = list()
        physChem_wt = list()
        physChem_mut_ownChain = list()
        physChem_wt_ownChain = list()
        os.chdir(self.PWD) # from os
        
        for item in repairedPDB_wt_list:
            # calculate the contact vector
            res_wt = [0, 0, 0, 0]
            res_wt_ownChain = [0, 0, 0, 0]
            i = -1
            for mut in mut_physChem:
                i += 1
                mut_snippet     = mutations_foldX[i][1].split('\n')[0]
                mut_snippet_pos = mutations_foldX[i][0]
                
                chains_complex = [ chain for chain in chains[0] ]
                atomicContactVector, atomicContactVector_ownChain, mutpos_tmp = get_atomicContactVector(item, mut[0], chains_complex, mut_snippet_pos, mut_snippet)
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
            i = -1
            for mut in mut_physChem:
                i += 1
                mut_snippet     = mutations_foldX[i][1].split('\n')[1]
                mut_snippet_pos = mutations_foldX[i][0]
                
                chains_complex = [ chain for chain in chains[0] ]
                atomicContactVector, atomicContactVector_ownChain, mutpos_tmp = get_atomicContactVector(item, mut[0], chains_complex, mut_snippet_pos, mut_snippet)
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
        # if self.templateFinding:
        mutpos_DSSP = mutpos_tmp
#        else:
#            mutpos_DSSP = mutations_modeller[0][3:-1]
        
        dssp_tmp_path = '/'.join(repairedPDB_wt_list[0].split('/')[:-1]) + '/'
        dssp_chain = [mutations_modeller[0][0],]
        dssp_position = mutpos_DSSP
        
        dssp_pdb_file_wt = repairedPDB_wt_list[0].split('/')[-1]
        dssp_residue_wt = mutations_modeller[0][2]
        analyze_structure_wt = class_analyze_structure.AnalyzeStructure(dssp_tmp_path, dssp_pdb_file_wt, dssp_chain, domain_defs=None, log=self.log)
        
        dssp_pdb_file_mut = repairedPDB_mut_list[0].split('/')[-1]
        dssp_residue_mut = mutations_modeller[0][-1]
        analyze_structure_mut = class_analyze_structure.AnalyzeStructure(dssp_tmp_path, dssp_pdb_file_mut, dssp_chain, domain_defs=None, log=self.log)

#        try:
        solvent_accessibility_wt, secondary_structure_wt = analyze_structure_wt.get_dssp(dssp_residue_wt, dssp_position)
        solvent_accessibility_mut, secondary_structure_mut = analyze_structure_mut.get_dssp(dssp_residue_mut, dssp_position)
#        except:
#            print "DSSP DID NOT WORK!"
#            self.log.error('DSSP did not work!')
#            solvent_accessibility_wt, secondary_structure_wt = '-1', '-1'
#            solvent_accessibility_mut, secondary_structure_mut = '-1', '-1'
#            mutation_errors += 'dssp did not work; '
            
        
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
        for mutation in mutations_modeller:
            fromAA = mutation.split('_')[1][0]
            toAA   = mutation.split('_')[1][-1]
            matrix_score += self.score_pairwise(fromAA, toAA, self.matrix, self.gap_s, self.gap_e)
        
        
        #######################################################################
        #
        if type(uniprot_domain) == sql.UniprotDomain:
            uniprot_mutation = sql.UniprotDomainMutation()
            uniprot_mutation.uniprot_domain_id = uniprot_template.uniprot_domain_id
            
        elif type(uniprot_domain) == sql.UniprotDomainPair:
            uniprot_mutation = sql.UniprotDomainPairMutation()
            uniprot_mutation.uniprot_domain_pair_id = uniprot_template.uniprot_domain_pair_id            
        
        uniprot_mutation.uniprot_id = mutated_uniprot_id
        uniprot_mutation.mutation = mutation
        uniprot_mutation.mutation_errors = mutation_errors
        
        uniprot_mutation.model_filename_wt = model_filename_wt
        uniprot_mutation.model_filename_mut = model_filename_mut        
        
        uniprot_mutation.chain_modeller = mutations_modeller[0][0]
        uniprot_mutation.mutation_modeller = mutations_modeller[0][2:]
        
        uniprot_mutation.AnalyseComplex_energy_wt = AnalyseComplex_energy_wt
        uniprot_mutation.Stability_energy_wt = Stability_energy_wt            
        uniprot_mutation.AnalyseComplex_energy_mut = AnalyseComplex_energy_mut
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
        
        
        #######################################################################
        # Save alignments and modeller models to output database for storage
        # Move template files to the output folder as a tar archives
#        self.make_tarfile(self.HOME + self.outputPath + save_path + '_' + mutation + '.tar.bz2', 
#                          save_path[:-1])
        
        #######################################################################
        self.log.info('Finished processing template:')
        self.log.info(save_path.split('/')[-2])
                
        return uniprot_mutation


    
    def make_tarfile(self, output_filename, source_dir):
        with tarfile.open(output_filename, "w:bz2") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))


    def map_to_pdb_sequence(self, alignment, alignmentID, position):
        """
        Given an alignment and a position of the uniprot sequence, this function
        maps the position of the uniprot sequence to the pdb sequence
        !! Note that the position numbering start with 1 !!
        
        input
        alignment       type class 'Bio.Align.MultipleSeqAlignment'
        alignmentID     type 'str'
        position        type 'int'
        
        
        pdb_position    type 'int'

        """
        if alignment[0].id == alignmentID:
            alignment_protein = alignment[0]
            alignment_uniprot = alignment[1]
        elif alignment[1].id == alignmentID:
            alignment_protein = alignment[1]
            alignment_uniprot = alignment[0]
        else:
            self.log.error('Could not assign the alignment to pdb and uniprot correctly!')
            return 1

        # now get the position
        pdb_position = 0
        uniprot_position = 0
        for index in range(len(alignment_uniprot)):
            if uniprot_position >= int(position)-1 and not alignment_uniprot[index] == '-':
                check = index
                break
            elif alignment_uniprot[index] == '-':
                pdb_position += 1
                continue
            elif alignment_protein[index] == '-':
                uniprot_position += 1
                continue
            else:
                pdb_position += 1
                uniprot_position += 1
        
        # check if the uniprot position falls into a gap
        if alignment_protein[check] == '-':
            return 'in gap'
        else:
            return pdb_position + 1

 

    def score_pairwise(self, seq1, seq2, matrix, gap_s, gap_e):
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


       
    def prepareMutationFoldX(self, sequence, mutations):
        """
        create the sequence snippets for the mutation with FoldX
        also, record the position of the mutation in the snippet
        """
        mutations_foldX = list()
        for mutation in mutations:
            if int(mutation[3:-1]) < 7:
                left = 0
                m_pos = int(mutation[3:-1]) # the position of the muation in the snippet
            else:
                left = int(mutation[3:-1])-7
                m_pos = 7
            if int(mutation[3:-1]) >= len(sequence) - 8:
                right = len(sequence)
            else:
                right = int(mutation[3:-1])+8

            wt_seq = sequence.seq[left:right]
            mut_seq = sequence.seq[left:int(mutation[3:-1])-1] + mutation[-1] + sequence.seq[int(mutation[3:-1]):right]

            self.log.debug('mutation:')
            self.log.debug(mutation)
            self.log.debug('sequence:')
            self.log.debug(sequence)
            self.log.debug(sequence.seq[int(mutation[3:-1])-1-5:int(mutation[3:-1])+5])
            assert( sequence.seq[int(mutation[3:-1])-1] == mutation[2] )
            
            mutations_foldX.append((m_pos, str(wt_seq) + '\n' + str(mut_seq)))

        return mutations_foldX


    
    def setChainID(self, chains, mutations, HETflag, SWITCH_CHAIN, do_modelling):
        """
        modeller renames the chains, in order to keep track of the chain IDs
        they are renamed here to match the modeller output. Modeller labels the
        chains subsequently, i.e. A, B, C etc., thus, depending on wether HETATMs
        are present as seperate chain, the chain labels are A, B or A, C
        """
        mut = list()
        for mutation in mutations:
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
                            self.log.error('could not assign the chain for mutation correctly!')
                            mutChain = None
                    
            mut.append( mutChain + '_' + mutation[2:] )
        return mut



###############################################################################
# Methods f

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
        save_path =  self.tmpPath + self.unique + '/' + saveFolder + '/'
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
        template['sequences'] = (sequence,)
        template['chains_modeller'] = chains_get_pdb
        template['switch_chain'] = False
                
        return template, template_mutations


    
    def _get_pdb_sequence(self, pdbCode, chain):
        """
        Return the pdb file sequence (not SEQRES)
        """
        pdbPath = self.pdbPath
        domains = [['Null', 'Null'], ] # extract the full sequence

        pdb = pdbTemplate(pdbPath, pdbCode, chain, domains, self.tmpPath + self.unique)
        
        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()
        chainNumberingDomain = pdb.getChainNumberingNOHETATMS(chain)
        if chainNumberingDomain == []:
            raise error.pdbError('Could not get the pdb numbering for ' + pdbCode + '_' + chain)
        
        return next(SeqIO.parse(self.tmpPath + self.unique + '/' + pdbCode + chain + '.seq.txt', 'fasta')), chainNumberingDomain



    def _getCrystalStructure(self, pdbCode, chains, savePDB, FULL_PATH=False):
        domains = [['Null', 'Null'] for i in range(len(chains)) ]
        pdb = pdbTemplate(self.pdbPath, pdbCode, chains, domains, savePDB)
        HETATMsInChain_PDBnumbering, HETflag, chains_pdb_order = pdb.extract()
        
        SWITCH_CHAIN = False
        if chains != chains_pdb_order:
            SWITCH_CHAIN = True
            chains = chains_pdb_order
        
        return '-', pdbCode.split('.')[0] + ''.join(chains) + '.pdb', SWITCH_CHAIN, chains        