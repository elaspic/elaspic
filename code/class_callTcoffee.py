# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 19:03:01 2012

@author: niklas
"""
import subprocess
from os import environ
from Bio import AlignIO
from class_error import TcoffeeError


class tcoffee_alignment:
    """
    Alignes sequences using t_coffee in expresso mode
    
    input:  tmpPath     type: string        set a unique path for t_coffee
                                            (see __call_tcoffee)
            seqFile     type: file          containing two sequences in fasta format
            seqID       type: string        sequence ID of the target sequence
                                            (in pipeline for modelling with modeller)
            seqFile2    type: file          containing two sequences in fasta format
            seqID2      type: string        sequence ID of the target sequence
                                            (in pipeline for modelling with modeller)
    
    return: Biopython alignment object
    

    """
    def __init__(self, tmpPath, unique, alnPath, seqFiles, seqIDs):

        self.blastPath = tmpPath
        self.tmpPath = tmpPath + unique + '/'
        
        self.seqFiles = seqFiles
        self.seqIDs = seqIDs
        
        self.alnPath = alnPath
        self.alnFormat = 'clustal'
        
    
    def align(self):
        """
        start t_coffee in expresso mode and return the alignment
        """
        alignments = list()
        for seq in self.seqFiles:
            alignments.append(self.__call_tcoffee(seq))
            self.__writeAlignment(alignments[-1], self.seqIDs[ self.seqFiles.index(seq) ])
        
        return alignments



    def __writeAlignment(self, alignment, seqID):
        """
        write the alignment in clustal format to the folder specified for the class instance
        """
        alignmentID = [aln.id for aln in alignment]
#        AlignIO.write(alignment, self.alnPath + seqID, self.alnFormat)
        AlignIO.write(alignment, self.alnPath + alignmentID[0], self.alnFormat) # AS changed from above so that the alignments with the same template are not overwritten


    def __call_tcoffee_system_command(self, alignInFile, out, mode):
        # to be able to run parallel instances of T_Coffee the environment
        # variables have to be set to unique paths for each T_Coffee call
        # also make sure to change the directory to a unique one bevore
        # calling T_Coffee
        my_env = environ.copy()
        my_env['HOME_4_TCOFFEE'] = self.tmpPath + 'tcoffee/'
        my_env['TMP_4_TCOFFEE'] = self.tmpPath + 'tcoffee/tmp/'
        my_env['CACHE_4_TCOFFEE'] = self.tmpPath + 'tcoffee/cache/'
        my_env['LOCKDIR_4_TCOFFEE'] = self.tmpPath + 'tcoffee/lck/'
        my_env['ERRORFILE_4_TCOFFEE'] = self.tmpPath + 't_coffee.ErrorReport'
        my_env['BLASTDB'] = self.blastPath + 'blast/pdbaa_db/'
        if mode == 'expresso':
            system_command = 'cd ' + self.tmpPath + ' && t_coffee' + \
                             ' -mode expresso' + \
                             ' -seq ' + alignInFile + \
                             ' -blast_server=LOCAL' + \
                             ' -method clustalw_pair slow_pair' + \
                             ' -pdb_db=pdbaa -protein_db=uniprot' + \
                             ' -outorder=input' + \
                             ' -output fasta_aln' + \
                             ' -quiet -no_warning' + \
                             ' -outfile=' + out
        if mode == 't_coffee':
            system_command = 'cd ' + self.tmpPath + ' && t_coffee' + \
                             ' -mode expresso' + \
                             ' -method TMalign_pair' + \
                             ' -seq ' + alignInFile + \
                             ' -blast_server=LOCAL' + \
                             ' -method clustalw_pair slow_pair' + \
                             ' -pdb_db=pdbaa -protein_db=pdbaa' + \
                             ' -outorder=input' + \
                             ' -output fasta_aln' + \
                             ' -quiet -no_warning' + \
                             ' -outfile=' + out

        return system_command, my_env
        
    def __call_tcoffee(self, alignInFile, GAPOPEN=-0.0, GAPEXTEND=-0.0):
        """
        calls t_coffee in expresso mode (make sure BLAST is installed locally)
        
        input:  alignInFile     type: string        file containing two sequences
                                                    in fasta format
                tmpPATH         type: string        used to set a unique path
                                                    for the t_coffee tmp, cache,
                                                    and lock directory (for 
                                                    paralellization)
                GAPOPEN         type: int or str    see t_coffee manual
                GAPEXTEND       type: int or str    see t_coffee manual
                
                
        return: Biopython multiple sequence alignment object
        """

        out = self.tmpPath + 'sequenceAlignment.aln'
        
        # try the alignment in expresso mode (structure based with sap alignment)
        system_command, my_env = self.__call_tcoffee_system_command(alignInFile, out, 'expresso')
        childProcess = subprocess.Popen(system_command, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        shell=True, 
                                        env=my_env
                                        )
        result, error = childProcess.communicate()
        rc = childProcess.returncode

        # check if tcoffee had an unexpected exit and if not, create and return 
        # the alignment object
        if rc == 0:
            alignment = AlignIO.read(out, 'fasta')
            return alignment
        else:
            # if aligning two identical sequences it may happen that for each
            # sequence the same pdb template is selected. If that happens
            # sap fails to align and the alignment does not work
            for line in error.split('\n'):
                if 'SAP failed to align' in line:
                    # if it happens because the same PDB was taken by blast
                    # it means that the sequences are fairly identical and 
                    # and normal t_coffee mode should be accurate enough
                    # the line looks like this: pid 5479 -- SAP failed to align: 1KU6A.pdb against 1KU6A.pdb [T-COFFEE:WARNING]
                    
                    # try running it with tmalign method
                    system_command, my_env = self.__call_tcoffee_system_command(alignInFile, out, 't_coffee')
                    childProcess = subprocess.Popen(system_command, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    shell=True, 
                                    env=my_env
                                    )
                    result, error = childProcess.communicate()
                    rc = childProcess.returncode
                    if rc == 0:
                        alignment = AlignIO.read(out, 'fasta')
                        if len(alignment) != 2:
                            print 'Alignment len not 2', out
                        return alignment
            # raise an error if it still didn't work
            print 'out', out
            print 'alignInFile', alignInFile
            print 'error', error
            raise TcoffeeError(result, error, alignInFile)
