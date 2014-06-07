# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 19:03:01 2012

@author: niklas
"""

from os import environ
from Bio import AlignIO
import errors
import helper_functions as hf

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
    def __init__(
            self, global_temp_path, tmpPath, alnPath, seqFiles, seqIDs, n_cores,
            pdb_path, mode, log):

        self.global_temp_path = global_temp_path
        self.tmpPath = tmpPath
        self.seqFiles = seqFiles
        self.seqIDs = seqIDs
        self.alnPath = alnPath
        self.alnFormat = 'clustal'
        self.n_cores = n_cores
        self.pdb_path = pdb_path
        self.mode = mode
        self.log = log


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
#        AlignIO.write(alignment, self.alnPath + seqID, self.alnFormat)
        try:
            AlignIO.write(alignment, self.alnPath + alignment[0].id + '_' + alignment[1].id + '.aln', self.alnFormat) # AS changed from above so that the alignments with the same template are not overwritten
        except IndexError as e:
            raise errors.EmptyPDBSequenceError(str(type(e)) + ': ' + e.message)


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
        my_env['BLASTDB'] = self.global_temp_path + 'blast/db/'
        my_env['PDB_DIR'] = '/home/kimlab1/database_data/pdb/'
#        my_env['NO_REMOTE_PDB_DIR'] = '1'
#        print (
#        'export tmp_path=`pwd` && '
#        'export HOME_4_TCOFFEE=$tmp_path/tcoffee/ && '
#        'export TMP_4_TCOFFEE=$tmp_path/tcoffee/ && '
#        'export CACHE_4_TCOFFEE=$tmp_path/tcoffee/ && '
#        'export LOCKDIR_4_TCOFFEE=$tmp_path/tcoffee/ && '
#        'export ERRORFILE_4_TCOFFEE=$tmp_path/tcoffee/ && '
#        'export BLASTDB=$tmp_path/../../blast/db/ && '
#        'export PDB_DIR=/home/kimlab1/database_data/pdb/data/data/structures/divided/pdb/ && '
#        'export NO_REMOTE_PDB_DIR=1')

        multi_core_option = '{}'.format(self.n_cores) if self.n_cores and self.n_cores > 1 else 'no'
        n_core_option = '{}'.format(self.n_cores) if self.n_cores else '1'
        if mode == 'expresso':
            system_command = (
                't_coffee' +
                ' -mode expresso' +
                ' -method sap_pair,mustang_pair' +
                ' -seq ' + alignInFile +
                ' -blast_server=LOCAL' +
                ' -pdb_db=pdbaa -protein_db=nr' +
                ' -outorder=input' +
                ' -output fasta_aln' +
                ' -quiet -no_warning' +
                ' -outfile=' + out +
#                ' -multi_core no' +
                ' -multi_core ' + multi_core_option + #AS changed !!!
                ' -n_core ' + n_core_option) #AS changed !!!
#
#        't_coffee -other_pg extract_from_pdb 32c2A.pdb > template.pdb '
#        't_coffee -mode 3dcoffee -method sap_pair,mustang_pair,TMalign_pair -blast_server=LOCAL -pdb_db=pdbaa -protein_db=nr -outorder=input -output fasta_aln -outfile tcoffee_output.aln -seq seqfiles.fasta -pdb_min_sim=20 -template_file seqfiles.template '
        if mode == 't_coffee':
            system_command = (
                't_coffee' +
                ' -mode expresso' +
                ' -method clustalw_pair,slow_pair' +
                ' -seq ' + alignInFile +
                ' -blast_server=LOCAL' +
                ' -pdb_db=pdbaa -protein_db=pdbaa' +
                ' -outorder=input' +
                ' -output fasta_aln' +
                ' -quiet -no_warning' +
                ' -outfile=' + out +
#                ' -multi_core no')
                ' -multi_core ' + multi_core_option + #AS changed !!!
                ' -n_core ' + n_core_option) #AS changed !!!

        if mode == 'quick':
            system_command = (
                't_coffee' +
                ' -mode quickaln' +
                ' -method clustalw_pair,slow_pair' +
                ' -seq ' + alignInFile +
                ' -blast_server=LOCAL' +
                ' -pdb_db=pdbaa -protein_db=nr' +
                ' -outorder=input' +
                ' -output fasta_aln' +
                ' -quiet -no_warning' +
                ' -outfile=' + out +
#                ' -multi_core no' +
                ' -multi_core ' + multi_core_option + #AS changed !!!
                ' -n_core ' + n_core_option) #AS changed !!!
        return system_command, my_env


    def __call_tcoffee(self, alignInFile, GAPOPEN=-0.0, GAPEXTEND=-0.0, recursion_counter=0):
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
        # mode should be 'expresso'
        out = self.tmpPath + 'sequenceAlignment.aln'

        # try the alignment in expresso mode (structure based with sap alignment)
        system_command, my_env = self.__call_tcoffee_system_command(alignInFile, out, self.mode)
        child_process = hf.run_subprocess_locally(self.tmpPath, system_command, env=my_env)
        result, error_message = child_process.communicate()
        return_code = child_process.returncode

        # check if tcoffee had an unexpected exit and if not, create and return
        # the alignment object
        if return_code == 0:
            alignment = AlignIO.read(out, 'fasta')
            return alignment
        else:
            self.log.error('Structural alignment failed with the following error:')
            self.log.error(error_message)
            self.log.error('Running quickaln alignment instead...')
            system_command, my_env = self.__call_tcoffee_system_command(alignInFile, out, 'quick')
            child_process = hf.run_subprocess_locally(self.tmpPath, system_command, env=my_env)
            result, error_message = child_process.communicate()
            return_code = child_process.returncode
            if return_code == 0:
                alignment = AlignIO.read(out, 'fasta')
                if len(alignment) != 2:
                    self.log.error('Alignment length not 2 for file %s' % out)
                return alignment
            self.log.error('Even quickaln didn\'t work. Cannot create an alignment. Giving up.')
            raise errors.TcoffeeError(result, error_message, alignInFile, system_command)
