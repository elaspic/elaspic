# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 19:03:01 2012

@author: kimlab
"""
from __future__ import absolute_import
from __future__ import unicode_literals
from builtins import object

from os import environ
import time
import six

from Bio import SeqIO
from Bio import AlignIO

from . import errors
from . import helper_functions as hf


class tcoffee_alignment(object):
    """
    Alignes sequences using t_coffee in expresso mode.
    """
    def __init__(
            self, global_temp_path, unique_temp_folder, uniprot_seqrecord,
            pdb_seqrecord, n_cores, pdb_path, mode, logger):

        self.global_temp_path = global_temp_path
        self.unique_temp_folder = unique_temp_folder
        self.uniprot_seqrecord = uniprot_seqrecord
        self.pdb_seqrecord = pdb_seqrecord
        self.alnFormat = 'clustal'
        self.n_cores = n_cores
        self.pdb_path = pdb_path
        self.mode = mode
        self.logger = logger


    def align(self):
        """ Start t_coffee in expresso mode and return the alignment
        """
        alignments = [self._call_tcoffee()]
        return alignments



    def __call_tcoffee_system_command(self, alignment_fasta_file, alignment_template_file, out, mode):
        # To be able to run parallel instances of t_coffee, the environment
        # variables have to be set to unique paths for each t_coffee call.
        # Also, make sure to change the directory to a unique one bevore
        # calling t_coffee.
        my_env = environ.copy()
        my_env['HOME_4_TCOFFEE'] = self.unique_temp_folder + 'tcoffee/'
        my_env['TMP_4_TCOFFEE'] = self.unique_temp_folder + 'tcoffee/tmp/'
        my_env['CACHE_4_TCOFFEE'] = self.unique_temp_folder + 'tcoffee/cache/'
        my_env['LOCKDIR_4_TCOFFEE'] = self.unique_temp_folder + 'tcoffee/lck/'
        my_env['ERRORFILE_4_TCOFFEE'] = self.unique_temp_folder + 't_coffee.ErrorReport'
        my_env['BLASTDB'] = self.global_temp_path + 'blast/db/'
        # my_env['PDB_DIR'] = '/home/kimlab1/database_data/pdb/'
        # my_env['PDB_DIR'] = self.pdb_path
        my_env['PDB_DIR'] = self.unique_temp_folder
        my_env['NO_REMOTE_PDB_DIR'] = '1'
        t_coffee_environment_variables = [
            'HOME_4_TCOFFEE', 'TMP_4_TCOFFEE', 'CACHE_4_TCOFFEE', 'LOCKDIR_4_TCOFFEE',
            'ERRORFILE_4_TCOFFEE', 'BLASTDB', 'PDB_DIR', 'NO_REMOTE_PDB_DIR']
        self.logger.debug("System command for setting environmental variables:\n" +
            ' && '.join(['export {}={}'.format(x, my_env.get(x, '$'+x)) for x in t_coffee_environment_variables]))

        # Use the following command to clean the pdb file (add headers, etc.)
        # 't_coffee -other_pg extract_from_pdb 32c2A.pdb > template.pdb '

        # Use the folllowing command to perform the most accurate alignment
        # 't_coffee -mode 3dcoffee -method sap_pair,mustang_pair,TMalign_pair
        # -blast_server=LOCAL -pdb_db=pdbaa -protein_db=nr -outorder=input
        # -output fasta_aln -outfile tcoffee_output.aln -seq seqfiles.fasta
        # -pdb_min_sim=20 -template_file seqfiles.template '

        multi_core_option = '{}'.format(self.n_cores) if self.n_cores and self.n_cores > 1 else 'no'
        n_core_option = '{}'.format(self.n_cores) if self.n_cores else '1'
        protein_db = self.global_temp_path + 'blast/db/nr'
        pdb_db = self.global_temp_path + 'blast/db/pdbaa'
        if mode == '3dcoffee':
            system_command = (
                "t_coffee " +
                " -seq " + alignment_fasta_file +
                " -method=sap_pair,mustang_pair,TMalign_pair " +
                " -blast_server=LOCAL " +
                " -protein_db=" + protein_db +
                " -pdb_db=" + pdb_db +
                " -outorder=input" +
                " -output=fasta_aln" +
                " -pdb_min_sim=30" +
                #" -quiet" +
                #" -no_warning" +
                " -outfile=" + out +
                " -multi_core=no" +
                " -n_core=" + n_core_option +
                " -template_file=" + alignment_template_file
            )
        if mode == 'expresso':
            system_command = (
                't_coffee' +
                ' -mode expresso' +
                ' -method sap_pair' +
                ' -seq ' + alignment_fasta_file +
                ' -blast_server=LOCAL' +
                " -protein_db=" + protein_db +
                " -pdb_db=" + pdb_db +
                ' -outorder=input' +
                ' -output fasta_aln' +
                ' -quiet -no_warning' +
                ' -outfile=' + out +
                ' -cache ignore ' +
                ' -multi_core ' + multi_core_option + #AS changed !!!
                ' -n_core ' + n_core_option #AS changed !!!
            )
        if mode == 't_coffee':
            system_command = (
                't_coffee' +
                ' -mode expresso' +
                ' -method clustalw_pair,slow_pair' +
                ' -seq ' + alignment_fasta_file +
                ' -blast_server=LOCAL' +
                " -protein_db=" + protein_db +
                " -pdb_db=" + pdb_db +
                ' -outorder=input' +
                ' -output fasta_aln' +
                ' -quiet -no_warning' +
                ' -outfile=' + out +
                ' -multi_core ' + multi_core_option + #AS changed !!!
                ' -n_core ' + n_core_option #AS changed !!!
            )
        if mode == 'quick':
            system_command = (
                't_coffee' +
                ' -mode quickaln' +
                ' -method clustalw_pair,slow_pair' +
                ' -seq ' + alignment_fasta_file +
                ' -blast_server=LOCAL' +
                " -protein_db=" + protein_db +
                " -pdb_db=" + pdb_db +
                ' -outorder=input' +
                ' -output fasta_aln' +
                ' -quiet -no_warning' +
                ' -outfile=' + out +
                ' -multi_core ' + multi_core_option + #AS changed !!!
                ' -n_core ' + n_core_option #AS changed !!!
            )
        return system_command, my_env






    def _call_tcoffee(self, GAPOPEN=-0.0, GAPEXTEND=-0.0):
        """ Calls t_coffee (make sure BLAST is installed locally)
        Parameters
        ----------
        alignment_fasta_file : string
            A file containing the fasta sequences to be aligned
        alignment_template_file : string
            A file containing the structural templates for the fasta sequences
            described above
        GAPOPEN : int or str
            See t_coffee manual
        GAPEXTEND : int or str
            See t_coffee manual

        return: Biopython multiple sequence alignment object
        """
        ### Spliced from another function

        # sequence_ids = [uniprot_seqrecord.id, pdb_seqrecord.id]
        alignment_fasta_file = self.unique_temp_folder + 'seqfiles.fasta'
        alignment_template_file = self.unique_temp_folder + 'seqfiles.template_list'

        # Write a fasta file with sequences to be aligned
        with open(alignment_fasta_file, 'w') as fh:
            SeqIO.write([self.uniprot_seqrecord, self.pdb_seqrecord], fh, 'fasta')

        # Write a template file for the sequences to be aligned
        with open(alignment_template_file, 'w') as fh:
            fh.writelines([
                ">" + self.uniprot_seqrecord.id + " _P_ " + self.pdb_seqrecord.id.upper() + "\n"
                ">" + self.pdb_seqrecord.id + " _P_ " + self.pdb_seqrecord.id.upper() + "\n"
            ])


        # try the alignment in expresso mode (structure based with sap alignment)
        out = self.unique_temp_folder + 'sequenceAlignment.aln'
        system_command, my_env = self.__call_tcoffee_system_command(
            alignment_fasta_file, alignment_template_file, out, self.mode)


        # Write a template PDB file in a format that is compatible with t_coffee
        self.logger.debug("Cleaning pdb {} to serve as a template for t_coffee...".format(self.pdb_seqrecord.id + '.pdb'))
        format_pdb_system_command = "t_coffee -other_pg extract_from_pdb {} > {}".format(
            self.pdb_seqrecord.id + '.pdb', self.pdb_seqrecord.id.upper() + '.pdb')
        child_process = hf.run_subprocess_locally(self.unique_temp_folder, format_pdb_system_command, env=my_env)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        if child_process.returncode:
            self.logger.error(
                "Error cleaning pdb!\nSystem command: '{}'\nResult: '{}'\nError message: '{}'"
                .format(system_command, result, error_message))
        time.sleep(0.2)


        # Perform t_coffee alignment
        self.logger.debug("t_coffee system command:\n{}".format(system_command))
        child_process = hf.run_subprocess_locally(self.unique_temp_folder, system_command, env=my_env)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        self.logger.debug("t_coffee results:\n{}".format(result.strip()))
        error_message_summary_idx = error_message.find('*                        MESSAGES RECAPITULATION')
        if error_message_summary_idx == -1:
            error_message_summary = ''
        else:
            error_message_summary = error_message[error_message_summary_idx:]
        self.logger.debug("t_coffee errors:\n{}".format(error_message_summary.strip()))
        return_code = child_process.returncode

        # Check if tcoffee had an unexpected exit and if not, create and return
        # the alignment object
        if return_code == 0:
            self.logger.info("Successfully made the alignment")
            alignment = AlignIO.read(out, 'fasta')
            return alignment
        else:
            self.logger.error('Structural alignment failed with the following error: {}'.format(error_message))
            self.logger.error('Running quickalign alignment instead...')
            system_command, my_env = self.__call_tcoffee_system_command(
                alignment_fasta_file, alignment_template_file, out, 'quick')
            child_process = hf.run_subprocess_locally(self.unique_temp_folder, system_command, env=my_env)
            result, error_message = child_process.communicate()
            if six.PY3:
                result = str(result, encoding='utf-8')
                error_message = str(error_message, encoding='utf-8')
            return_code = child_process.returncode
            if return_code == 0:
                alignment = AlignIO.read(out, 'fasta')
                if len(alignment) != 2:
                    self.logger.error('Alignment length not 2 for file %s' % out)
                return alignment
            self.logger.error('Even quickaln didn\'t work. Cannot create an alignment. Giving up.')
            raise errors.TcoffeeError(result, error_message, alignment_fasta_file, system_command)
