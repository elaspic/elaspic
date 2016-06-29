import os.path as op
from os import environ
import time
import logging
import shutil
from Bio import SeqIO
from . import conf, errors, helper, structure_tools

logger = logging.getLogger(__name__)


class TCoffee(object):
    """Alignes sequences using t_coffee in expresso mode."""

    def __init__(self, alignment_fasta_file, mode, pdb_file=None):
        """
        """
        self.alignment_id = op.splitext(op.basename(alignment_fasta_file))[0]
        self.alignment_fasta_file = op.join(
            conf.CONFIGS['tcoffee_dir'], self.alignment_id + '.fasta')
        if alignment_fasta_file != self.alignment_fasta_file:
            shutil.copy(alignment_fasta_file, self.alignment_fasta_file)

        self.target_seqrecord, self.template_seqrecord = (
            list(SeqIO.parse(self.alignment_fasta_file, 'fasta'))
        )
        self.mode = mode

        if pdb_file is not None:
            self.pdb_id, self.pdb_file = self._clean_pdb(pdb_file)

            # Write a template file for the sequences to be aligned
            self.alignment_template_file = (
                op.join(conf.CONFIGS['tcoffee_dir'], '{}.template_list'.format(self.alignment_id))
            )
            with open(self.alignment_template_file, 'w') as fh:
                fh.writelines([
                    ">" + self.target_seqrecord.id + " _P_ " + self.pdb_id.upper() + "\n"
                    ">" + self.template_seqrecord.id + " _P_ " + self.pdb_id.upper() + "\n"
                ])

    def _clean_pdb(self, pdb_file):
        """Write a template PDB file in a format that is compatible with t_coffee."""
        message = (
            "Cleaning pdb {} to serve as a template for t_coffee...".format(pdb_file)
        )
        logger.debug(message)

        pdb_id = structure_tools.get_pdb_id(pdb_file)
        pdb_file_new = op.join(
            conf.CONFIGS['tcoffee_dir'], pdb_id + '.pdb')

        system_command = (
            "t_coffee -other_pg extract_from_pdb {} > {}".format(pdb_file, pdb_file_new)
        )
        p = helper.run(system_command, cwd=conf.CONFIGS['tcoffee_dir'])
        if p.returncode != 0:
            logger.error("Error cleaning pdb!")
            logger.error("System command: '{}'".format(system_command))
            logger.error("Result:\n{}".format(p.stdout))
            logger.error("Error message:\n{}".format(p.stderr))
        time.sleep(0.2)
        return pdb_id, pdb_file_new

    def _get_tcoffee_system_command(
            self, alignment_fasta_file, alignment_template_file, alignment_output_file, mode):
        """.

        Parameters
        ----------
        alignment_fasta_file : str
            Name of the file that contains the sequences to be aligned in fasta format.
        alignment_template_file : str
            Name of the file that contains the structureal templates to be used for
            structure-assisted alignments.
        alignment_output_file : str
            Name of the file where the alignment should be saved.
        mode : str
            T-coffee mode the should be run.

        Returns
        --------
        system_command : str
            System call that runs tcoffee.
        tcoffee_env : str
            A dictionary of environment variables to run tcoffee in a thread-safe manner.
        """
        # Environment variables
        # To be able to run parallel instances of t_coffee, the environment
        # variables have to be set to unique paths for each t_coffee call.
        # Also, make sure to change the directory to a unique one bevore
        # calling t_coffee.
        tcoffee_env = environ.copy()
        tcoffee_env['HOME_4_TCOFFEE'] = op.join(conf.CONFIGS['tcoffee_dir'])
        tcoffee_env['TMP_4_TCOFFEE'] = op.join(conf.CONFIGS['tcoffee_dir'], 'tmp')
        tcoffee_env['CACHE_4_TCOFFEE'] = op.join(conf.CONFIGS['tcoffee_dir'], 'cache')
        tcoffee_env['LOCKDIR_4_TCOFFEE'] = op.join(conf.CONFIGS['tcoffee_dir'], 'lck')
        tcoffee_env['ERRORFILE_4_TCOFFEE'] = (
            op.join(conf.CONFIGS['tcoffee_dir'], 't_coffee.ErrorReport')
        )
        tcoffee_env['BLASTDB'] = conf.CONFIGS['blast_db_dir']
        tcoffee_env['PDB_DIR'] = conf.CONFIGS['pdb_dir']
        tcoffee_env['NO_REMOTE_PDB_DIR'] = '1'

        # Print a command that can be used to set environmental variables
        t_coffee_environment_variables = [
            'HOME_4_TCOFFEE', 'TMP_4_TCOFFEE', 'CACHE_4_TCOFFEE', 'LOCKDIR_4_TCOFFEE',
            'ERRORFILE_4_TCOFFEE', 'BLASTDB', 'PDB_DIR', 'NO_REMOTE_PDB_DIR']
        exports = [
            'export {}={}'.format(x, tcoffee_env.get(x, '$' + x))
            for x in t_coffee_environment_variables
        ]
        message = (
            "\nSystem command for setting environmental variables:\n" + ' && '.join(exports)
        )
        logger.debug(message)

        # ### System command
        # Use the following command to clean the pdb file (add headers, etc.)
        # 't_coffee -other_pg extract_from_pdb 32c2A.pdb > template.pdb '

        # Use the folllowing command to perform the most accurate alignment
        # 't_coffee -mode 3dcoffee -method sap_pair,mustang_pair,TMalign_pair
        # -blast_server=LOCAL -pdb_db=pdbaa -protein_db=nr -outorder=input
        # -output fasta_aln -outfile tcoffee_output.aln -seq seqfiles.fasta
        # -pdb_min_sim=20 -template_file seqfiles.template '

        multi_core_option = (
            '{}'.format(conf.CONFIGS['n_cores'])
            if conf.CONFIGS['n_cores'] and int(conf.CONFIGS['n_cores']) > 1 else 'no'
        )
        n_core_option = '{}'.format(conf.CONFIGS['n_cores']) if conf.CONFIGS['n_cores'] else '1'
        protein_db = op.join(conf.CONFIGS['blast_db_dir'], 'nr')
        pdb_db = op.join(conf.CONFIGS['blast_db_dir'], 'pdbaa')
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
                # " -quiet" +
                # " -no_warning" +
                " -outfile=" + alignment_output_file +
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
                ' -outfile=' + alignment_output_file +
                ' -cache ignore ' +
                ' -multi_core ' + multi_core_option +  # AS changed !!!
                ' -n_core ' + n_core_option  # AS changed !!!
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
                ' -outfile=' + alignment_output_file +
                ' -multi_core ' + multi_core_option +  # AS changed !!!
                ' -n_core ' + n_core_option  # AS changed !!!
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
                ' -outfile=' + alignment_output_file +
                ' -multi_core ' + multi_core_option +  # AS changed !!!
                ' -n_core ' + n_core_option  # AS changed !!!
            )
        return system_command, tcoffee_env

    def align(self, GAPOPEN=-0.0, GAPEXTEND=-0.0):
        """Call t_coffee (make sure BLAST is installed locally!).

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

        Returns
        --------
        alignment_output_file : str
            Name of file which contains the alignment in fasta format.
        """
        # try the alignment in expresso mode (structure based with sap alignment)
        alignment_output_file = op.join(conf.CONFIGS['tcoffee_dir'], self.alignment_id + '.aln')
        system_command, tcoffee_env = self._get_tcoffee_system_command(
            self.alignment_fasta_file, self.alignment_template_file, alignment_output_file,
            self.mode)

        # Perform t_coffee alignment
        logger.debug("\nTCoffee system command:\n{}".format(system_command))
        p = helper.run(system_command, cwd=conf.CONFIGS['tcoffee_dir'], env=tcoffee_env)
        logger.debug("t_coffee results:\n{}".format(p.stdout))
        error_message_summary_idx = p.stderr.find(
            '*                        MESSAGES RECAPITULATION')
        if error_message_summary_idx == -1:
            error_message_summary = ''
        else:
            error_message_summary = p.stderr[error_message_summary_idx:]
        logger.debug("t_coffee errors:\n{}".format(error_message_summary.strip()))

        # Check if tcoffee had an unexpected exit and if not, create and return
        # the alignment object
        if p.returncode == 0:
            logger.info("Successfully made the alignment")
            return alignment_output_file
        else:
            logger.error(
                'Structural alignment failed with the following error: {}'.format(p.stderr))
            logger.error('Running quickalign alignment instead...')
            system_command, tcoffee_env = self._get_tcoffee_system_command(
                self.alignment_fasta_file, self.alignment_template_file, alignment_output_file,
                'quick')
            p = helper.run(system_command, cwd=conf.CONFIGS['tcoffee_dir'], env=tcoffee_env)
            if p.returncode == 0:
                return alignment_output_file
            else:
                logger.error('Even quickaln didn\'t work. Cannot create an alignment. Giving up.')
                raise errors.TcoffeeError(
                    p.stdout, p.stderr, self.alignment_fasta_file, system_command)
