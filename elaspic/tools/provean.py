import functools
import os
import os.path as op
import re
import shlex
import subprocess
import time

import psutil

import Bio.SeqIO
import elaspic
from kmtools import py_tools, system_tools

from ._abc import SequenceAnalyzer, ToolError

logger = py_tools.get_logger(__name__)


CANONICAL_AMINO_ACIDS = 'ARNDCEQGHILKMFPSTWYV'


class ProveanError(ToolError):
    pass


class ProveanResourceError(ToolError):
    def __init__(self, message, child_process_group_id):
        Exception.__init__(self, message)
        self.child_process_group_id = child_process_group_id


class Provean(SequenceAnalyzer):

    _result_slots = ['supporting_set_file', 'supporting_set_length']

    def build(self):
        logger.info("{} {}", self.sequence, self.sequence.id)
        # Write sequence file
        self.sequence_file = op.join(
            elaspic.CONFIGS['sequence_dir'],
            "{}.fasta".format(self.sequence.id))
        with open(self.sequence_file, 'wt') as ofh:
            Bio.SeqIO.write(self.sequence, ofh, 'fasta')
        # Build the supporting set
        if self.done:
            logger.debug('Provean supset is already calculated!')
        else:
            self._build_provean_supset()
            assert self.done

    @functools.lru_cache(maxsize=512)
    def analyze(self, mutation):
        return {'provean_score': self.run_provean(mutation)}

    # Helper
    @property
    def supporting_set_file(self):
        return op.join(self.tempdir, self.sequence.id + '_provean_supset')

    @property
    def supporting_set_length(self):
        provean_supset_length = 0
        with open(self.supporting_set_file) as fh:
            for line in fh:
                if line and not line.startswith('#'):
                    provean_supset_length += 1
        assert provean_supset_length > 0
        return provean_supset_length

    def _build_provean_supset(self, mutation=None):
        logger.debug('Building a Provean supporting set. This might take a while...')

        # Get the required parameters
        any_position = 0
        while self.sequence[any_position] not in CANONICAL_AMINO_ACIDS:
            any_position += 1
        first_aa = self.sequence[any_position]
        if mutation is None:
            mutation = '{0}{1}{0}'.format(first_aa, any_position)

        # Run provean to generate supporting set files
        self._run_provean(mutation, save_supporting_set=True, check_mem_usage=True)

    def run_provean(self, mutation, *args, **kwargs):
        n_tries = 0
        provean_score = None
        while n_tries < 5:
            n_tries += 1
            try:
                provean_score = self._run_provean(mutation, *args, **kwargs)
                break
            except ProveanError as e:
                bad_ids = re.findall("Entry not found in BLAST database: '(.*)'", e.args[0])
                provean_supset_data = []
                with open(self.result['supporting_set_file'], 'rt') as ifh:
                    for line in ifh:
                        if any([(gi_id in line) for gi_id in bad_ids]):
                            logger.debug(
                                "Removing line '{}' from the provean supset file..."
                                .format(line.strip()))
                        else:
                            provean_supset_data.append(line)
                with open(self.result['supporting_set_file'], 'wt') as ofh:
                    ofh.writelines(provean_supset_data)
        if provean_score is None:
            # Recalculate provean supporting set
            provean_score = self._build_provean_supset(mutation)
        return provean_score

    def _run_provean(self, mutation, save_supporting_set=False, check_mem_usage=False):
        """.

        Provean results look something like this::

            #[23:28:34] clustering subject sequences...
            #[23:28:34] selecting clusters...
            #[23:28:34] 0 subject sequences in 0 clusters were selected for supporting sequences.
            #[23:28:34] use the query itself as a supporting sequence
            #[23:28:34] loading subject sequences from a FASTA file...
            #[23:28:34] scores were computed based on the query sequence itself.
            ## Number of clusters:	1
            ## Number of supporting sequences used:	1
            #[23:28:34] computing delta alignment scores...
            #[23:28:34] printing PROVEAN scores...
            ### PROVEAN scores ##
            ## VARIATION	SCORE
            #M1A	-6.000

        Parameters
        ----------
        domain_mutation : string
            Mutation in domain coordinates (i.e. relative to the start of the domain)

        Returns
        -------
        list
            [result, error_message, return_code] --
            The output from running a provean system command.

        Raises
        ------
        ProveanError
            Can raise this exception only if ``check_mem_usage`` is set to ``True``.
        """
        if check_mem_usage:
            # Get initial measurements of how much virtual memory and disk space is availible
            disk_space_availible = psutil.disk_usage(self.tempdir).free / (1024**3)
            logger.debug('Disk space availible: {:.2f} GB'.format(disk_space_availible))
            if disk_space_availible < 5:
                raise ProveanError(
                    'Not enough disk space ({:.2f} GB) to run provean'
                    .format(disk_space_availible))
            memory_availible = psutil.virtual_memory().available / float(1024)**3
            logger.debug('Memory availible: {:.2f} GB'.format(memory_availible))
            if memory_availible < 0.5:
                raise ProveanError(
                    'Not enough memory ({:.2f} GB) to run provean'
                    .format(memory_availible))

        # Create a file with mutation
        mutation_file = op.join(elaspic.CONFIGS['sequence_dir'], '{}.var'.format(mutation))
        with open(mutation_file, 'w') as ofh:
            # Need to add one because Provean accepts 1-based mutations
            ofh.write(mutation[0] + str(int(mutation[1:-1]) + 1) + mutation[-1])

        # Run provean
        system_command = (
            "provean " +
            " -q '{}' ".format(self.sequence_file) +
            " -v '{}' ".format(mutation_file) +
            " -d " + op.join(elaspic.CONFIGS['blast_db_dir'], 'nr') +
            " --tmp_dir '{}' ".format(self.tempdir) +
            " --num_threads {} ".format(elaspic.CONFIGS['n_cores']) +
            " --psiblast '{}' ".format(system_tools.which('psiblast')) +
            " --blastdbcmd '{}' ".format(system_tools.which('blastdbcmd')) +
            " --cdhit '{}' ".format(system_tools.which('cd-hit'))
        )

        if self.done:
            # use supporting set
            system_command += " --supporting_set '{}' ".format(self.result['supporting_set_file'])
        else:
            system_command += " --save_supporting_set '{}' ".format(self.supporting_set_file)

        logger.debug(system_command)
        p = subprocess.Popen(
            shlex.split(system_command),
            cwd=elaspic.CONFIGS['sequence_dir'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )

        logger.debug('Parent group id: {}'.format(os.getpgrp()))
        child_process_group_id = os.getpgid(p.pid)
        logger.debug('Child group id: {}'.format(child_process_group_id))

        # Keep an eye on provean to make sure it doesn't do anything crazy
        time.sleep(5)
        while check_mem_usage and p.poll() is None:
            disk_space_availible_now = psutil.disk_usage(self.tempdir).free / float(1024)**3
            if disk_space_availible_now < 5:  # less than 5 GB of free disk space left
                raise ProveanResourceError(
                    'Ran out of disk space and provean had to be terminated ({} GB used)'
                    .format(disk_space_availible - disk_space_availible_now),
                    child_process_group_id)
            memory_availible_now = (
                psutil.virtual_memory().available / float(1024)**3
            )
            if memory_availible_now < 0.5:
                raise ProveanResourceError(
                    'Ran out of RAM and provean had to be terminated ({} GB left)'
                    .format(memory_availible - memory_availible_now),
                    child_process_group_id)
            time.sleep(60)  # Wait for 1 minute before checking again

        # Collect the results and check for errors
        stdout, stderr = p.communicate()
        stdout = stdout.strip()
        stderr = stderr.strip()
        logger.debug(stdout)

        # Extract provean score from the results message
        provean_score = None
        result_list = stdout.split('\n')
        for i in range(len(result_list)):
            if re.findall('# VARIATION\s*SCORE', result_list[i]):
                provean_score = float(result_list[i + 1].split()[-1])
                break

        if p.returncode != 0 or provean_score is None:
            logger.error('return_code: {}'.format(p.returncode))
            logger.error('provean_score: {}'.format(provean_score))
            logger.error('error_message: {}'.format(stderr))
            raise ProveanError(stderr)

        if not self.done:
            self.result['supporting_set_file'] = self.supporting_set_file
            self.result['supporting_set_length'] = self.supporting_set_length
            assert self.done

        return provean_score
