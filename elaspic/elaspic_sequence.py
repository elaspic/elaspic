import re
import os
import os.path as op
import requests
import psutil
import time
import shutil
import logging
import atexit
import subprocess
import shlex
import six

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo

from . import conf, helper, errors

logger = logging.getLogger(__name__)

CANONICAL_AMINO_ACIDS = 'ARNDCEQGHILKMFPSTWYV'


# %% Sequence tools
def download_uniport_sequence(uniprot_id, output_dir):
    """
    """
    output_file = op.join(output_dir, uniprot_id + '.fasta')

    # If the file already exists, do nothing...
    if op.isfile(output_file):
        logger.debug('Sequence file {} already exists...'.format(output_file))
        return output_file

    logger.debug('Downloading sequence {}...'.format(uniprot_id + '.fasta'))
    address = 'http://www.uniprot.org/uniprot/{}.fasta'.format(uniprot_id)
    r = requests.get(address)
    if r.status_code == 200:
        with open(output_file, 'w') as ofh:
            ofh.write(r.text)
        return output_file


def convert_basestring_to_seqrecord(sequence, sequence_id='id'):
    if any([isinstance(sequence, string_type) for string_type in six.string_types]):
        seqrec = SeqRecord(Seq(sequence), id=str(sequence_id))
    elif isinstance(sequence, Seq):
        seqrec = SeqRecord(sequence, id=str(sequence_id))
    elif isinstance(sequence, SeqRecord):
        seqrec = sequence
    else:
        raise Exception("Wrong class type %s for ``sequence``" % str(type(sequence)))
    return seqrec


# %%
class Sequence:
    """Class for calculating sequence level features."""

    def __init__(self, sequence_file, provean_supset_file=None):
        """.

        Parameters
        ----------
        sequence_file : str
            Full filename of the file containing the protein sequence.
            Only fasta format is supported.
        provean_supset_file : str
            Full path and
        """
        logger.debug('Initialising a Sequence instance with parameters:')
        logger.debug('sequence_file: {}'.format(sequence_file))
        logger.debug('provean_supset_file: {}'.format(provean_supset_file))

        self.sequence_file = sequence_file
        self.seqrecord = SeqIO.read(self.sequence_file, 'fasta')
        self.protein_id = helper.slugify(self.seqrecord.id)
        self.sequence = str(self.seqrecord.seq)

        # Provean supset
        if provean_supset_file is not None and provean_supset_file != self.provean_supset_file:
            shutil.copy(provean_supset_file, self.provean_supset_file)
            shutil.copy(provean_supset_file, self.provean_supset_file + '.fasta')
        if not self.provean_supset_exists:
            logger.debug('Calculating provean supset...')
            self._build_provean_supset()
        else:
            logger.debug('Provean supset is already calculated!')
        self.provean_supset_length = self._get_provean_supset_length()

        # Mutations
        self.mutations = {}

    def mutate(self, mutation):
        if mutation in self.mutations:
            return self.mutations[mutation]

        if mutation[0] != self.sequence[int(mutation[1:-1]) - 1]:
            logger.error('sequence: {}'.format(self.sequence))
            logger.error('mutation: {}'.format(mutation))
            raise errors.MutationMismatchError()

        results = dict(
            protein_id=self.protein_id,
            mutation=mutation,
            provean_score=self.run_provean(mutation),
            matrix_score=self.score_pairwise(mutation[0], mutation[-1])
        )
        return results

    @property
    def provean_supset_file(self):
        return op.join(
            conf.CONFIGS['sequence_dir'],
            helper.slugify(self.protein_id + '_provean_supset')
        )

    @property
    def provean_supset_exists(self):
        return (
            op.isfile(self.provean_supset_file) and
            op.isfile(self.provean_supset_file + '.fasta')
        )

    @property
    def result(self):
        result = dict(
            protein_id=self.protein_id,
            sequence=self.sequence,
            sequence_file=op.relpath(self.sequence_file, conf.CONFIGS['unique_temp_dir']),
            provean_supset_exists=self.provean_supset_exists,
            provean_supset_file=op.relpath(
                self.provean_supset_file, conf.CONFIGS['unique_temp_dir']),
            provean_supset_length=self.provean_supset_length,
            mutations=self.mutations,
        )
        return result

    def _build_provean_supset(self):
        """
        """
        logger.debug('Building Provean supporting set. This might take a while...')
        atexit.register(_clear_provean_temp)

        # Get the required parameters
        any_position = 0
        while self.sequence[any_position] not in CANONICAL_AMINO_ACIDS:
            any_position += 1
        first_aa = self.sequence[any_position]
        mutation = '{0}{1}{0}'.format(first_aa, any_position + 1)

        # Run provean
        provean_score = self._run_provean(
            mutation, save_supporting_set=True, check_mem_usage=True
        )
        return provean_score
#        provean_supset_length = None
#        for line in result.split('\n'):
#            if 'Number of supporting sequences used:' in line:
#                provean_supset_length = int(line.split()[-1])
#        if provean_supset_length is None:
#            logger.error('Provean supporting set length could not be estimated.
#            This is an error!')
#            logger.error('Provean result: {}'.format(result))
#            logger.error('Provean error_message: {}'.format(error_message))
#            logger.error('Provean return_code: {}'.format(return_code))
#            logger.error('Protein sequence: {}'.format(self.sequence))
#            logger.error('Protein mutation: {}'.format(mutation))
#
#        logger.info('Provean supset length: {}'.format(provean_supset_length))
#        return provean_supset_length

    def _get_provean_supset_length(self):
        provean_supset_length = 0
        with open(self.provean_supset_file) as fh:
            for line in fh:
                if line and not line.startswith('#'):
                    provean_supset_length += 1
        return provean_supset_length

    def run_provean(self, mutation, *args, **kwargs):
        n_tries = 0
        while n_tries < 5:
            n_tries += 1
            try:
                provean_score = self._run_provean(mutation, *args, **kwargs)
                break
            except errors.ProveanError as e:
                bad_ids = re.findall("Entry not found in BLAST database: '(.*)'", e.args[0])
                provean_supset_data = []
                with open(self.provean_supset_file, 'rt') as ifh:
                    for line in ifh:
                        if any([(gi_id in line) for gi_id in bad_ids]):
                            logger.debug(
                                "Removing line '{}' from the provean supset file..."
                                .format(line.strip()))
                        else:
                            provean_supset_data.append(line)
                with open(self.provean_supset_file, 'wt') as ofh:
                    ofh.writelines(provean_supset_data)
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
        errors.ProveanError
            Can raise this exception only if ``check_mem_usage`` is set to ``True``.
        """
        if check_mem_usage:
            # Get initial measurements of how much virtual memory and disk space is availible
            disk_space_availible = psutil.disk_usage(
                conf.CONFIGS['provean_temp_dir']).free / (1024**3)
            logger.debug('Disk space availible: {:.2f} GB'.format(disk_space_availible))
            if disk_space_availible < 5:
                raise errors.ProveanError(
                    'Not enough disk space ({:.2f} GB) to run provean'
                    .format(disk_space_availible))
            memory_availible = psutil.virtual_memory().available / float(1024)**3
            logger.debug('Memory availible: {:.2f} GB'.format(memory_availible))
            if memory_availible < 0.5:
                raise errors.ProveanError(
                    'Not enough memory ({:.2f} GB) to run provean'
                    .format(memory_availible))

        # Create a file with mutation
        mutation_file = op.join(conf.CONFIGS['sequence_dir'], '{}.var'.format(mutation))
        with open(mutation_file, 'w') as ofh:
            ofh.write(mutation)

        # Run provean
        system_command = (
            "provean " +
            " -q '{}' ".format(self.sequence_file) +
            " -v '{}' ".format(mutation_file) +
            " -d " + op.join(conf.CONFIGS['blast_db_dir'], 'nr') +
            " --tmp_dir '{}' ".format(conf.CONFIGS['provean_temp_dir']) +
            " --num_threads {} ".format(conf.CONFIGS['n_cores']) +
            " --psiblast '{}' ".format(helper.get_which('psiblast')) +
            " --blastdbcmd '{}' ".format(helper.get_which('blastdbcmd')) +
            " --cdhit '{}' ".format(helper.get_which('cd-hit'))
        )

        if self.provean_supset_exists:
            # use supporting set
            system_command += " --supporting_set '{}' ".format(self.provean_supset_file)
        else:
            system_command += " --save_supporting_set '{}' ".format(self.provean_supset_file)

        logger.debug(system_command)
        p = subprocess.Popen(
            shlex.split(system_command),
            cwd=conf.CONFIGS['sequence_dir'],
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
            disk_space_availible_now = (
                psutil.disk_usage(conf.CONFIGS['provean_temp_dir']).free / float(1024)**3
            )
            if disk_space_availible_now < 5:  # less than 5 GB of free disk space left
                raise errors.ProveanResourceError(
                    'Ran out of disk space and provean had to be terminated ({} GB used)'
                    .format(disk_space_availible - disk_space_availible_now),
                    child_process_group_id)
            memory_availible_now = (
                psutil.virtual_memory().available / float(1024)**3
            )
            if memory_availible_now < 0.5:
                raise errors.ProveanResourceError(
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
            raise errors.ProveanError(stderr)

        return provean_score

    # === Other sequence scores ===

    def score_pairwise(self, seq1, seq2, matrix=None, gap_s=None, gap_e=None):
        """Get the BLOSUM (or what ever matrix is given) score."""
        matrix = matrix or getattr(MatrixInfo, conf.CONFIGS['matrix_type'])
        gap_s = gap_s or conf.CONFIGS['gap_start']
        gap_e = gap_e or conf.CONFIGS['gap_extend']

        score = 0
        gap = False
        for i in range(len(seq1)):
            pair = (seq1[i], seq2[i])
            if not gap:
                if '-' in pair:
                    gap = True
                    score += gap_s
                else:
                    score += self._score_match(pair, matrix)
            else:
                if '-' not in pair:
                    gap = False
                    score += self._score_match(pair, matrix)
                else:
                    score += gap_e
        return score

    def _score_match(self, pair_match, matrix_match):
        """
        """
        if pair_match not in matrix_match:
            return matrix_match[(tuple(reversed(pair_match)))]
        else:
            return matrix_match[pair_match]


def _clear_provean_temp():
    provean_temp_dir = conf.CONFIGS['provean_temp_dir']
    logger.info("Clearning provean temporary files from '{}'...".format(provean_temp_dir))
    for filename in os.listdir(provean_temp_dir):
        file_path = os.path.join(provean_temp_dir, filename)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(e)
