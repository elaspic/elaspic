import logging

from Bio import SeqIO

import elaspic
from kmtools import py_tools, system_tools

logger = logging.getLogger(__name__)


class Sequence:
    """Class for calculating sequence level features."""

    @py_tools.log_calls
    def __init__(self, sequence_file, results=None, mutations=None):
        self.sequence_file = sequence_file
        seqrecord = SeqIO.read(self.sequence_file, 'fasta')
        self.protein_id = system_tools.slugify(seqrecord.id)
        self.protein_sequence = str(seqrecord.seq)
        # Results
        self.results = {} if results is None else results.copy()
        self.mutations = {} if mutations is None else mutations.copy()
        # if provean_supset_file:
        #     self.provean_supset_file = provean_supset_file
        # else:
        #     self.provean_supset_file = op.join(
        #         elaspic.CONFIGS['provean_temp_dir'],
        #         system_tools.slugify(self.protein_id + '_provean_supset')
        #     )
        # self.provean_supset_length = None

    def mutate(self, mutation):
        if mutation in self.mutations:
            return self.mutations[mutation]

        if mutation[0] != self.protein_sequence[int(mutation[1:-1]) - 1]:
            logger.error('sequence: {}'.format(self.protein_sequence))
            logger.error('mutation: {}'.format(mutation))
            raise elaspic.exc.MutationMismatchError()

        results = dict(
            protein_id=self.protein_id,
            mutation=mutation,
            provean_score=self.run_provean(mutation),
            matrix_score=self.score_pairwise(mutation[0], mutation[-1])
        )
        return results
