import itertools
import logging

from Bio.SubsMat import MatrixInfo

from ._abc import Mutator

logger = logging.getLogger(__name__)


class AlignmentAnalyzer(Mutator):

    _results_keys = []

    def __init__(self, sequence1, sequence2):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        super().__init__()

    def build(self):
        pass

    def mutate(self, mutation):
        if mutation in self.mutations:
            return self.mutations[mutation]
        else:
            self.mutations[mutation] = {
                'score': self.run_provean(mutation)
            }
            return self.mutations[mutation]

    # Helper
    def score_alignment(self, matrix='blosum80', gap_s=-16, gap_e=-4):
        """Get the BLOSUM (or what ever matrix is given) score."""
        matrix = getattr(MatrixInfo, matrix)

        score = 0
        gap = False
        for aa1, aa2 in itertools.zip_longest(self.sequence1, self.sequence2, fillvalue='-'):
            if not gap:
                if aa1 == '-' or aa2 == '-':
                    gap = True
                    score += gap_s
                else:
                    score += get_matrix_score(aa1, aa2, matrix)
            else:
                if aa1 == '-' or aa2 == '-':
                    gap = False
                    score += get_matrix_score(aa1, aa2, matrix)
                else:
                    score += gap_e
        return score


def get_matrix_score(self, seq1, seq2, matrix):
    if (seq1, seq2) in matrix:
        return matrix[(seq1, seq2)]
    else:
        return matrix[(seq2, seq1)]
