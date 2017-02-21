import functools
import itertools
import logging

from Bio.SubsMat import MatrixInfo

from ._abc import SequenceAnalyzer

import elaspic

logger = logging.getLogger(__name__)


class Conservation(SequenceAnalyzer):

    _result_slots = []

    def build(self):
        self.matrix_type = elaspic.CONFIGS.get('matrix_type', 'blosum80')
        self.matrix = getattr(MatrixInfo, self.matrix_type)
        self.gap_s = elaspic.CONFIGS.get('gap_start', -16)
        self.gap_e = elaspic.CONFIGS.get('gap_extend', -4)

    @functools.lru_cache(maxsize=512)
    def analyze(self, mutation):
        score = self.get_matrix_score(mutation[0], mutation[-1])
        return {'{}_score'.format(self.matrix_type): score}

    def analyze_alignment(self, alignment):
        """Get the BLOSUM (or what ever matrix is given) score."""
        assert len(alignment) == 2
        score = 0
        gap = False
        for aa1, aa2 in itertools.zip_longest(alignment[0], alignment[1], fillvalue='-'):
            if not gap:
                if aa1 == '-' or aa2 == '-':
                    gap = True
                    score += self.gap_s
                else:
                    score += self.get_matrix_score(aa1, aa2)
            else:
                if aa1 == '-' or aa2 == '-':
                    gap = False
                    score += self.get_matrix_score(aa1, aa2)
                else:
                    score += self.gap_e
        return {'{}_score'.format(self.matrix_type): score}

    def get_matrix_score(self, seq1, seq2):
        if (seq1, seq2) in self.matrix:
            return self.matrix[(seq1, seq2)]
        else:
            return self.matrix[(seq2, seq1)]
