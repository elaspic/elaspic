import elaspic.elaspic_model
import pytest


@pytest.mark.parametrize("alignment, scores", [
    [('AAAAA', 'AAAAA'), (1.0, 1.0, None, None)],
])
def test_analyze_alignment(alignment, scores):
    assert elaspic.elaspic_model.analyze_alignment(alignment) == scores
