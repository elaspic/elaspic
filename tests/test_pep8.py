import os.path as op
import pep8
import pytest

REPO_ROOT_DIR = op.dirname(op.dirname(op.abspath(__file__)))


def test_pep8_compliance():
    """Test that we conform to the PEP8 standard."""
    pep8style = pep8.StyleGuide(path=REPO_ROOT_DIR)
    list_of_files = [
        op.join(REPO_ROOT_DIR, f) for f in op.join(REPO_ROOT_DIR, 'elaspic') if f.endswith('.py')
    ]
    result = pep8style.check_files(list_of_files)
    assert result.total_errors < 60  # don't ask for much


if __name__ == '__main__':
    pytest.main([__file__, '-svx'])
