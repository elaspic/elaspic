import os.path as op
import glob
import logging
try:
    # flake8 < 3
    from flake8.engine import get_style_guide
except ImportError:
    # flake8 >= 3
    from flake8.api.legacy import get_style_guide

logger = logging.getLogger(__name__)
PKG_ROOT_DIR = op.dirname(op.dirname(op.abspath(__file__)))


def test_pep8_compliance():
    """Test that we conform to the PEP8 standard."""
    # Hack to get flake8 to shut up
    logging.getLogger('flake8').setLevel(logging.WARNING)
    logging.getLogger('flake8.api').setLevel(logging.WARNING)
    logging.getLogger('flake8.main').setLevel(logging.WARNING)
    logging.getLogger('flake8.checker').setLevel(logging.WARNING)
    logging.getLogger('flake8.plugins').setLevel(logging.WARNING)
    logging.getLogger('flake8.options').setLevel(logging.WARNING)

    pep8style = get_style_guide(
        config_file=op.join(PKG_ROOT_DIR, 'setup.cfg'),
        path=PKG_ROOT_DIR)
    list_of_files = glob.glob(op.join(PKG_ROOT_DIR, '**', '*.py'), recursive=True)
    logger.debug("Testing {} files for PEP8 compience.".format(len(list_of_files)))
    result = pep8style.check_files(list_of_files)
    assert result.total_errors == 0
