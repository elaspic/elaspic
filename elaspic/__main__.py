import argparse
import logging
import logging.config

import coloredlogs

from elaspic.cli.elaspic_database import configure_database_parser
from elaspic.cli.elaspic_run import configure_run_parser
from elaspic.cli.elaspic_train import configure_train_parser

logger = logging.getLogger(__name__)

LOGGING_LEVELS = {
    None: logging.ERROR,
    0: logging.ERROR,
    1: logging.WARNING,    # -v
    2: logging.INFO,       # -vv
    3: logging.DEBUG,      # -vvv
}


def main():
    parser = argparse.ArgumentParser(
        prog='elaspic',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '-v', '--verbose', action='count',
        help=('Specify verbosity level'))
    sub_parsers = parser.add_subparsers(
        title='command',
        help=''
    )
    configure_run_parser(sub_parsers)
    configure_database_parser(sub_parsers)
    configure_train_parser(sub_parsers)
    args = parser.parse_args()
    # Configure logging
    default_format = '%(levelname)-7s %(name)s[%(lineno)d] %(message)s'
    logger.debug("Before setting level")
    if True:
        coloredlogs.install(
            fmt=default_format,
            level=LOGGING_LEVELS[args.verbose])
    else:
        logging.basicConfig(
            format=default_format,
            level=LOGGING_LEVELS[args.verbose])
    logger.warning("After setting level")
    #
    if 'func' not in args.__dict__:
        args = parser.parse_args(['--help'])
    args.func(args)


if __name__ == '__main__':
    import sys
    sys.exit(main())
