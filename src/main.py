'''
A Bioinformatics pipeline based on Ruffus.

Author: Khalid Mahmood (khalid.mahmood@unimelb.edu.au).

Copyright: 2015

See ../README.md for details about how to use the program.

Repository: https://github.com/khalidm/hiplexpipe
'''

from __future__ import print_function
from ruffus import *
import ruffus.cmdline as cmdline
from version import version
from name import program_name
import sys
from config import Config
from state import State
from logger import Logger
from pipeline import make_pipeline
import error_codes

# default place to save cluster job scripts
# (mostly useful for post-mortem debugging)
DEFAULT_JOBSCRIPT_DIR = 'jobscripts'
# default name of the pipeline configuration file
DEFAULT_CONFIG_FILE = 'pipeline.config'


def parse_command_line():
    '''Parse the command line arguments of the pipeline'''
    parser = cmdline.get_argparse(description='A variant discovery pipeline',
        ignored_args = ["version"] )
    parser.add_argument('--config', type=str, default=DEFAULT_CONFIG_FILE,
        help='Pipeline configuration file in YAML format, defaults to {}' \
            .format(DEFAULT_CONFIG_FILE))
    parser.add_argument('--jobscripts', type=str,
        default=DEFAULT_JOBSCRIPT_DIR,
        help='Directory to store cluster job scripts created by the ' \
             'pipeline, defaults to {}'.format(DEFAULT_JOBSCRIPT_DIR))
    parser.add_argument('--version', action='version',
        version='%(prog)s ' + version)
    return parser.parse_args()

def main():
    '''Initialise the pipeline, then run it'''
    # Parse command line arguments
    options = parse_command_line()
    # Initialise the logger
    logger = Logger(__name__, options.log_file, options.verbose)
    # Log the command line used to run the pipeline
    logger.info(' '.join(sys.argv))
    drmaa_session = None
    try:
        # Set up the DRMAA session for running cluster jobs
        import drmaa
        drmaa_session = drmaa.Session()
        drmaa_session.initialize()
    except Exception as e:
        print("{progname} error using DRMAA library".format(progname=program_name), file=sys.stdout)
        print("Error message: {msg}".format(msg=e.message, file=sys.stdout))
        exit(error_codes.DRMAA_ERROR)
    # Parse the configuration file, and initialise global state
    config = Config(options.config)
    config.validate()
    state = State(options=options, config=config, logger=logger,
                  drmaa_session=drmaa_session)
    # Build the pipeline workflow
    pipeline = make_pipeline(state)
    # Run (or print) the pipeline
    cmdline.run(options)
    if drmaa_session is not None:
        # Shut down the DRMAA session
        drmaa_session.exit()


if __name__ == '__main__':
    main()
