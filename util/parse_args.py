#!/usr/bin/env python
import os
import argparse
import logging
import pprint
import sys
from .helper_functions import check_if_fasta


def config_rundir(runid: str) -> str:
    """Identifies existing (or not) run directory from input runid

    input  - runid as given as input
    output - full path to the run directory
    """
    rundir: str = ''
    if os.path.exists('../{}'.format(runid)):
        rundir = os.path.abspath('../{}'.format(runid))
    else:
        rundir = os.path.abspath('./{}'.format(runid))
        if not os.path.exists(rundir):
            try:
                os.mkdir(rundir)
                os.mkdir(os.path.join(rundir, '.autoMLSA'))
                os.mkdir(os.path.join(rundir, '.autoMLSA', 'backups'))
                os.mkdir(os.path.join(rundir, '.autoMLSA', 'checkpoint'))
            except OSError as ose:
                msg = 'Unable to generate rundir "{}" : {}'
                sys.stderr.write(msg.format(rundir, ose) + '\n')
                sys.stderr.write('Check your path and try again.\n')
                exit(71)
    os.chdir(rundir)
    return rundir


def extant_file(x: str) -> str:
    """
    'Type' for argparse - checks that file exists
    input  - path to file
    return - absolute path to file
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{} does not exist".format(x))
    return os.path.abspath(x)


def get_fasta_files(x: str) -> str:
    """
    'Type' for argparse - checks if file is FASTA format
    input  - path to file
    return - absolute path to FASTA file
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    elif check_if_fasta(x):
        return os.path.abspath(x)
    else:
        raise argparse.ArgumentTypeError("{0} is not a FASTA file".format(x))


def init_logger(debug: bool, quiet: bool, rundir: str, runid: str) -> \
        None:
    """Sets up logging system to file and stderr

    input  - parsed arguments
    return - None; use logger = logging.getLogger(__name__)
    """

    # initialize logger with program name
    # must set logging level to DEBUG to print DEBUG to file
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(name)s - %(asctime)s - '
                                  '%(levelname)s - %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')

    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(formatter)
    if debug:
        stderr_handler.setLevel(logging.DEBUG)
    elif quiet:
        stderr_handler.setLevel(logging.WARNING)
    else:
        stderr_handler.setLevel(logging.INFO)

    # file_handler = logging.FileHandler('{}/{}.log'.format(rundir, runid))
    if rundir and runid:
        file_handler = logging.FileHandler(
            os.path.join(rundir, '{}.log'.format(runid)))
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    logger.addHandler(stderr_handler)


def run_argparse() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='This is a rewrite of autoMLSA.pl. Generates automated '
                    'multi-locus sequence analyses.')
    parser.add_argument(
        'runid', help='Name of the run directory.', type=str)
    parser.add_argument(
        '--query', help='Path to file with input seq(s).', type=extant_file,
        nargs='+')
    parser.add_argument(
        '--files', help='Path to the target genome FASTA files.',
        type=get_fasta_files, nargs='+')
    parser.add_argument(
        '--dir', help='Path to the target genome directory with FASTA files.',
        type=extant_file, nargs='+')
    parser.add_argument(
        '-e', '--evalue', help='E-value cutoff for BLAST searches. [1e-5]',
        type=float)
    parser.add_argument(
        '-c', '--coverage', help='Sets the coverage cut-off threshold. [50]',
        type=int)
    parser.add_argument(
        '-i', '--identity', help='Sets the identity cut-off threshold. [30]',
        type=int)
    parser.add_argument(
        '-p', '--program', help='Which BLAST program to run. [tblastn]',
        choices=['blastn', 'tblastn'], type=str)
    parser.add_argument(
        '--config', help='Path to configuration json file to copy.',
        type=extant_file)
    parser.add_argument(
        '--missing_check', help='Use this to confirm that settings have been '
                                'checked when genes are missing.',
        action='store_true', default=False)
    parser.add_argument(
        '-t', '--threads', help='Number of threads to use. [1]', type=int)
    parser.add_argument(
        '--dups', help='Allow for duplicate query names for more sequence '
                       'coverage across disparate organisms.',
        action='store_true', default=False)
    parser.add_argument(
        '--allow_missing', help='Allow for N missing genes per genome. [0]',
        type=int)
    # parser.add_argument(
    #     '-pb', '--printblast', help='Print BLAST searches and exit. For '
    #                                 'submitting jobs to clusters.',
    #     action='store_true', default=False)
    parser.add_argument(
        '--outgroup', help='Name of outgroup file or strain to root on.',
        type=str)
    parser.add_argument(
        '--protect', help='Save files from getting overwritten. By default, '
        'as input files update, older alignments and trees are deleted. '
        '*NOT CURRENTLY IMPLEMENTED*',
        action='store_true')
    parser.add_argument(
        '--checkpoint', help='Name of stage to stop computing on. [none]',
        type=str, choices=['validate', 'preblast', 'filtering', 'prealign',
                           'postalign', 'nexus', 'none'])
    # validate, newfastas, queries, blastout, blastfilt, alignment, nexus
    parser.add_argument(
        '--debug', help='Turn on debugging messages.', action='store_true',
        default=False)
    parser.add_argument(
        '--quiet', help='Turn off progress messages.', action='store_true',
        default=False)
    args = parser.parse_args()
    args.rundir = config_rundir(args.runid)
    args.configfile = os.path.join(args.rundir, 'config.json')
    init_logger(args.debug, args.quiet, args.rundir, args.runid)
    logger = logging.getLogger(__name__)
    logger.debug('Started autoMLSA.py for run: {}'.format(args.runid))
    logger.debug(pprint.pformat(vars(args)))
    return args