#!/usr/bin/env python3
from __future__ import print_function
import sys
import argparse
import time
# import csv
import os
import logging
import json
import pprint
import subprocess
import re
import glob
import shlex
from multiprocessing import Pool
from hashlib import blake2b
from collections import defaultdict
from shutil import copy
import util.blast_parser as bp
# from timeit import default_timer as timer
from Bio import SeqIO
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

__version__ = '2.2.0'


def validate_requirements(logger):
    """Identifies installed software and sets paths for running

    input  - logger
    output - exes dict with paths to executables
    """

    # Set local path to blast executable DIR
    BLASTPATH = ''

    # Expects in $PATH, can hard-code to full path here
    tblastn = 'tblastn'
    blastn = 'blastn'
    makeblastdb = 'makeblastdb'
    mafft = 'mafft-linsi'
    iqtree = 'iqtree2'

    if os.path.exists(os.path.join(BLASTPATH, tblastn)):
        tblastn = os.path.join(BLASTPATH, tblastn)
    if os.path.exists(os.path.join(BLASTPATH, blastn)):
        blastn = os.path.join(BLASTPATH, blastn)
    if os.path.exists(os.path.join(BLASTPATH, makeblastdb)):
        makeblastdb = os.path.join(BLASTPATH, makeblastdb)

    try:
        tblastn_ver = subprocess.check_output([tblastn, '-version'],
                                              stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find tblastn executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(tblastn))
        msg = 'Please add tblastn to your $PATH or supply a BLASTPATH in {}'
        logger.error(msg.format(__file__))
        end_program(logger, 78)
    else:
        logger.debug('tblastn found: {}'.format(tblastn_ver.strip().decode()))

    try:
        blastn_ver = subprocess.check_output([blastn, '-version'],
                                             stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find blastn executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(blastn))
        msg = 'Please add blastn to your $PATH or supply a BLASTPATH in {}'
        logger.error(msg.format(__file__))
        end_program(logger, 78)
    else:
        logger.debug('blastn found: {}'.format(blastn_ver.strip().decode()))

    try:
        makeblastdb_ver = subprocess.check_output([makeblastdb, '-version'],
                                                  stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find makeblastdb executable in path or in provided '\
              'dir ({})'
        logger.error(msg.format(makeblastdb))
        msg = 'Please add makeblastdb to your $PATH or supply a BLASTPATH in' \
            '{}'
        logger.error(msg.format(__file__))
        end_program(logger, 78)
    else:
        logger.debug('makeblastdb found: {}'.format(
            makeblastdb_ver.strip().decode()))

    # MAFFT testing
    try:
        mafft_ver = subprocess.check_output([mafft, '--version'],
                                            stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find mafft executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(mafft))
        msg = 'Please add mafft to your $PATH or supply a full mafft path in '\
              '{}'
        logger.error(msg.format(__file__))
        end_program(logger, 78)
    else:
        logger.debug('mafft found: {}'.format(mafft_ver.strip().decode()))

    # iqtree testing
    try:
        iqtree_ver = subprocess.check_output([iqtree, '--version'],
                                             stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find iqtree2 executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(iqtree))
        msg = 'Please add iqtree2 to your $PATH or supply a full iqtree path '\
            ' in {}'
        logger.error(msg.format(__file__))
        end_program(logger, 78)
    else:
        logger.debug('iqtree2 found: {}'.format(iqtree_ver.strip().decode()))

    exes = {}
    exes['tblastn'] = tblastn
    exes['blastn'] = blastn
    exes['makeblastdb'] = makeblastdb
    exes['mafft'] = mafft
    exes['iqtree'] = iqtree
    return(exes)


def end_program(logger, code):
    """
    Program message including success or failure of the program.
    """
    if code:
        msg = 'Program exiting with code ({}) indicating failure.'
        logger.info(msg.format(code))
        msg = 'Check error messages to resolve the problem.'
        logger.info(msg)
        exit(code)
    else:
        msg = 'Program was a success! Congratulations!'
        logger.info(msg)
        exit(code)


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return os.path.abspath(x)


def get_fasta_files(x):
    """
    'Type' for argparse - checks if file is FASTA format
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    elif check_if_fasta(x):
        return os.path.abspath(x)
    else:
        raise argparse.ArgumentTypeError("{0} is not a FASTA file".format(x))


def get_fasta_dir(x, logger, rundir):
    """
    Finds FASTA files from directory
    """
    if not os.path.exists(x):
        fastadir = os.path.join(rundir, 'fasta')
        if os.path.exists(fastadir):
            msg = 'Original genome directory ({}) does not exist.'
            logger.warning(msg.format(x))
            msg = 'Attempting to run from already renamed FASTA files.'
            logger.warning(msg)
            x = fastadir
        else:
            raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    fasta_list = []
    for fa in glob.iglob('{}/*'.format(x)):
        # sys.stderr.write('Checking if {} is FASTA.\n'.format(fa))
        if os.path.isdir(fa):
            continue
        if check_if_fasta(fa):
            # sys.stderr.write('{} is FASTA.\n'.format(fa))
            fasta_list.append(os.path.abspath(fa))
    return fasta_list


def check_if_fasta(fa):
    """
    Checks first line of file for '>' to determine if FASTA, returns bool
    """
    with open(fa, 'r') as fah:
        if fah.readline().startswith('>'):
            return True
        else:
            return False


def validate_email(x):
    """
    'Type' for argparse - checks that email appears somewhat valid
    """
    if not re.fullmatch(r"[^@]+@[^@]+\.[^@]+", x):
        msg = '{} does not appear to be a valid email.'
        raise argparse.ArgumentTypeError(msg.format(x))
    return x


def sanitize_path(s):
    """
    Turns string into filename without symbols
    """
    s = s.replace(' ', '_')
    filename = ''.join(x for x in s if (x.isalnum() or x in '._-'))
    return filename


def generate_hash(s):
    """
    Generates blake2b hash for sequence
    input  - nucleotide sequence
    return - blake2b hash digest
    """
    seqhash = blake2b(s.encode(), digest_size=16)
    return seqhash.hexdigest()


def make_blast_database(logger, makeblastdb, fasta):
    """Generates blast DB if necessary per fasta

    input  - fasta file
    return - null
    """
    suffixes = ['.nsq', '.nin', '.nhr', '.nto', '.not', '.ndb', '.ntf']
    exist = [os.path.exists(fasta + x) for x in suffixes]
    if not all(exist):
        subprocess.run([makeblastdb, '-dbtype', 'nucl', '-in', fasta])


def init_logger(args):
    """Sets up logging system to file and stderr

    input  - parsed arguments
    output - logger object
    """

    # initialize logger with program name
    # must set logging level to DEBUG to print DEBUG to file
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(name)s - %(asctime)s - '
                                  '%(levelname)s - %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')

    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(formatter)
    if args.debug:
        stderr_handler.setLevel(logging.DEBUG)
    elif args.quiet:
        stderr_handler.setLevel(logging.WARNING)
    else:
        stderr_handler.setLevel(logging.INFO)

    file_handler = logging.FileHandler('{}/{}.log'.format(args.rundir,
                                                          args.runid))
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(stderr_handler)
    logger.addHandler(file_handler)

    return(logger)


def run_argparse():
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
    parser.add_argument(
        '-pb', '--printblast', help='Print BLAST searches and exit. For '
                                    'submitting jobs to clusters.',
        action='store_true', default=False)
    parser.add_argument(
        '--email', help='E-mail is required for downloading data from '
                        'NCBI. You can also set the $EMAIL environment '
                        'variable.',
        type=validate_email)
    parser.add_argument(
        '--debug', help='Turn on debugging messages.', action='store_true',
        default=False)
    parser.add_argument(
        '--quiet', help='Turn off progress messages.', action='store_true',
        default=False)
    args = parser.parse_args()
    args.rundir = config_rundir(args.runid)
    args.configfile = os.path.join(args.rundir, 'config.json')
    args.logger = init_logger(args)
    args.logger.debug('Started autoMLSA.py for run: {}'.format(args.runid))
    args.logger.debug(pprint.pformat(vars(args)))
    return args


def config_rundir(runid):
    """Identifies existing (or not) run directory from input runid

    input  - runid as given as input
    output - full path to the run directory
    """
    if os.path.exists('../{}'.format(runid)):
        rundir = os.path.abspath('../{}'.format(runid))
    else:
        rundir = os.path.abspath('./{}'.format(runid))
        if not os.path.exists(rundir):
            try:
                os.mkdir(rundir)
                os.mkdir(os.path.join(rundir, '.autoMLSA'))
            except OSError as ose:
                msg = 'Unable to generate rundir "{}" : {}'
                sys.stderr.write(msg.format(rundir, ose) + '\n')
                sys.stderr.write('Check your path and try again.\n')
                exit(71)
    os.chdir(rundir)
    return(rundir)


def read_config(logger, configfile):
    """Finds settings in config file (json format). path/rundir/config.json

    input  - rundir path, logger
    return - dict with configuration settings
    """

    logger.info('Reading from configuration file: {}'.format(configfile))

    with open(configfile, 'r') as cf:
        config = json.load(cf)
    return(config)


def validate_arguments(args, config):
    """Checks for executables as well as makes sure files exist.

    input  - args from argparse; config from config file
    return - validated args from argparse
    """
    # Determine if defaults need to be set
    # e = 1e-5
    # m = 500
    # c = 50
    # p = tblastn
    # t = 1

    args.logger.debug('Validating & reconciling arguments.')

    # Reconcile config file and command line options
    if args.evalue is None:
        if 'evalue' in config.keys() and config['evalue'] is not None:
            args.evalue = float(config['evalue'])
        else:
            args.evalue = 1e-5
    # if args.target is None:
    #     if 'target' in config.keys() and config['target'] is not None:
    #         args.target = int(config['target'])
    #     else:
    #         args.target = 500
    if args.coverage is None:
        if 'coverage' in config.keys() and config['coverage'] is not None:
            args.coverage = int(config['coverage'])
        else:
            args.coverage = 50
    if args.identity is None:
        if 'identity' in config.keys() and config['identity'] is not None:
            args.identity = int(config['identity'])
        else:
            args.identity = 30
    if args.allow_missing is None:
        if 'allow_missing' in config.keys() and \
                config['allow_missing'] is not None:
            args.allow_missing = int(config['allow_missing'])
        else:
            args.allow_missing = 0
    if not args.missing_check:
        if 'missing_check' in config.keys() and \
                config['missing_check']:
            args.allow_missing = True
        else:
            args.missing_check = False
    if args.program is None:
        if 'program' in config.keys():
            if config['program'] in ('tblastn', 'blastn'):
                args.program = config['program']
            else:
                msg = 'Program specified in config file "{}" is not valid.'
                args.logger.error(msg.format(config['program']))
                msg = 'Please give either "tblastn" or "blastn" and try again.'
                args.logger.error(msg)
                end_program(args.logger, 78)
        else:
            args.program = 'tblastn'
    if args.threads is None:
        if 'threads' in config.keys() and config['threads'] is not None:
            args.threads = int(config['threads'])
        else:
            args.threads = 1
    if args.email is None:
        try:
            args.email = os.environ['EMAIL']
        except KeyError:
            msg = 'No email was provided. Give one using --email or set ' \
                  '$EMAIL as an environment variable in your shell.'
            args.logger.error(msg)
            end_program(args.logger, 78)

    if args.evalue > 10:
        msg = 'Specified evalue "{}" is greater than 10.'
        args.logger.error(msg.format(args.evalue))
        msg = 'Please specify an evalue < 10 and try again.'
        args.logger.error(msg)
        end_program(args.logger, 78)
    # if args.target not in range(1, 20000):
    #     msg = 'Specified target value "{}" is outside of range 1-20000.'
    #     args.logger.error(msg.format(args.target))
    #     msg = 'Please specify a target between 1 and 20000 and try again.'
    #     args.logger.error(msg)
    #     end_program(args.logger, 78)
    if args.coverage not in range(0, 100):
        msg = 'Coverage value is not between 0 and 100.'
        args.logger.error(msg.format(args.coverage))
        msg = 'Please specify a coverage between 0 and 100 and try again.'
        args.logger.error(msg)
        end_program(args.logger, 78)

    cmd = 'lscpu | grep -G "^CPU(s):" | grep -o -E "[0-9]+"'
    try:
        maxthreads = subprocess.check_output(cmd,
                                             shell=True,
                                             stderr=subprocess.STDOUT
                                             ).decode().strip()
    except subprocess.CalledProcessError:
        msg = 'Unable to check number of available threads.'
        args.logger.warning(msg)
        msg = 'Make sure the number of specified threads is correct.'
        args.logger.warning(msg)
    else:
        args.logger.debug('Max threads found {}'.format(maxthreads))
        if args.threads > int(maxthreads):
            msg = 'Threads specified {} greater than number of available ' \
                  'threads {}'
            args.logger.error(msg.format(args.threads, maxthreads))
            msg = 'Specify threads less than or equal to {} and try again.'
            args.logger.error(msg.format(maxthreads))
            end_program(args.logger, 78)

    # Set up lists for files, dirs, and queries
    if args.files is None:
        args.files = []
    if 'files' in config.keys() and config['files'] is not None:
        args.files = config['files'] + args.files

    if args.dir is None:
        args.dir = []
    if 'dir' in config.keys() and config['dir'] is not None:
        args.dir = config['dir'] + args.dir

    if args.query is None:
        args.query = []
    if 'query' in config.keys() and config['query'] is not None:
        args.query = config['query'] + args.query

    args.fasta = []
    # Combine lists
    for fastadir in args.dir:
        fastas = get_fasta_dir(fastadir, args.logger, args.rundir)
        args.fasta.extend(fastas)
    args.fasta.extend(args.files)

    # Remove duplicates
    args.fasta = list(dict.fromkeys(args.fasta))
    args.dir = list(dict.fromkeys(args.dir))
    args.query = list(dict.fromkeys(args.query))

    seen = {}
    for fasta in args.fasta:
        base = os.path.basename(fasta)
        if base in seen:
            msg = 'Same genome name found in two locations'
            args.logger.error(msg)
            msg = '++++++++++++Offending files++++++++++++'
            args.logger.error(msg)
            args.logger.error('{}'.format(fasta))
            args.logger.error('{}'.format(seen[base]))
            msg = '+++++++++++++++++++++++++++++++++++++++'
            args.logger.error(msg)
            msg = 'Remove or rename one of these to continue.'
            args.logger.error(msg)
            end_program(args.logger, 65)
        else:
            seen[base] = fasta

    if 'dups' in config.keys() and config['dups']:
        args.dups = True

    args.logger.debug('Validated arguments:')
    args.logger.debug(pprint.pformat(vars(args)))

    return(args)


def write_config(args):
    """Writes new config file, overwriting old if necessary

    input  - valiated & compiled arguments
    return - NULL
    """
    args.logger.debug('Writing config file {}'.format(args.configfile))
    configdict = vars(args).copy()
    configdict.pop('config')
    configdict.pop('debug')
    configdict.pop('logger')
    configdict.pop('quiet')
    configdict.pop('configfile')
    configdict.pop('rundir')
    bp.json_writer(args.configfile, configdict)
    # with open(args.configfile, 'w') as cf:
    #     json.dump(configdict, cf, indent=4)
    #     cf.write('\n')


def get_labels(logger, rundir, fastas):
    """Determines the arbitrary labels per fasta file
    labels are index of list

    input - fasta files
    return - label list
    """
    labelsf = os.path.join(rundir, '.autoMLSA', 'labels.json')
    bases = [os.path.basename(x) for x in fastas]
    if os.path.exists(labelsf):
        with open(labelsf, 'r') as lf:
            labels = json.load(lf)
        for base in bases:
            if base not in labels:
                # sys.stderr.write(
                #     'Base {} not found in bases list'.format(base))
                labels.append(base)
    else:
        labels = bases
    bp.json_writer(labelsf, labels)
    # with open(labelsf, 'w') as lf:
    #     json.dump(labels, lf, indent=4)
    #     lf.write('\n')
    return(labels)


def convert_fasta(logger, rundir, fastas, labels, makeblastdb):
    """Converts fasta file with placeholder label names

    input - unformatted fasta files
    return - Set of formatted fasta files
    also generates BLAST dbs per fasta
    """
    fastadir = os.path.join(rundir, 'fasta')
    new_fastas = []
    renamef = os.path.join(rundir, '.autoMLSA', 'rename.json')
    if os.path.exists(renamef):
        with open(renamef, 'r') as rf:
            rename_info = json.load(rf)
    else:
        rename_info = defaultdict(dict)
    if not os.path.exists(fastadir):
        os.mkdir(fastadir)
    for fasta in fastas:
        base = os.path.basename(fasta)
        label = labels.index(base)
        labelfastaf = os.path.join(fastadir, base)
        if base not in rename_info:
            try:
                rename_info[base]['index'] = label
            except KeyError:
                rename_info[base] = {}
                rename_info[base]['index'] = label
        new_fastas.append(labelfastaf)
        if not os.path.exists(labelfastaf):
            msg = 'Writing renamed fasta for {}'
            logger.debug(msg.format(base))
            with open(fasta, 'r') as ff, open(labelfastaf, 'w') as lff:
                for rec in SeqIO.parse(ff, 'fasta'):
                    lff.write('>{} {}\n'.format(label, rec.id))
                    lff.write('{}\n'.format(rec.seq))
                    if 'recid' not in rename_info[base]:
                        rename_info[base]['recid'] = rec.id
        make_blast_database(logger, makeblastdb, labelfastaf)
    bp.json_writer(renamef, rename_info)

    expected_fastas_fn = os.path.join(rundir, '.autoMLSA',
                                      'expected_fastas.json')
    if os.path.exists(expected_fastas_fn):
        with open(expected_fastas_fn, 'r') as eff:
            expected_fastas = json.load(eff)
    else:
        expected_fastas = []

    updated = set(expected_fastas) != set(new_fastas)
    if updated:
        logger.debug('Found new genome files')
    if not updated:
        logger.debug('Found no new genome files')
    bp.json_writer(expected_fastas_fn, new_fastas)

    return(new_fastas, updated)


def get_queries(logger, rundir, dups, queries):
    """Converts query fasta file(s) into individual fastas for BLAST

    input - FASTA files with queries as individual entries
    return - individual fasta files, bool showing changes in queries
    """
    querydir = os.path.join(rundir, 'queries')
    backups = os.path.join(rundir, '.autoMLSA', 'backups')

    new_queries = []
    seen = {}
    hashes = {}
    if not os.path.exists(querydir):
        os.mkdir(querydir)
    if not os.path.exists(backups):
        os.mkdir(backups)
    for query_file in queries:
        i = 1
        query_base = os.path.basename(query_file)
        query_backup = os.path.join(backups, query_base)
        if not os.path.exists(query_file):
            if os.path.exists(query_backup):
                msg = 'Using query backup file as {} is missing.'
                logger.warning(msg.format(query_file))
                query_file = query_backup
            else:
                msg = 'Query file {} seems to have been removed and/or lost.'
                logger.critical(msg.format(query_file))
                msg = 'No backup was found in the backups directory.'
                logger.critical(msg)
                msg = 'Unable to continue without the queries. Either replace'\
                    ' the file or start again from a new analysis.'
                logger.critical(msg)
                end_program(logger, 66)
        with open(query_file, 'r') as qf:
            for rec in SeqIO.parse(qf, 'fasta'):
                safeid = sanitize_path(rec.id)
                seqhash = generate_hash(str(rec.seq))

                if seqhash in hashes and hashes[seqhash] == safeid:
                    # Skip duplicate header + seq
                    continue
                elif seqhash in hashes and hashes[seqhash] != safeid:
                    # Same seq but different names, throw error
                    mismatchid = hashes[seqhash]
                    msg = 'Same sequence found in two query inputs with '\
                        'different sequence IDs:'
                    logger.error(msg)
                    msg = '+++++++++++++++Offending queries+++++++++++++++'
                    logger.error(msg)
                    if query_file == seen[mismatchid]['path']:
                        msg = '{} - TWICE - header ID {} & {}'
                        logger.error(
                            msg.format(os.path.basename(query_file),
                                       seen[mismatchid]['unsafeid'],
                                       rec.id))
                    else:
                        msg = '{} header ID {}'
                        logger.error(
                            msg.format(os.path.basename(query_file), rec.id))
                        logger.error(
                            msg.format(
                                os.path.basname(seen[mismatchid]['path']),
                                seen[mismatchid]['unsafeid']))
                    msg = '+++++++++++++++++++++++++++++++++++++++++++++++'
                    logger.error(msg)
                    msg = 'Check your sequences to make sure they aren\'t '\
                        'misnamed, fix the problem, and try again.'
                    logger.error(msg)
                    end_program(logger, 65)
                else:
                    # Expected outcome; new id, new seq
                    hashes[seqhash] = safeid

                # Check for same id, allow if --dups flag, otherwise error
                if safeid in seen and not dups:
                    msg = 'Same query name ({}) found in two query inputs:'
                    logger.error(msg.format(safeid))
                    msg = '+++++++++++++++Offending queries+++++++++++++++'
                    logger.error(msg)
                    if query_file == seen[safeid]['path']:
                        msg = '{} - TWICE - sequences {} & {}'
                        logger.error(
                            msg.format(os.path.basename(query_file),
                                       seen[safeid]['index'], i))
                    else:
                        msg = '{} sequence {}'
                        logger.error(
                            msg.format(os.path.basename(query_file), i))
                        logger.error(
                            msg.format(os.path.basename(seen[safeid]['path']),
                                       seen[safeid]['index']))
                    msg = '+++++++++++++++++++++++++++++++++++++++++++++++'
                    logger.error(msg)
                    msg = 'Remove or rename one of these to continue.'
                    logger.error(msg)
                    msg = 'Alternatively, if this is intended, use the --dups'\
                        ' flag to include both copies.'
                    logger.error(msg)
                    end_program(logger, 65)
                else:
                    if safeid in seen and dups:
                        msg = 'Keeping additional query {} with duplicate ID '\
                            'from file {}.'
                        logger.info(msg.format(safeid,
                                               os.path.basename(query_file)))
                    seen[safeid] = {}
                    seen[safeid]['path'] = query_file
                    seen[safeid]['index'] = i
                    seen[safeid]['unsafeid'] = rec.id

                fn = os.path.join(querydir, safeid + '_' + seqhash + '.fas')
                new_queries.append(fn)
                if not os.path.exists(fn):
                    msg = 'Writing {} for seq {} in {}'
                    logger.debug(msg.format(os.path.basename(fn), i,
                                            os.path.basename(query_file)))
                    with open(fn, 'w') as fh:
                        fh.write('>{}\n'.format(safeid))
                        fh.write('{}\n'.format(rec.seq))
                i += 1
        if query_file != query_backup:
            copy(query_file, backups)

    # Check if expected equals new_queries
    expected_queries_fn = os.path.join(rundir, '.autoMLSA',
                                       'expected_queries.json')
    if os.path.exists(expected_queries_fn):
        with open(expected_queries_fn, 'r') as eqf:
            expected_queries = json.load(eqf)
    else:
        expected_queries = []

    updated = set(expected_queries) != set(new_queries)
    if updated:
        logger.debug('Found new query sequences')
    if not updated:
        logger.debug('Found no new query sequences')
    bp.json_writer(expected_queries_fn, new_queries)

    return(new_queries, updated)


def generate_blast_list(logger, rundir, exe, queries,
                        targets, evalue, threads, pb):
    """
    Generates list of BLAST searches that need to be run

    input  - blast type, exes, queries, and targets
    output - list of blast output names
    """
    cmds = []
    # blastout = defaultdict(list)
    blastout = []
    blastdir = os.path.join(rundir, 'blast')
    outfmt = '7 qseqid sseqid saccver pident qlen length evalue qcovhsp '\
        'stitle sseq'
    base_cmd = [exe, '-evalue', str(evalue), '-outfmt', outfmt]
    cmd = []
    if not os.path.exists(blastdir):
        os.mkdir(blastdir)
    for target in targets:
        target_base = os.path.basename(target)
        for query in queries:
            query_base = os.path.splitext(os.path.basename(query))[0]
            outname = '{}_vs_{}.tab'.format(query_base, target_base)
            outpath = os.path.join(blastdir, outname)
            # blastout[target].append(os.path.join(outpath))
            blastout.append(os.path.join(outpath))
            cmd = base_cmd + ['-db', target,
                              '-query', query,
                              '-out', outpath]
            if os.path.exists(outpath):
                size = os.path.getsize(outpath)
                if size == 0:
                    cmds.append(cmd)
            else:
                cmds.append(cmd)

    if cmds:
        if pb:
            blastfile = os.path.join(rundir, 'blastcmds.txt')
            with open(blastfile, 'w') as fh:
                for cmd in cmds:
                    fh.write(shlex.join(cmd) + '\n')
            msg = 'BLAST commands written to {}. Exiting.'
            logger.info(msg.format(blastfile))
            exit(0)
        else:
            msg = 'Running {} BLAST searches using {} CPUs.'
            logger.info(msg.format(len(cmds), threads))
            p = Pool(threads)
            with p:
                p.map(run_blast_search, cmds)
    else:
        logger.info('No BLAST searches remaining. Moving to parse.')
    return(blastout)


def run_blast_search(cmd):
    subprocess.run(cmd)


def run_mafft(logger, threads, mafft, unaligned, updated):
    """
    input  - unaligned fasta files per query
    return - list of aligned files per query
    """
    base_cmd = [mafft, '--thread', str(threads)]
    aligneddir = 'aligned'
    aligned = []
    logger.info('Aligning FASTA sequences, if necessary.')
    if not os.path.exists(aligneddir):
        os.mkdir(aligneddir)
    for unalign in unaligned:
        outname = '{}/{}.aln'.format(
            aligneddir, os.path.basename(os.path.splitext(unalign)[0]))
        logname = '{}.log'.format(outname)
        logger.debug('Aligning and writing to {}'.format(outname))
        aligned.append(outname)
        if os.path.exists(outname) and not updated:
            msg = 'Aligned file {} already found, skipping.'
            logger.debug(msg.format(outname))
        else:
            if os.path.exists(outname):
                msg = 'Overwriting {} as genome was updated.'
                logger.debug(msg.format(outname))
            cmd = base_cmd + [unalign]
            logger.debug(shlex.join(cmd) + ' > {}'.format(outname))
            with open(outname, 'w') as fh,\
                    open(logname, 'w') as logfh:
                subprocess.run(cmd, stdout=fh, stderr=logfh, text=True)
    return(aligned)


def generate_nexus(logger, runid, aligned, updated):
    """
    Generate nexus file containing partitions per aligned gene

    input  - list of aligned FASTA filenames
    return - path to nexus file
    """
    # Example nexus file
    # #nexus
    # begin sets;
    #         charset part1 = rpoB.all.aln: *;
    #         charset part2 = ppsA.all.aln: *;
    #         charset part3 = gyrB.all.aln: *;
    #         charset part4 = dnaE.all.aln: *;
    #         charset part5 = pyrC.all.aln: *;
    #         charset part6 = aroE.all.aln: *;
    #         charset part7 = acsA.all.aln: *;
    #         charset part8 = recA.all.aln: *;
    # end;

    nexus = '{}.nex'.format(runid)
    if os.path.exists(nexus) and not updated:
        msg = 'Nexus file {} already found and no new queries were added.'
        logger.info(msg.format(nexus))
        msg = 'Skipping nexus file generation.'
        logger.info(msg)
    else:
        if os.path.exists(nexus):
            msg = 'Nexus file is being overwritten due to new query.'
            logger.debug(msg)
        with open(nexus, 'w') as nex:
            nex.write('#nexus\n')
            nex.write('begin sets;\n')
            for fn in aligned:
                base = os.path.basename(os.path.splitext(fn)[0])
                nex.write('\tcharset {} = {}: *;\n'.format(base, fn))
            nex.write('end;\n')
    return(nexus)


def run_iqtree(logger, threads, iqtree, nexus, updated):
    """
    Runs external iqtree command to generate phylogeny

    input  - nexus file
    return - path to output file
    """
    out_tree = '{}.treefile'.format(nexus)
    cmd = [iqtree,
           '-p', nexus,
           '-B', '1000',
           '-alrt', '1000',
           '-m', 'MFP+MERGE',
           '-rcluster', '10',
           '-nt', str(threads)]
    if os.path.exists(out_tree) and not updated:
        logger.info('Treefile {} already found. Skipping iqtree.'
                    .format(out_tree))
    else:
        if os.path.exists(out_tree):
            logger.info('Re-doing iqtree search.')
            cmd.append('-redo')
        logger.info('Running: {}'.format(shlex.join(cmd)))
        subprocess.run(cmd)
        if not os.path.exists(out_tree):
            msg = 'iqtree2 seems to have failed.'
            logger.critical(msg)
            msg = 'Check the log files for error messages to see if they can '\
                'be resolved.'
            logger.critical(msg)
            end_program(logger, 73)
    return out_tree


def exit_successfully(logger, rundir, treefile):
    """Temporary command

    input  - rundir and treefile to print to log
    return - None
    """
    msg = 'Your treefile is ready: {}/{}'
    logger.info(msg.format(rundir, treefile))
    end_program(logger, 0)


def main():
    args = run_argparse()
    args.logger.info('Welcome to autoMLSA.py version {}'.format(__version__))
    time.sleep(1)
    exes = validate_requirements(args.logger)

    # Reconcile configuration settings
    args.logger.info('Reconciling configuration settings.')
    config = {}
    if os.path.exists(args.configfile):
        config = read_config(args.logger, args.configfile)
        if args.config is not None:
            msg = 'Config file specified {} is ignored as one is found in {}.'
            args.logger.warning(msg.format(args.config, args.rundir))
    elif args.config:
        config = read_config(args.logger, args.config)
    args = validate_arguments(args, config)
    write_config(args)

    args.logger.info('Converting genome FASTA files for BLAST if necessary.')
    # Convert FASTA files
    labels = get_labels(args.logger, args.rundir, args.fasta)
    newfastas, updated_genome = convert_fasta(args.logger, args.rundir,
                                              args.fasta, labels,
                                              exes['makeblastdb'])

    args.logger.info('Extracting query FASTA files if necessary.')
    queries, updated_query = get_queries(args.logger, args.rundir, args.dups,
                                         args.query)

    args.logger.info('Generating list of BLAST searches and outputs.')
    blastout = generate_blast_list(args.logger, args.rundir,
                                   exes[args.program], queries, newfastas,
                                   args.evalue, args.threads, args.printblast)

    # BLAST output results, summary, and files
    blastres = bp.read_blast_results(args.logger, blastout, args.coverage,
                                     args.identity)
    blastfilt, updated_filt = bp.print_blast_summary(args.logger, blastres,
                                                     labels,
                                                     args.allow_missing,
                                                     args.missing_check)
    unaligned = bp.print_fasta_files(args.logger, blastfilt, labels,
                                     updated_genome or updated_filt)

    # Run external programs
    aligned = run_mafft(args.logger, args.threads, exes['mafft'], unaligned,
                        updated_genome or updated_filt)
    nexus = generate_nexus(args.logger, args.runid, aligned, updated_query)
    treefile = run_iqtree(args.logger, args.threads, exes['iqtree'], nexus,
                          updated_query or updated_genome or updated_filt)
    exit_successfully(args.logger, args.rundir, treefile)


if __name__ == '__main__':
    main()
