#!/usr/bin/env python3
from __future__ import print_function
import sys
import argparse
import time
import os
import logging
import json
import pprint
import subprocess
import glob
import shlex
from util.validate_requirements import validate_requirements
from util.helper_functions import end_program, sanitize_path, generate_hash,\
    checkpoint_reached, make_blast_database, json_writer, check_if_fasta,\
    remove_update_tracker
from util.blast_parser import read_blast_results, print_blast_summary,\
    print_fasta_files
from multiprocessing import Pool
from collections import defaultdict
from shutil import copy
from Bio import SeqIO
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

__version__ = '2.2.0'

# ARGPARSE SECTION ############################################################


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{} does not exist".format(x))
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

    return logger


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
    # parser.add_argument(
    #     '-pb', '--printblast', help='Print BLAST searches and exit. For '
    #                                 'submitting jobs to clusters.',
    #     action='store_true', default=False)
    parser.add_argument(
        '--outgroup', help='Name of outgroup file or strain to root on.',
        type=str)
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
    args.logger = init_logger(args)
    args.logger.debug('Started autoMLSA.py for run: {}'.format(args.runid))
    args.logger.debug(pprint.pformat(vars(args)))
    return args

# END ARGPARSE SECTION ########################################################
# CONFIGURATION SECTION #######################################################


def get_fasta_list(dirpath, logger, rundir):
    """
    Finds FASTA files from directory
    """
    if not os.path.exists(dirpath):
        fastadir = os.path.join(rundir, 'fasta')
        if os.path.exists(fastadir):
            msg = 'Original genome directory ({}) does not exist.'
            logger.warning(msg.format(dirpath))
            msg = 'Attempting to run from already renamed FASTA files.'
            logger.warning(msg)
            dirpath = fastadir
    fasta_list = []
    for fa in glob.iglob('{}/*'.format(dirpath)):
        # sys.stderr.write('Checking if {} is FASTA.\n'.format(fa))
        if os.path.isdir(fa):
            continue
        if check_if_fasta(fa):
            # sys.stderr.write('{} is FASTA.\n'.format(fa))
            fasta_list.append(os.path.abspath(fa))
        else:
            msg = '{} does not appear to be FASTA file, skipping.'
            logger.debug(msg.format(fa))
    if fasta_list:
        msg = 'Found FASTA files in {}'
        logger.info(msg.format(dirpath))
    return fasta_list


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
                os.mkdir(os.path.join(rundir, '.autoMLSA', 'backups'))
                os.mkdir(os.path.join(rundir, '.autoMLSA', 'updated'))
                os.mkdir(os.path.join(rundir, '.autoMLSA', 'checkpoint'))
            except OSError as ose:
                msg = 'Unable to generate rundir "{}" : {}'
                sys.stderr.write(msg.format(rundir, ose) + '\n')
                sys.stderr.write('Check your path and try again.\n')
                exit(71)
    os.chdir(rundir)
    return rundir


def read_config(logger, configfile):
    """Finds settings in config file (json format). path/rundir/config.json

    input  - rundir path, logger
    return - dict with configuration settings
    """

    logger.info('Reading from configuration file: {}'.format(configfile))

    with open(configfile, 'r') as cf:
        config = json.load(cf)
    return config


def validate_arguments(args, config, checkpoint):
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
    error = False

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
    if args.outgroup is None:
        if 'outgroup' in config.keys() and config['outgroup'] is not None:
            args.outgroup = config['outgroup']
        else:
            args.outgroup = ''
    if args.checkpoint is None:
        if 'checkpoint' in config.keys() and config['checkpoint'] is not None:
            args.checkpoint = config['checkpoint']
        else:
            args.checkpoint = 'none'

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
        for f in config['files']:
            if f:
                args.files.append(f)

    if args.dir is None:
        args.dir = []
    if 'dir' in config.keys() and config['dir'] is not None:
        for d in config['dir']:
            if d:
                args.dir.append(d)

    if args.query is None:
        args.query = []
    if 'query' in config.keys() and config['query'] is not None:
        for q in config['query']:
            if q:
                args.query.append(q)
    for query in args.query.copy():
        if not check_if_fasta(query):
            msg = 'Specified query file {} does not appear to be FASTA file.'
            args.logger.error(msg.format(query))
            args.logger.debug('Removing {} from query list'.format(query))
            args.query.remove(query)
            error = True

    args.fasta = []
    # Combine lists
    for fastadir in args.dir:
        fastas = get_fasta_list(fastadir, args.logger, args.rundir)
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

    if not args.fasta:
        msg = 'No valid FASTA files provided as --dir or --files.'
        args.logger.error(msg)
        error = True

    args.logger.debug('Validated arguments:')
    args.logger.debug(pprint.pformat(vars(args)))

    write_config(args)
    if error:
        end_program(args.logger, -1)

    if checkpoint:
        checkpoint_reached(args.logger, 'at validate arguments')

    return args


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
    configdict.pop('checkpoint')
    json_writer(args.configfile, configdict)
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
                logger.debug('Adding {} to list of labels'.format(base))
                labels.append(base)
    else:
        logger.debug('Generating '.format(base))
        labels = bases
    json_writer(labelsf, labels)
    return(labels)

# END CONFIGURATION SECTION ###################################################
# FORMATTING SECTION ##########################################################


def convert_fasta(logger, rundir, fastas, labels, makeblastdb):
    """Converts fasta file with placeholder label names

    input - unformatted fasta files
    return - Set of formatted fasta files
    also generates BLAST dbs per fasta
    """
    fastadir = os.path.join(rundir, 'fasta')
    new_fastas = []
    renamef = os.path.join(rundir, '.autoMLSA', 'rename.json')
    deleted = False
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
        else:
            if check_hash(logger, fasta, base, rename_info[base]['hash']):
                deleted = True
        new_fastas.append(labelfastaf)
        if not os.path.exists(labelfastaf):
            msg = 'Writing renamed fasta for {}'
            logger.debug(msg.format(base))
            seq = ''
            with open(fasta, 'r') as ff, open(labelfastaf, 'w') as lff:
                for rec in SeqIO.parse(ff, 'fasta'):
                    lff.write('>{} {}\n'.format(label, rec.id))
                    lff.write('{}\n'.format(rec.seq))
                    seq += str(rec.seq)
            seqhash = generate_hash(seq)
            rename_info[base]['hash'] = seqhash
        make_blast_database(logger, makeblastdb, labelfastaf)
    json_writer(renamef, rename_info)

    expected_fastas_fn = os.path.join(rundir, '.autoMLSA',
                                      'expected_fastas.json')
    if os.path.exists(expected_fastas_fn):
        with open(expected_fastas_fn, 'r') as eff:
            expected_fastas = json.load(eff)
    else:
        expected_fastas = []

    updated = set(expected_fastas) != set(new_fastas)
    if updated or deleted:
        logger.debug('Found new genome files')
        open(os.path.join('.autoMLSA', 'updated', 'genome'), 'w').close()
    else:
        logger.debug('Found no new genome files')
    json_writer(expected_fastas_fn, new_fastas)

    return(new_fastas, updated or deleted)


def check_hash(logger, fasta, base, seqhash):
    """
    Checks hash, deletes files related to genome if necessary
    return - value if something is deleted
    """
    deleted = False
    with open(fasta, 'r') as ff:
        newseqhash = generate_hash(
            ''.join([str(rec.seq) for rec in SeqIO.parse(ff, 'fasta')]))
    if newseqhash != seqhash:
        deleted = True
        msg = 'Old genome file {} sequence has changed. Updating.'
        logger.info(msg.format(base))
        for fa in glob.iglob(os.path.join('fasta', base + '*')):
            os.remove(fa)
        for res in glob.iglob(os.path.join('blast', '*_' + base + '.tab')):
            os.remove(res)
    else:
        logger.debug('Genome file {} sequence has not changed.'.format(base))
    return deleted


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
        open(os.path.join('.autoMLSA', 'updated', 'query'), 'w').close()
        logger.debug('Found new query sequences')
    if not updated:
        logger.debug('Found no new query sequences')
    json_writer(expected_queries_fn, new_queries)

    if not new_queries:
        msg = 'No query sequences found. Check your query file and try again.'
        logger.error(msg)
        end_program(logger, -1)

    open(os.path.join('.autoMLSA', 'checkpoint', 'get_queries'), 'w').close()

    return(new_queries, updated)

# END FORMATTING SECTION ######################################################
# BLAST SECTION ###############################################################


def generate_blast_list(logger, rundir, exe, queries,
                        targets, evalue, threads, checkpoint):
    """
    Generates list of BLAST searches that need to be run

    input  - blast type, exes, queries, and targets
    output - list of blast output names
    """
    cmds = []
    # blastout = defaultdict(list)
    blastout = []
    blastdir = os.path.join(rundir, 'blast')
    outfmt = '7 qseqid sseqid saccver pident qlen length bitscore qcovhsp ' \
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
            blastout.append(os.path.join(outpath))
            cmd = base_cmd + ['-db', target,
                              '-query', query,
                              '-out', outpath]
            if not os.path.exists(outpath) or os.path.getsize(outpath) == 0:
                cmds.append(cmd)

    if cmds:
        json_writer(os.path.join('.autoMLSA', 'new_blasts.json'))
        if checkpoint:
            blastfile = os.path.join(rundir, 'blastcmds.txt')
            with open(blastfile, 'w') as fh:
                for cmd in cmds:
                    fh.write(' '.join([shlex.quote(x) for x in cmd]) + '\n')
            msg = 'BLAST commands written to {}. Exiting.'
            logger.info(msg.format(blastfile))
            checkpoint_reached(logger, 'prior to BLAST searches')
        else:
            msg = 'Running {} BLAST searches using {} CPUs.'
            logger.info(msg.format(len(cmds), threads))
            p = Pool(threads)
            with p:
                p.map(subprocess.run, cmds)
    else:
        logger.info('No BLAST searches remaining. Moving to parse.')

    open(os.path.join('.autoMLSA', 'checkpoint', 'generate_blast_list'),
         'w').close()
    return(blastout)

# END BLAST SECTION ###########################################################
# ALIGNMENT SECTION ###########################################################


def run_mafft(logger, threads, mafft, unaligned, updated, checkpoint):
    """
    input  - unaligned fasta files per query
    return - list of aligned files per query
    """
    base_cmd = [mafft, '--thread', str(threads)]
    aligneddir = 'aligned'
    aligned = []
    cmdstrs = []
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
            cmdstr = ' '.join([shlex.quote(x) for x in cmd]) + \
                ' > {}'.format(shlex.quote(outname))
            logger.debug(cmdstr)
            if checkpoint == 'prealign':
                cmdstrs.append(cmdstr)
            else:
                with open(outname, 'w') as fh,\
                        open(logname, 'w') as logfh:
                    subprocess.run(cmd, stdout=fh, stderr=logfh, text=True)
    if cmdstrs:
        with open('mafftcmds.txt', 'w') as fh:
            fh.write('\n'.join(cmdstrs))
        logger.info('MAFFT alignment commands written to mafftcmds.txt.')
        logger.info('Run these commands and resubmit to continue.')
        checkpoint_reached(logger, 'prior to MAFFT alignments')
    elif checkpoint == 'postalign':
        checkpoint_reached(logger, 'after MAFFT alignments')

    checkpath = os.path.join('.autoMLSA', 'checkpoint', 'run_mafft')
    if not os.path.exists(checkpath):
        remove_update_tracker(['genome', 'filt'])
        open(os.path.join(checkpath, 'w')).close()

    return(aligned)

# END ALIGNMENT SECTION #######################################################
# PHYLOGENY SECTION ###########################################################


def generate_nexus(logger, runid, aligned, updated, checkpoint):
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
            msg = 'Nexus file is being overwritten due to updated query.'
            logger.debug(msg)
        with open(nexus, 'w') as nex:
            nex.write('#nexus\n')
            nex.write('begin sets;\n')
            for fn in aligned:
                base = os.path.basename(os.path.splitext(fn)[0])
                nex.write('\tcharset {} = {}: *;\n'.format(base, fn))
            nex.write('end;\n')
    if checkpoint:
        checkpoint_reached(logger,
                           'after generating nexus file {}'.format(nexus))

    open(os.path.join('.autoMLSA', 'checkpoint', 'generate_nexus'),
         'w').close()
    checkpath = os.path.join('.autoMLSA', 'checkpoint', 'run_mafft')
    if not os.path.exists(checkpath):
        remove_update_tracker(['query'])
        open(os.path.join(checkpath, 'w')).close()
    return(nexus)


def run_iqtree(logger, threads, iqtree, nexus, updated, outgroup):
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
    if outgroup:
        cmd.extend(['-o', outgroup])
    if os.path.exists(out_tree) and not updated:
        logger.info('Treefile {} already found. Skipping iqtree.'
                    .format(out_tree))
    else:
        if os.path.exists(out_tree):
            logger.info('Re-doing iqtree search.')
            cmd.append('-redo')
        logger.info('Running: {}'.format(' '.join(cmd)))
        subprocess.run(cmd)
        if not os.path.exists(out_tree):
            msg = 'iqtree2 seems to have failed.'
            logger.critical(msg)
            msg = 'Check the log files for error messages to see if they can '\
                'be resolved.'
            logger.critical(msg)
            end_program(logger, 73)

    checkpath = os.path.join('.autoMLSA', 'checkpoint', 'run_iqtree')
    if not os.path.exists(checkpath):
        remove_update_tracker(['query', 'genome', 'filt'])
        open(os.path.join(checkpath, 'w')).close()
    return out_tree


def exit_successfully(logger, rundir, treefile):
    """Temporary command

    input  - rundir and treefile to print to log
    return - None
    """
    remove_update_tracker(['query', 'genome', 'filt'])

    msg = 'Your treefile is ready: {}/{}'
    logger.info(msg.format(rundir, treefile))
    end_program(logger, 0)

# END PHYLOGENY SECTION #######################################################


def main():
    # ARGPARSE SECTION
    args = run_argparse()
    args.logger.info('Welcome to autoMLSA.py version {}'.format(__version__))
    time.sleep(1)

    # CONFIGURATION SECTION
    exes = validate_requirements(args.logger)

    args.logger.info('Reconciling configuration settings.')
    config = {}
    if os.path.exists(args.configfile):
        config = read_config(args.logger, args.configfile)
        if args.config is not None:
            msg = 'Config file specified {} is ignored as one is found in {}.'
            args.logger.warning(msg.format(args.config, args.rundir))
    elif args.config:
        config = read_config(args.logger, args.config)
    args = validate_arguments(args, config, args.checkpoint == 'validate')

    # FORMATTING SECTION
    args.logger.info('Converting genome FASTA files for BLAST if necessary.')
    labels = get_labels(args.logger, args.rundir, args.fasta)
    newfastas, updated_genome = convert_fasta(
        args.logger, args.rundir, args.fasta, labels, exes['makeblastdb'])

    args.logger.info('Extracting query FASTA files if necessary.')
    queries, updated_query = get_queries(
        args.logger, args.rundir, args.dups, args.query)

    # BLAST SECTION
    args.logger.info('Generating list of BLAST searches and outputs.')
    blastout = generate_blast_list(
        args.logger, args.rundir, exes[args.program], queries, newfastas,
        args.evalue, args.threads, args.checkpoint == 'preblast')

    # BLAST output results, summary, and files
    blastres = read_blast_results(
        args.logger, blastout, args.coverage, args.identity, updated_genome)
    blastfilt, updated_filt = print_blast_summary(
        args.logger, blastres, labels, args.allow_missing, args.missing_check,
        args.checkpoint == 'filtering')
    unaligned = print_fasta_files(
        args.logger, blastfilt, labels, updated_genome or updated_filt)

    # ALIGNMENT SECTION
    aligned = run_mafft(
        args.logger, args.threads, exes['mafft'], unaligned,
        updated_genome or updated_filt, args.checkpoint)

    # PHYLOGENY SECTION
    nexus = generate_nexus(
        args.logger, args.runid, aligned, updated_query,
        args.checkpoint == 'nexus')
    treefile = run_iqtree(
        args.logger, args.threads, exes['iqtree'], nexus,
        updated_query or updated_genome or updated_filt,
        args.outgroup)
    exit_successfully(args.logger, args.rundir, treefile)


if __name__ == '__main__':
    main()
