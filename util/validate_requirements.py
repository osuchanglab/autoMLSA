#!/usr/bin/env python3
from __future__ import print_function
import os.path
import subprocess
import logging
import platform
import urllib.request
import shutil
import inspect
from pathlib import Path as p
from hashlib import md5
from typing import Dict
from util.helper_functions import end_program
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)


def validate_requirements() -> Dict[str, str]:
    """Identifies installed software and sets paths for running

    input  - logger
    output - exes dict with paths to executables
    """
    logger = logging.getLogger('validate_requirements')

    # Set local path to blast executable DIR
    BLASTPATH: str = ''

    # Expects in $PATH, can hard-code to full path here
    tblastn: str = 'tblastn'
    blastn: str = 'blastn'
    makeblastdb: str = 'makeblastdb'
    mafft: str = 'mafft-linsi'
    iqtree: str = 'iqtree2'

    if os.path.exists(os.path.join(BLASTPATH, tblastn)):
        tblastn = os.path.join(BLASTPATH, tblastn)
    if os.path.exists(os.path.join(BLASTPATH, blastn)):
        blastn = os.path.join(BLASTPATH, blastn)
    if os.path.exists(os.path.join(BLASTPATH, makeblastdb)):
        makeblastdb = os.path.join(BLASTPATH, makeblastdb)

    try:
        tblastn_ver: str = subprocess.check_output([tblastn, '-version'],
                                                   stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find tblastn executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(tblastn))
        msg = 'Please add tblastn to your $PATH or supply a BLASTPATH in {}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.debug('tblastn found: {}'.format(tblastn_ver.strip().decode()))

    try:
        blastn_ver: str = subprocess.check_output([blastn, '-version'],
                                                  stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find blastn executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(blastn))
        msg = 'Please add blastn to your $PATH or supply a BLASTPATH in {}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.debug('blastn found: {}'.format(blastn_ver.strip().decode()))

    try:
        makeblastdb_ver: str = subprocess.check_output(
            [makeblastdb, '-version'],
            stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find makeblastdb executable in path or in provided '\
              'dir ({})'
        logger.error(msg.format(makeblastdb))
        msg = 'Please add makeblastdb to your $PATH or supply a BLASTPATH in' \
            '{}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.debug('makeblastdb found: {}'.format(
            makeblastdb_ver.strip().decode()))

    # MAFFT testing
    try:
        mafft_ver: str = subprocess.check_output([mafft, '--version'],
                                                 stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find mafft executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(mafft))
        msg = 'Please add mafft to your $PATH or supply a full mafft path in '\
              '{}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.debug('mafft found: {}'.format(mafft_ver.strip().decode()))

    # iqtree testing
    try:
        iqtree_ver: str = subprocess.check_output([iqtree, '--version'],
                                                  stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find iqtree2 executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(iqtree))
        msg = 'Please add iqtree2 to your $PATH or supply a full iqtree path '\
              'in {}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.debug('iqtree2 found: {}'.format(iqtree_ver.strip().decode()))

    exes: Dict[str, str] = {}
    exes['tblastn'] = tblastn
    exes['blastn'] = blastn
    exes['makeblastdb'] = makeblastdb
    exes['mafft'] = mafft
    exes['iqtree'] = iqtree
    return(exes)


def install_blast():
    logger = logging.getLogger(__name__)
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    scriptpath = os.path.dirname(os.path.abspath(filename))
    externalpath = os.path.join(os.path.dirname(scriptpath), 'external')
    if not os.path.exists(externalpath):
        os.mkdir(externalpath)
    version = '2.10.1'
    base_url = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{}'\
        .format(version)
    base_file = 'ncbi-blast-{}+-x64-'.format(version)
    if platform.system() == 'Linux':
        base_out_file = base_file + 'linux.tar.gz'
    elif platform.system() == 'Darwin':
        base_out_file = base_file + 'macosx.tar.gz'
    elif platform.system() == 'Windows':
        base_out_file = base_file + 'win64.tar.gz'
    out_file = os.path.join(externalpath, base_out_file)
    blast_tar_url = '/'.join([base_url, base_out_file])
    md5_url = '/'.join([base_url, base_out_file + '.md5'])
    logger.info('Downloading BLAST from {}'.format(blast_tar_url))
    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        with urllib.request.urlopen(blast_tar_url) as response, \
                open(out_file, 'wb') as blastfh:
            shutil.copyfileobj(response, blastfh)
    else:
        logging.info('BLAST already found.')
    blast_md5_digest = md5(p(out_file).read_bytes()).hexdigest()
    logging.info('calculated md5 digest: {}'.format(blast_md5_digest))
    with urllib.request.urlopen(md5_url) as response, \
            open(out_file + '.md5', 'wb') as md5fh:
        shutil.copyfileobj(response, md5fh)
    with open(out_file + '.md5', 'r') as md5fh:
        md5_hash = md5fh.readline().split(' ', 1)[0]
    logging.info('given md5 digest: {}'.format(md5_hash))
    if blast_md5_digest != md5_hash:
        logging.info('BLAST download was unsuccessful. Delete file and try '
                     'again.')
        exit(-1)
    logging.info('Unpacking BLAST tar.gz file.')
    shutil.unpack_archive(out_file, externalpath, 'gztar')
