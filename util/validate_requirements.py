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

BLASTVER = '2.10.1'
MAFFTVER = '7.471'
IQTREEVER = '2.1.1'


def get_external_path() -> str:
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    scriptpath = os.path.dirname(os.path.abspath(filename))
    externalpath = os.path.join(os.path.dirname(scriptpath), 'external')
    if not os.path.exists(externalpath):
        os.mkdir(externalpath)
    return externalpath


def validate_requirements() -> Dict[str, str]:
    """Identifies installed software and sets paths for running

    input  - logger
    output - exes dict with paths to executables
    """
    logger = logging.getLogger(__name__)

    # Set local path to blast executable DIR
    BLASTPATH: str = ''
    externalpath = get_external_path()
    blastcheck = os.path.join(externalpath, 'ncbi-blast-{}+'.format(BLASTVER),
                              'bin')
    if os.path.exists(blastcheck):
        logger.debug('BLAST found in external dir.')
        BLASTPATH = blastcheck

    if platform.system() == 'Linux':
        mafftsuffix = 'linux64'
        iqtreesuffix = os.path.join('{}-Linux'.format(IQTREEVER), 'bin')
    elif platform.system() == 'Darwin':
        mafftsuffix = 'mac'
        iqtreesuffix = os.path.join('{}-MacOSX'.format(IQTREEVER), 'bin')
    elif platform.system() == 'Windows':
        mafftsuffix = ''
        iqtreesuffix = os.path.join('{}-Windows'.format(IQTREEVER), 'bin')

    mafftcheck = os.path.join(externalpath, 'mafft-{}'.format(mafftsuffix))
    iqtreecheck = os.path.join(externalpath, 'iqtree-{}'.format(iqtreesuffix))

    # Expects in $PATH, can hard-code to full path here
    tblastn: str = 'tblastn'
    blastn: str = 'blastn'
    makeblastdb: str = 'makeblastdb'
    mafft: str = 'mafft'
    if os.path.exists(mafftcheck):
        logger.debug('MAFFT found in external dir.')
        mafft = os.path.join(mafftcheck, 'mafft.bat')
    iqtree: str = 'iqtree2'
    if os.path.exists(iqtreecheck):
        logger.debug('iqtree found in external dir.')
        iqtree = os.path.join(iqtreecheck, 'iqtree2')

    if os.path.exists(os.path.join(BLASTPATH, tblastn)):
        tblastn = os.path.join(BLASTPATH, tblastn)
    if os.path.exists(os.path.join(BLASTPATH, blastn)):
        blastn = os.path.join(BLASTPATH, blastn)
    if os.path.exists(os.path.join(BLASTPATH, makeblastdb)):
        makeblastdb = os.path.join(BLASTPATH, makeblastdb)

    logger.debug('Checking tblastn: {}'.format(tblastn))
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

    logger.debug('Checking blastn: {}'.format(blastn))
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

    logger.debug('Checking makeblastdb: {}'.format(makeblastdb))
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
    logger.debug('Checking mafft: {}'.format(mafft))
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
    logger.debug('Checking iqtree: {}'.format(iqtree))
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


def download_file(out_file: str, url: str, program: str) -> None:
    logger = logging.getLogger(__name__)
    logger.info('Downloading {} from {}.'.format(program, url))
    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        with urllib.request.urlopen(url) as response, \
                open(out_file, 'wb') as tarfh:
            shutil.copyfileobj(response, tarfh)
    else:
        logger.info('{} already downloaded.'.format(program))


def init_logger(debug: bool = False) -> None:
    logger = logging.getLogger()
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # formatter = logging.Formatter('%(name)s - %(asctime)s - '
    #                               '%(levelname)s - %(message)s',
    #                               datefmt='%Y-%m-%d %H:%M:%S')
    formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(formatter)

    logger.addHandler(stderr_handler)


def install_blast() -> None:
    logger = logging.getLogger(__name__)
    externalpath = get_external_path()
    base_url = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{}'\
        .format(BLASTVER)
    base_file = 'ncbi-blast-{}+-x64-'.format(BLASTVER)
    if platform.system() == 'Linux':
        base_out_file = base_file + 'linux.tar.gz'
    elif platform.system() == 'Darwin':
        base_out_file = base_file + 'macosx.tar.gz'
    elif platform.system() == 'Windows':
        base_out_file = base_file + 'win64.tar.gz'
    out_file = os.path.join(externalpath, base_out_file)
    blast_tar_url = '/'.join([base_url, base_out_file])
    md5_url = '/'.join([base_url, base_out_file + '.md5'])
    download_file(out_file, blast_tar_url, 'BLAST')

    blast_md5_digest = md5(p(out_file).read_bytes()).hexdigest()
    logger.debug('calculated md5 digest: {}'.format(blast_md5_digest))

    download_file(out_file + '.md5', md5_url, 'BLAST md5')
    with open(out_file + '.md5', 'r') as md5fh:
        md5_hash = md5fh.readline().split(' ', 1)[0]
    logger.debug('given md5 digest: {}'.format(md5_hash))

    if blast_md5_digest != md5_hash:
        logger.info('BLAST download was unsuccessful. Delete file and try '
                    'again.')
        end_program(-1)
    tblastn = os.path.join(externalpath, 'ncbi-blast-{}+'.format(BLASTVER),
                           'bin', 'tblastn')
    if not os.path.exists(tblastn):
        logger.info('Unpacking BLAST tar.gz file.')
        shutil.unpack_archive(out_file, externalpath, 'gztar')
    logger.info('Checking tblastn version.')
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
        logger.info('tblastn found: {}'.format(tblastn_ver.strip().decode()))
        logger.info('Add {} to your $PATH to have access to this program '
                    .format(os.path.dirname(tblastn)) +
                    'outside of autoMLSA.')


def install_mafft() -> None:
    logger = logging.getLogger(__name__)

    base_url = 'https://mafft.cbrc.jp/alignment/software'
    base_file = 'mafft-{}'.format(MAFFTVER)
    if platform.system() == 'Linux':
        base_out_file = base_file + '-linux.tgz'
        md5_hash = '641cf15f90fc6634f9b2b2c7cf67e27d'
    elif platform.system() == 'Darwin':
        base_out_file = base_file + '-mac.zip'
        md5_hash = 'f57cd93ea546e47fbad975aecb58d303'
    elif platform.system() == 'Windows':
        logger.error('Windows is not currently supported for MAFFT install.')
        logger.error('Look at install options here: '
                     'https://mafft.cbrc.jp/alignment/software/windows.html')
        return

    externalpath = get_external_path()
    out_file = os.path.join(externalpath, base_out_file)
    tar_url = '/'.join([base_url, base_out_file])
    download_file(out_file, tar_url, 'MAFFT')

    mafft_md5_digest = md5(p(out_file).read_bytes()).hexdigest()
    logger.debug('calculated md5 digest: {}'.format(mafft_md5_digest))
    logger.debug('given md5 digest: {}'.format(md5_hash))
    if mafft_md5_digest != md5_hash:
        logger.info('MAFFT download was unsuccessful. Delete file and try '
                    'again.')
        exit(-1)
    if platform.system() == 'Linux':
        mafft = os.path.join(externalpath, 'mafft-linux64', 'mafft.bat')
        if not os.path.exists(mafft):
            logger.info('Unpacking MAFFT tar.gz file.')
            shutil.unpack_archive(out_file, externalpath, 'gztar')
    elif platform.system() == 'Darwin':
        mafft = os.path.join(externalpath, 'mafft-mac', 'mafft.bat')
        if not os.path.exists(mafft):
            logger.info('Unpacking MAFFT zip file.')
            shutil.unpack_archive(out_file, externalpath, 'zip')
    logger.info('Checking mafft version.')
    try:
        mafft_ver: str = subprocess.check_output([mafft, '--version'],
                                                 stderr=subprocess.STDOUT)
    except FileNotFoundError:
        msg = 'Unable to find mafft executable in path or in provided dir ' \
              '({})'
        logger.error(msg.format(mafft))
        msg = 'Please add mafft to your $PATH or supply a BLASTPATH in {}'
        logger.error(msg.format(__file__))
        end_program(78)
    else:
        logger.info('mafft found: {}'.format(mafft_ver.strip().decode()))
        logger.info('Add {} to your $PATH to have access to this program '
                    .format(os.path.dirname(mafft)) +
                    'ouside of autoMLSA.')


def install_iqtree() -> None:
    logger = logging.getLogger(__name__)

    base_url = 'https://github.com/iqtree/iqtree2/releases/download/v{}'\
        .format(IQTREEVER)
    base_file = 'iqtree-{}'.format(IQTREEVER)
    if platform.system() == 'Linux':
        base_out_file = base_file + '-Linux.tar.gz'
        md5_hash = '2e2dfeeb2a1e9e123f542b200f18340d'
    elif platform.system() == 'Darwin':
        base_out_file = base_file + '-MacOSX.zip'
        md5_hash = '45314fba262f200c32c37deceefbf5b1'
    elif platform.system() == 'Windows':
        base_out_file = base_file + '-Windows.zip'
        md5_hash = 'f1a25e89a5c31437125ae0c666d29ee8'

    externalpath = get_external_path()
    out_file = os.path.join(externalpath, base_out_file)
    tar_url = '/'.join([base_url, base_out_file])
    download_file(out_file, tar_url, 'iqtree')

    iqtree_md5_digest = md5(p(out_file).read_bytes()).hexdigest()
    logger.debug('calculated md5 digest: {}'.format(iqtree_md5_digest))
    logger.debug('given md5 digest: {}'.format(md5_hash))
    if iqtree_md5_digest != md5_hash:
        logger.info('iqtree2 download was unsuccessful. Delete file and try '
                    'again.')
        exit(-1)
    if platform.system() == 'Linux':
        iqtree = os.path.join(externalpath,
                              'iqtree-{}-Linux'.format(IQTREEVER),
                              'bin', 'iqtree2')
        if not os.path.exists(iqtree):
            logger.info('Unpacking iqtree tar.gz file.')
            shutil.unpack_archive(out_file, externalpath, 'gztar')
    elif platform.system() == 'Darwin':
        iqtree = os.path.join(externalpath,
                              'iqtree-{}-MacOSX'.format(IQTREEVER),
                              'bin', 'iqtree2')
        if not os.path.exists(iqtree):
            logger.info('Unpacking iqtree zip file.')
            shutil.unpack_archive(out_file, externalpath, 'zip')
    elif platform.system() == 'Windows':
        iqtree = os.path.join(externalpath,
                              'iqtree-{}-Windows'.format(IQTREEVER),
                              'bin', 'iqtree2.exe')
        if not os.path.exists(iqtree):
            logger.info('Unpacking iqtree zip file.')
            shutil.unpack_archive(out_file, externalpath, 'zip')
    logger.info('Checking iqtree version.')
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
        logger.info('iqtree2 found: {}'.format(iqtree_ver.strip().decode()))
        logger.info('Add {} to your $PATH to have access to this program '
                    .format(os.path.dirname(iqtree)) +
                    'outside of autoMLSA.')
