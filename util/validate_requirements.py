#!/usr/bin/env python3
from __future__ import print_function
import sys
import os.path
import subprocess
from util.helper_functions import end_program
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)


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
              'in {}'
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


def main():
    sys.stderr.write('This program is not intended to be run on its own.\n')
    sys.stderr.write('Provides functionality to autoMLSA.py.\n')


if __name__ == '__main__':
    main()
