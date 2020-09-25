#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import subprocess
import json
import shutil
import glob
from typing import List
from hashlib import blake2b
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

SUFFIXES = ['.nsq', '.nin', '.nhr', '.nto', '.not', '.ndb', '.ntf']


def remove_intermediates(runid: str, intermediates: List[str]) -> None:
    genome_tmpfiles = ['blast_results.tsv', 'keepsidx.json',
                       'single_copy.json', 'blast_filtered.tsv',
                       'expected_filt.json']
    if 'genome' in intermediates:
        for deldir in ['unaligned', 'aligned']:
            if os.path.exists(deldir):
                shutil.rmtree(deldir)
        for fullname in [os.path.join('.autoMLSA', fn) for
                         fn in genome_tmpfiles]:
            if os.path.exists(fullname):
                os.remove(fullname)
    for delfile in glob.iglob(runid + '.nex*'):
        os.remove(delfile)


def end_program(logger, code):
    """
    Program message including success or failure of the program.
    """
    if code == 1:
        msg = 'Program was stopped at an intermediate stage.'
        logger.info(msg)
    elif code == 0:
        msg = 'Program was a success! Congratulations!'
        logger.info(msg)
    else:
        msg = 'Program exiting with code ({}) indicating failure.'
        logger.info(msg.format(code))
        msg = 'Check error messages to resolve the problem.'
        logger.info(msg)
    exit(code)


def check_if_fasta(fa):
    """
    Checks first line of file for '>' to determine if FASTA, returns bool
    """
    for suffix in SUFFIXES:
        if suffix in fa:
            return False
    with open(fa, 'r') as fah:
        if fah.readline().startswith('>'):
            return True
        else:
            return False


def json_writer(fn, x):
    """
    Writes object to file as json format.
    """
    with open(fn, 'w') as fh:
        json.dump(x, fh, indent=4)
        fh.write('\n')


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
    exist = [os.path.exists(fasta + x) for x in SUFFIXES]
    if not all(exist):
        subprocess.run([makeblastdb, '-dbtype', 'nucl', '-in', fasta])


def checkpoint_reached(logger, stage):
    msg = 'Checkpoint reached {}. Stopping...'
    logger.info(msg.format(stage))
    end_program(logger, 1)


def checkpoint_tracker(fn_name):
    checkpath = os.path.join('.autoMLSA', 'checkpoint', fn_name)
    if not os.path.exists(checkpath):
        open(os.path.join(checkpath), 'w').close()


def main():
    sys.stderr.write('This program is not intended to be run on its own.\n')
    sys.stderr.write('Provides functionality to autoMLSA.py.\n')


if __name__ == '__main__':
    main()
