#!/usr/bin/env python3
from __future__ import print_function
import time
import os
import logging
import argparse
import pandas as pd  # type: ignore
from typing import List, Dict, Any
from util.parse_args import run_argparse
from util.validate_requirements import validate_requirements
from util.helper_functions import exit_successfully
from util.blast_functions import \
    (read_blast_results, print_blast_summary, print_fasta_files,
     generate_blast_list)
from util.configuration import read_config, validate_arguments, get_labels
from util.formatting import convert_fasta, get_queries
from util.mafft import run_mafft
from util.phylogeny import generate_nexus, run_iqtree
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

__version__ = '2.2.0'


def main() -> None:
    # ARGPARSE SECTION --------------------------------------------------------
    args: argparse.Namespace = run_argparse()
    logger = logging.getLogger(__name__)
    logger.info('Welcome to autoMLSA.py version {}'.format(__version__))
    time.sleep(1)

    # CONFIGURATION SECTION ---------------------------------------------------
    exes: Dict[str, str] = validate_requirements()

    logger.info('Reconciling configuration settings.')
    config: Dict[str, Any] = {}
    if os.path.exists(args.configfile):
        config = read_config(args.configfile)
        if args.config is not None:
            msg = 'Config file specified {} is ignored as one is found in {}.'
            logger.warning(msg.format(args.config, args.rundir))
    elif args.config:
        config = read_config(args.config)
    args = validate_arguments(args, config, args.checkpoint == 'validate')

    # FORMATTING SECTION ------------------------------------------------------
    logger.info('Converting genome FASTA files for BLAST if necessary.')
    labels: List[str] = get_labels(args.rundir, args.fasta)
    newfastas: List[str] = convert_fasta(args.rundir, args.fasta, labels,
                                         exes['makeblastdb'], args.runid)

    logger.info('Extracting query FASTA files if necessary.')
    queries: List[str] = get_queries(args.runid, args.rundir, args.dups,
                                     args.query)

    # BLAST SECTION -----------------------------------------------------------
    logger.info('Generating list of BLAST searches and outputs.')
    blastout: List[str] = generate_blast_list(
        args.rundir, exes[args.program], queries, newfastas, args.evalue,
        args.threads, args.checkpoint == 'preblast')

    # BLAST output results, summary, and files --------------------------------
    blastres: pd.DataFrame = read_blast_results(blastout, args.coverage,
                                                args.identity)
    blastfilt: pd.DataFrame = print_blast_summary(
        args.runid, blastres, labels, args.allow_missing,
        args.missing_check, args.checkpoint == 'filtering')
    unaligned: List[str] = print_fasta_files(blastfilt, labels)

    # ALIGNMENT SECTION -------------------------------------------------------
    aligned: List[str] = run_mafft(
        args.threads, exes['mafft'], unaligned, args.checkpoint)

    # PHYLOGENY SECTION -------------------------------------------------------
    nexus: str = generate_nexus(args.runid, aligned, args.checkpoint ==
                                'nexus')
    treefile: str = run_iqtree(args.threads, exes['iqtree'], nexus,
                               args.outgroup)
    exit_successfully(args.rundir, treefile)


if __name__ == '__main__':
    main()
