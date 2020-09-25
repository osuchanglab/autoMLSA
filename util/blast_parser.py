#!/usr/bin/env python3
from __future__ import print_function
import sys
import numpy as np
import pandas as pd
import json
import os
# import csv
# import os.path
# import logging
from util.helper_functions import json_writer, checkpoint_reached,\
    checkpoint_tracker, remove_intermediates
from collections import defaultdict
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

SINGLE_COPY_ESTIMATE = 0.90


def read_blast_results(logger, blastfiles, coverage, identity):
    """
    Parses set of BLAST results. Expects list of files.

    input  - list of BLAST output filenames.
    return - pandas data frame of BLAST output
    """
    headers = ['qseqid', 'sseqid', 'saccver', 'pident', 'qlen', 'length',
               'bitscore', 'qcovhsp', 'stitle', 'sseq']
    # skip = 5
    logger.info('Reading BLAST results.')
    dtypes = {'qseqid': 'category',
              'sseqid': 'category',
              'saccver': 'string',
              'pident': np.float,
              'qlen': np.int,
              'length': np.int,
              'bitscore': np.float,
              'qcovhsp': np.int,
              'stitle': 'string',
              'sseq': 'string'}
    results_fn = os.path.join('.autoMLSA', 'blast_results.tsv')
    if os.path.exists(results_fn):
        results = pd.read_csv(results_fn, dtype=dtypes, sep='\t')
    else:
        results = pd.DataFrame(columns=headers)
        for blastfile in blastfiles:
            data = pd.read_csv(blastfile, header=None, names=headers,
                               comment='#', dtype=dtypes, sep='\t')
            results = results.append(data, ignore_index=True)
        results.query('(pident >= @identity) & (qcovhsp >= @coverage)',
                      inplace=True)
        results.to_csv(results_fn, sep='\t')
    checkpoint_tracker('read_blast_results')
    return(results)


def print_blast_summary(logger, runid, blastout, labels, nallowed,
                        missing_check, checkpoint):
    """
    Generates summary of BLAST results.

    input  - pandas data frame with BLAST results
    return - pandas data frame with BLAST results to keep
    """
    logger.info('Summarizing and filtering BLAST hits.')
    summary = defaultdict(dict)
    missing_count = defaultdict(list)
    genome_missing_dict = {}
    keeps = []
    keepsidx = []
    genome_map = {k: i for i, k in enumerate(labels)}

    grouped = blastout.groupby(['sseqid', 'qseqid'])
    presence_matrix = grouped.size().unstack('qseqid', fill_value=0)
    label_indexes = list(map(int, presence_matrix.index.values))
    presence_matrix.index = [labels[i] for i in label_indexes]
    indexes = presence_matrix.index.tolist()
    presence_matrix.to_csv('presence_matrix.tsv', sep='\t')
    summary['queries']['names'] = presence_matrix.columns.tolist()
    summary['queries']['count'] = len(presence_matrix.columns)
    summary['genomes']['names'] = presence_matrix.index.tolist()
    summary['genomes']['indexes'] = indexes
    summary['genomes']['count'] = len(presence_matrix.index)

    # Get zero counts for filtering
    # Queries
    summary['queries']['missing'] = defaultdict(dict)
    zero_counts_query = presence_matrix.apply(
        lambda x: x[x == 0].index.tolist(), axis=0)
    zero_counts_query_dict = zero_counts_query.to_dict()
    for k, v in zero_counts_query_dict.items():
        summary['queries']['missing'][k]['genomes'] = v
        summary['queries']['missing'][k]['count'] = len(v)
        summary['queries']['missing'][k]['percent'] = \
            '{:.2f}'.format(len(v) / summary['genomes']['count'] * 100)
    summary['queries']['missing'] = dict(summary['queries']['missing'])

    # Genomes
    summary['genomes']['missing'] = defaultdict(dict)
    zero_counts_genome = presence_matrix.apply(
        lambda x: x[x == 0].index.tolist(), axis=1)
    zero_counts_genome_dict = zero_counts_genome.to_dict()
    for k, v in zero_counts_genome_dict.items():
        summary['genomes']['missing'][k]['queries'] = v
        summary['genomes']['missing'][k]['count'] = len(v)
        summary['genomes']['missing'][k]['percent'] = \
            '{:.2f}'.format(len(v) / summary['queries']['count'] * 100)
        missing_count[len(v)].append(k)
        if v != []:
            genome_missing_dict[k] = v
        else:
            keeps.append(k)
    summary['genomes']['missing'] = dict(summary['genomes']['missing'])

    for k in genome_missing_dict:
        nmissing = len(genome_missing_dict[k])
        if nmissing > nallowed:
            msg = 'Genome {} is going to be removed due to missing queries.'
            logger.warning(msg.format(k))
            msg = 'Increase --allow_missing to {} from {} to keep this genome.'
            logger.warning(msg.format(nmissing, nallowed))
        else:
            msg = 'Keeping genome {} missing {} queries due to --allow_missing'
            logger.debug(msg.format(k, nmissing))
            keeps.append(k)
    msg = 'Keeping these genomes ({}):\n         - {}'
    logger.info(msg.format(len(keeps), '\n         - '.join(keeps)))
    keepsidx = [str(genome_map[x]) for x in keeps]
    json_writer(os.path.join('.autoMLSA', 'keepsidx.json'), keepsidx)
    blastout_keeps = blastout.query('sseqid in @keepsidx')\
        .sort_values(['bitscore', 'qcovhsp'], ascending=False)\
        .drop_duplicates(['qseqid', 'sseqid']).sort_index()

    blastout_keeps.to_csv(os.path.join('.autoMLSA', 'blast_filtered.tsv'),
                          sep='\t')

    # Estimate single copy
    # Currently UNUSED #
    single_copy = presence_matrix.apply(
        lambda x: ((x.values == 1).sum() / len(presence_matrix)) >
        SINGLE_COPY_ESTIMATE,
        axis=0)
    json_writer(os.path.join('.autoMLSA', 'single_copy.json'),
                single_copy.to_dict())
    json_writer('missing_counts.json', missing_count)
    json_writer('blast_summary.json', summary)
    json_writer('missing_by_genome.json', genome_missing_dict)
    if len(keeps) < summary['genomes']['count']:
        logger.warning('Some genomes ({}) are missing one or more genes'
                       .format(len(genome_missing_dict)))
        logger.warning('Check out the presence_matrix.tsv, and '
                       'missing_by_genome.json files to review.')
        if not missing_check:
            logger.warning('Use --missing_check to continue.')
            exit(1)
    expected_filt_fn = os.path.join('.autoMLSA', 'expected_filt.json')
    if os.path.exists(expected_filt_fn):
        with open(expected_filt_fn, 'r') as eqf:
            expected_filt = json.load(eqf)
    else:
        expected_filt = []

    updated = set(expected_filt) != set(keeps)
    if updated:
        logger.debug('Found new filtered sequences')
        logger.debug('Removing downstream files, if present')
        remove_intermediates(runid, ['genome'])
    if not updated:
        logger.debug('Found no new filtered sequences')
    json_writer(expected_filt_fn, keeps)

    if checkpoint:
        checkpoint_reached(logger, 'after BLAST result filtering')

    checkpoint_tracker('print_blast_summary')

    return(blastout_keeps)


def print_fasta_files(logger, blastout, labels):
    fastdir = 'unaligned'
    labels = [os.path.splitext(x)[0] for x in labels]
    unaligned = []
    if not os.path.exists(fastdir):
        os.mkdir(fastdir)
    msg = 'Writing unaligned FASTA sequences.'
    logger.info(msg)
    for name, group in blastout.groupby('qseqid'):
        fasta = os.path.join(fastdir, '{}.fas'.format(name))
        unaligned.append(fasta)
        if os.path.exists(fasta):
            logger.debug('Unaligned file {} already found, skipping.'
                         .format(fasta))
        else:
            if os.path.exists(fasta):
                text = 'Overwriting'
            else:
                text = 'Writing'
            logger.debug('{} {} FASTA file'.format(text, fasta))
            with open(fasta, 'w') as fh:
                for row in group.itertuples():
                    fh.write('>{}\n'.format(labels[int(row.sseqid)]))
                    fh.write('{}\n'.format(row.sseq.replace('-', '')))

    checkpoint_tracker('print_fasta_files')

    return(unaligned)


def main():
    sys.stderr.write('This program is not intended to be run on its own.\n')
    sys.stderr.write('Provides functionality to autoMLSA.py.\n')


if __name__ == '__main__':
    main()
