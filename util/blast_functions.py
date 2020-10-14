#!/usr/bin/env python3
from __future__ import print_function
import subprocess
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import json
import os
import logging
import multiprocessing
import shlex
from mutltiprocessing import Pool
from typing import List, Dict, Any, DefaultDict
from .helper_functions import \
    (json_writer, checkpoint_reached, checkpoint_tracker, remove_intermediates,
     end_program, SUFFIXES)
from collections import defaultdict
from tqdm import tqdm  # type: ignore
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

SINGLE_COPY_ESTIMATE = 0.90


def make_blast_database(makeblastdb: str, fasta: str) -> None:
    """Generates blast DB if necessary per fasta

    input  - fasta file
    return - null
    """
    exist: List[bool] = [os.path.exists(fasta + x) for x in SUFFIXES]
    if not all(exist):
        with open(os.path.join('fasta', 'makeblastdb.log'), 'a') as textlog:
            subprocess.run([makeblastdb, '-dbtype', 'nucl', '-in', fasta],
                           stdout=textlog, stderr=subprocess.STDOUT)


def read_blast_results(blastfiles: List[str], coverage: int, identity: int) ->\
        pd.DataFrame:
    """
    Parses set of BLAST results. Expects list of files.

    input  - list of BLAST output filenames.
    return - pandas data frame of BLAST output
    """
    logger: logging.Logger = logging.getLogger(__name__)
    headers: List[str] = ['qseqid', 'sseqid', 'saccver', 'pident', 'qlen',
                          'length', 'bitscore', 'qcovhsp', 'stitle', 'sseq']

    logger.info('Reading BLAST results.')
    dtypes: Dict[str, Any] = {'qseqid': 'category',
                              'sseqid': 'category',
                              'saccver': 'string',
                              'pident': np.float,
                              'qlen': np.int,
                              'length': np.int,
                              'bitscore': np.float,
                              'qcovhsp': np.int,
                              'stitle': 'string',
                              'sseq': 'string'}

    results_fn: str = os.path.join('.autoMLSA', 'blast_results.tsv')
    if os.path.exists(results_fn):
        logger.debug('Reading from existing BLAST results.')
        results: pd.DataFrame = pd.read_csv(results_fn, dtype=dtypes, sep='\t')
    else:
        results = pd.DataFrame(columns=headers)
        for blastfile in tqdm(blastfiles, 'Reading BLAST Result'):
            data: pd.DataFrame = pd.read_csv(blastfile, header=None,
                                             names=headers, comment='#',
                                             dtype=dtypes, sep='\t')
            results = results.append(data, ignore_index=True)
        results.query('(pident >= @identity) & (qcovhsp >= @coverage)',
                      inplace=True)
        results.to_csv(results_fn, sep='\t')
    checkpoint_tracker('read_blast_results')
    return(results)


def print_blast_summary(runid: str, blastout: pd.DataFrame, labels: List[str],
                        nallowed: int, missing_check: bool,
                        checkpoint: bool) -> pd.DataFrame:
    """
    Generates summary of BLAST results.

    input  - pandas data frame with BLAST results
    return - pandas data frame with BLAST results to keep
    """
    logger: logging.Logger = logging.getLogger(__name__)
    logger.info('Summarizing and filtering BLAST hits.')
    summary: DefaultDict[str, Dict[str, Any]] = defaultdict(dict)
    missing_count: DefaultDict[int, List[str]] = defaultdict(list)
    genome_missing_dict: Dict = {}
    keeps: List[str] = []
    keepsidx: List[str] = []
    genome_map: Dict[str, int] = {k: i for i, k in enumerate(labels)}

    grouped: pd.DataFrame = blastout.groupby(['sseqid', 'qseqid'])
    presence_matrix: pd.DataFrame = grouped.size().unstack('qseqid',
                                                           fill_value=0)
    label_indexes = list(map(int, presence_matrix.index.values))
    indexes: List[int] = presence_matrix.index.tolist()
    presence_matrix.index = [labels[i] for i in label_indexes]
    presence_matrix.to_csv('presence_matrix.tsv', sep='\t')
    summary['queries']['names'] = presence_matrix.columns.tolist()
    summary['queries']['count'] = len(presence_matrix.columns)
    summary['genomes']['names'] = \
        presence_matrix.index.tolist()  # type: ignore
    summary['genomes']['indexes'] = indexes
    summary['genomes']['count'] = len(presence_matrix.index)

    # Get zero counts for filtering
    # Queries
    summary['queries']['missing'] = defaultdict(dict)
    zero_counts_query: pd.DataFrame = presence_matrix.apply(
        lambda x: x[x == 0].index.tolist(), axis=0)
    zero_counts_query_dict: Dict[str, List[str]] = zero_counts_query.to_dict()
    for k, v in zero_counts_query_dict.items():
        summary['queries']['missing'][k]['genomes'] = v
        summary['queries']['missing'][k]['count'] = len(v)
        summary['queries']['missing'][k]['percent'] = \
            '{:.2f}'.format(len(v) / summary['genomes']['count'] * 100)
    summary['queries']['missing'] = dict(summary['queries']['missing'])

    # Genomes
    summary['genomes']['missing'] = defaultdict(dict)
    zero_counts_genome: pd.DataFrame = presence_matrix.apply(
        lambda x: x[x == 0].index.tolist(), axis=1)
    zero_counts_genome_dict: Dict[str, List[str]] = \
        zero_counts_genome.to_dict()
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
        nmissing: int = len(genome_missing_dict[k])
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
    blastout_keeps: pd.DataFrame = blastout.query('sseqid in @keepsidx')\
        .sort_values(['bitscore', 'qcovhsp'], ascending=False)\
        .drop_duplicates(['qseqid', 'sseqid']).sort_index()

    blastout_keeps.to_csv(os.path.join('.autoMLSA', 'blast_filtered.tsv'),
                          sep='\t')

    # Estimate single copy
    # Currently UNUSED #
    single_copy: pd.DefaultDict = presence_matrix.apply(
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
            end_program(1)
    expected_filt_fn: str = os.path.join('.autoMLSA', 'expected_filt.json')
    if os.path.exists(expected_filt_fn):
        with open(expected_filt_fn, 'r') as eqf:
            expected_filt: List[str] = json.load(eqf)
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
        checkpoint_reached('after BLAST result filtering')

    checkpoint_tracker('print_blast_summary')

    return(blastout_keeps)


def print_fasta_files(blastout: pd.DataFrame, labels: List[str]) -> List[str]:
    """
    Writes unaligned FASTA files as output from BLAST results

    input  - pandas DataFrame with blast output, list of labels
    return - List of unaligned FASTA files
    """
    logger: logging.Logger = logging.getLogger(__name__)
    fastdir: str = 'unaligned'
    labels = [os.path.splitext(x)[0] for x in labels]
    unaligned: List[str] = []
    if not os.path.exists(fastdir):
        os.mkdir(fastdir)
    msg = 'Writing unaligned FASTA sequences, if necessary.'
    logger.info(msg)
    for name, group in blastout.groupby('qseqid'):
        fasta: str = os.path.join(fastdir, '{}.fas'.format(name))
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


def generate_blast_list(rundir: str, exe: str, queries: List[str], targets:
                        List[str], evalue: float, threads: int,
                        checkpoint: bool) -> List[str]:
    """
    Generates list of BLAST searches that need to be run

    input  - blast type, exes, queries, and targets
    output - list of blast output names
    """
    logger = logging.getLogger(__name__)
    cmds: List[List[str]] = []
    # blastout = defaultdict(list)
    blastout: List[str] = []
    blastdir: str = os.path.join(rundir, 'blast')
    outfmt: str = '7 qseqid sseqid saccver pident qlen length bitscore '\
        'qcovhsp stitle sseq'
    base_cmd: List[str] = [exe, '-evalue', str(evalue), '-outfmt', outfmt]
    cmd: List[str] = []
    if not os.path.exists(blastdir):
        os.mkdir(blastdir)
    for target in targets:
        target_base: str = os.path.basename(target)
        for query in queries:
            query_base: str = os.path.splitext(os.path.basename(query))[0]
            outname: str = '{}_vs_{}.tab'.format(query_base, target_base)
            outpath: str = os.path.join(blastdir, outname)
            blastout.append(os.path.join(outpath))
            cmd = base_cmd + ['-db', target,
                              '-query', query,
                              '-out', outpath]
            if not os.path.exists(outpath) or os.path.getsize(outpath) == 0:
                cmds.append(cmd)

    if cmds:
        if checkpoint:
            blastfile: str = os.path.join(rundir, 'blastcmds.txt')
            with open(blastfile, 'w') as fh:
                for cmd in cmds:
                    fh.write(' '.join([shlex.quote(x) for x in cmd]) + '\n')
            msg = 'BLAST commands written to {}. Exiting.'
            logger.info(msg.format(blastfile))
            checkpoint_reached('prior to BLAST searches')
        else:
            msg = 'Running {} BLAST searches using {} CPUs.'
            ncmds: int = len(cmds)
            logger.info(msg.format(ncmds, threads))
            p: multiprocessing.pool.Pool = Pool(threads)
            with p:
                with tqdm(total=ncmds, desc='BLAST Search') as pbar:
                    for i, _ in enumerate(
                            p.imap_unordered(subprocess.run, cmds)):
                        pbar.update()
                # p.map(subprocess.run, cmds)
    else:
        logger.info('No BLAST searches remaining. Moving to parse.')

    open(os.path.join('.autoMLSA', 'checkpoint', 'generate_blast_list'),
         'w').close()
    return blastout