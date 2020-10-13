#!/usr/bin/env python3
from __future__ import print_function
import time
import os
import logging
import json
import subprocess
import glob
import shlex
from typing import List, Dict, Any
from util.parse_args import run_argparse
from util.validate_requirements import validate_requirements
from util.helper_functions import \
    (end_program, sanitize_path, generate_hash, checkpoint_reached,
     json_writer, checkpoint_tracker, remove_intermediates)
from util.blast_functions import \
    (read_blast_results, print_blast_summary, print_fasta_files,
     make_blast_database)
from util.configuration import read_config, validate_arguments, get_labels
from multiprocessing import Pool
from tqdm import tqdm
from collections import defaultdict
from shutil import copy
from Bio import SeqIO
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

__version__ = '2.2.0'

# CONFIGURATION SECTION #######################################################
# END CONFIGURATION SECTION ###################################################
# FORMATTING SECTION ##########################################################


def convert_fasta(rundir: str, fastas: List[str], labels: List[str],
                  makeblastdb: str, runid: str) -> List[str]:
    """Converts fasta file with placeholder label names

    input - unformatted fasta files
    return - Set of formatted fasta files
    also generates BLAST dbs per fasta
    """
    logger = logging.getLogger(__name__)
    fastadir: str = os.path.join(rundir, 'fasta')
    new_fastas: List[str] = []
    renamef: str = os.path.join(rundir, '.autoMLSA', 'rename.json')
    deleted: bool = False
    if os.path.exists(renamef):
        with open(renamef, 'r') as rf:
            rename_info: Dict[Any, Any] = json.load(rf)
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
            if check_hash(fasta, base, rename_info[base]['hash']):
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
    for labelfastaf in tqdm(new_fastas, desc='makeblastdb'):
        make_blast_database(makeblastdb, labelfastaf)
    json_writer(renamef, rename_info)

    expected_fastas_fn: str = os.path.join(rundir, '.autoMLSA',
                                           'expected_fastas.json')
    if os.path.exists(expected_fastas_fn):
        with open(expected_fastas_fn, 'r') as eff:
            expected_fastas: List[str] = json.load(eff)
    else:
        expected_fastas = []

    updated: bool = set(expected_fastas) != set(new_fastas)
    if updated or deleted:
        logger.debug('Found new genome files')
        remove_intermediates(runid, ['genome'])
    else:
        logger.debug('Found no new genome files')
    json_writer(expected_fastas_fn, new_fastas)

    return new_fastas


def check_hash(fasta: str, base: str, seqhash: str) -> bool:
    """
    Checks hash, deletes files related to genome if necessary
    return - value if something is deleted
    """
    logger = logging.getLogger(__name__)
    deleted: bool = False
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


def get_queries(runid: str, rundir: str, dups: bool, queries: List[str]) -> \
        List[str]:
    """Converts query fasta file(s) into individual fastas for BLAST

    input - FASTA files with queries as individual entries
    return - individual fasta files
    """
    logger = logging.getLogger(__name__)
    querydir: str = os.path.join(rundir, 'queries')
    backups: str = os.path.join(rundir, '.autoMLSA', 'backups')

    new_queries: List[str] = []
    seen: Dict[str, Dict[str, str]] = {}
    hashes: Dict[str, str] = {}
    if not os.path.exists(querydir):
        os.mkdir(querydir)
    for query_file in queries:
        i = 1
        query_base: str = os.path.basename(query_file)
        query_backup: str = os.path.join(backups, query_base)
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
                end_program(66)
        with open(query_file, 'r') as qf:
            for rec in SeqIO.parse(qf, 'fasta'):
                safeid: str = sanitize_path(rec.id)
                seqhash: str = generate_hash(str(rec.seq))

                if seqhash in hashes and hashes[seqhash] == safeid:
                    # Skip duplicate header + seq
                    continue
                elif seqhash in hashes and hashes[seqhash] != safeid:
                    # Same seq but different names, throw error
                    mismatchid: str = hashes[seqhash]
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
                                os.path.basename(seen[mismatchid]['path']),
                                seen[mismatchid]['unsafeid']))
                    msg = '+++++++++++++++++++++++++++++++++++++++++++++++'
                    logger.error(msg)
                    msg = 'Check your sequences to make sure they aren\'t '\
                        'misnamed, fix the problem, and try again.'
                    logger.error(msg)
                    end_program(65)
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
                    end_program(65)
                else:
                    if safeid in seen and dups:
                        msg = 'Keeping additional query {} with duplicate ID '\
                            'from file {}.'
                        logger.info(msg.format(safeid,
                                               os.path.basename(query_file)))
                    seen[safeid] = {}
                    seen[safeid]['path'] = query_file
                    seen[safeid]['index'] = str(i)
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
            expected_queries: List[str] = json.load(eqf)
    else:
        expected_queries = []

    updated: bool = set(expected_queries) != set(new_queries)
    if updated:
        if all(query in new_queries for query in expected_queries):
            logger.debug('Found new query sequences')
            remove_intermediates(runid, ['query'])
        else:
            logger.debug('Query sequences have been removed')
            remove_intermediates(runid, ['query', 'genome'])
    else:
        logger.debug('Found no new query sequences')
    json_writer(expected_queries_fn, new_queries)

    if not new_queries:
        msg = 'No query sequences found. Check your query file and try again.'
        logger.error(msg)
        end_program(-1)

    open(os.path.join('.autoMLSA', 'checkpoint', 'get_queries'), 'w').close()

    return new_queries

# END FORMATTING SECTION ######################################################
# BLAST SECTION ###############################################################


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
        if checkpoint:
            blastfile = os.path.join(rundir, 'blastcmds.txt')
            with open(blastfile, 'w') as fh:
                for cmd in cmds:
                    fh.write(' '.join([shlex.quote(x) for x in cmd]) + '\n')
            msg = 'BLAST commands written to {}. Exiting.'
            logger.info(msg.format(blastfile))
            checkpoint_reached('prior to BLAST searches')
        else:
            msg = 'Running {} BLAST searches using {} CPUs.'
            ncmds = len(cmds)
            logger.info(msg.format(ncmds, threads))
            p = Pool(threads)
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

# END BLAST SECTION ###########################################################
# ALIGNMENT SECTION ###########################################################


def run_mafft(threads: int, mafft: str, unaligned: List[str],
              checkpoint: bool) -> List[str]:
    """
    input  - unaligned fasta files per query
    return - list of aligned files per query
    """
    logger = logging.getLogger(__name__)
    base_cmd: List[str] = [mafft, '--localpair', '--maxiterate', str(1000),
                           '--thread', str(threads)]
    aligneddir: str = 'aligned'
    aligned: List[str] = []
    cmdstrs: List[str] = []
    logger.info('Aligning FASTA sequences, if necessary.')
    if not os.path.exists(aligneddir):
        os.mkdir(aligneddir)
    for unalign in tqdm(unaligned, desc='mafft'):
        outname: str = '{}/{}.aln'.format(
            aligneddir, os.path.basename(os.path.splitext(unalign)[0]))
        logname: str = '{}.log'.format(outname)
        logger.debug('Aligning and writing to {}'.format(outname))
        aligned.append(outname)
        if os.path.exists(outname):
            msg = 'Aligned file {} already found, skipping.'
            logger.debug(msg.format(outname))
        else:
            cmd: List[str] = base_cmd + [unalign]
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
        checkpoint_reached('prior to MAFFT alignments')
    elif checkpoint == 'postalign':
        checkpoint_reached('after MAFFT alignments')

    checkpoint_tracker('run_mafft')

    return aligned

# END ALIGNMENT SECTION #######################################################
# PHYLOGENY SECTION ###########################################################


def generate_nexus(runid: str, aligned: List[str], checkpoint: bool) -> str:
    """
    Generate nexus file containing partitions per aligned gene

    input  - list of aligned FASTA filenames
    return - path to nexus file
    """
    logger = logging.getLogger(__name__)
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

    nexus: str = '{}.nex'.format(runid)
    if os.path.exists(nexus):
        msg = 'Nexus file {} already found, skipping nexus file generation.'
        logger.info(msg.format(nexus))
    else:
        with open(nexus, 'w') as nex:
            nex.write('#nexus\n')
            nex.write('begin sets;\n')
            for fn in aligned:
                base = os.path.basename(os.path.splitext(fn)[0])
                nex.write('\tcharset {} = {}: *;\n'.format(base, fn))
            nex.write('end;\n')
    if checkpoint:
        checkpoint_reached('after generating nexus file {}'.format(nexus))

    checkpoint_tracker('generate_nexus')
    return(nexus)


def run_iqtree(threads: int, iqtree: str, nexus: str, outgroup: str) -> str:
    """
    Runs external iqtree command to generate phylogeny

    input  - nexus file
    return - path to output file
    """
    logger = logging.getLogger(__name__)
    out_tree: str = '{}.treefile'.format(nexus)
    logger.info('Generating phylogenetic tree ({}).'.format(out_tree))
    cmd: List[str] = [
        iqtree,
        '-p', nexus,
        '-B', '1000',
        '-alrt', '1000',
        '-m', 'MFP+MERGE',
        '-rcluster', '10',
        '-nt', str(threads)]
    if outgroup:
        cmd.extend(['-o', outgroup])
    if os.path.exists(out_tree):
        logger.info('Treefile {} already found. Skipping iqtree.'
                    .format(out_tree))
    else:
        logger.info('Running: {}'.format(' '.join(cmd)))
        subprocess.run(cmd)
        if not os.path.exists(out_tree):
            msg = 'iqtree2 seems to have failed.'
            logger.critical(msg)
            msg = 'Check the log files for error messages to see if they can '\
                'be resolved.'
            logger.critical(msg)
            end_program(73)

    checkpoint_tracker('run_iqtree')

    return out_tree


def exit_successfully(rundir: str, treefile: str) -> None:
    """Temporary command

    input  - rundir and treefile to print to log
    return - None
    """
    logger = logging.getLogger(__name__)

    msg = 'Your treefile is ready: {}/{}'
    logger.info(msg.format(rundir, treefile))
    end_program(0)

# END PHYLOGENY SECTION #######################################################


def main() -> None:
    # ARGPARSE SECTION
    args = run_argparse()
    logger = logging.getLogger(__name__)
    logger.info('Welcome to autoMLSA.py version {}'.format(__version__))
    time.sleep(1)

    # CONFIGURATION SECTION
    exes = validate_requirements()

    logger.info('Reconciling configuration settings.')
    config = {}
    if os.path.exists(args.configfile):
        config = read_config(args.configfile)
        if args.config is not None:
            msg = 'Config file specified {} is ignored as one is found in {}.'
            logger.warning(msg.format(args.config, args.rundir))
    elif args.config:
        config = read_config(args.config)
    args = validate_arguments(args, config, args.checkpoint == 'validate')

    # FORMATTING SECTION
    logger.info('Converting genome FASTA files for BLAST if necessary.')
    labels = get_labels(args.rundir, args.fasta)
    newfastas = convert_fasta(args.rundir, args.fasta, labels,
                              exes['makeblastdb'], args.runid)

    logger.info('Extracting query FASTA files if necessary.')
    queries = get_queries(args.runid, args.rundir, args.dups, args.query)

    # BLAST SECTION
    logger.info('Generating list of BLAST searches and outputs.')
    blastout = generate_blast_list(
        args.rundir, exes[args.program], queries, newfastas, args.evalue,
        args.threads, args.checkpoint == 'preblast')

    # BLAST output results, summary, and files
    blastres = read_blast_results(blastout, args.coverage, args.identity)
    blastfilt = print_blast_summary(
        args.runid, blastres, labels, args.allow_missing,
        args.missing_check, args.checkpoint == 'filtering')
    unaligned = print_fasta_files(blastfilt, labels)

    # ALIGNMENT SECTION
    aligned = run_mafft(
        args.threads, exes['mafft'], unaligned, args.checkpoint)

    # PHYLOGENY SECTION
    nexus = generate_nexus(args.runid, aligned, args.checkpoint == 'nexus')
    treefile = run_iqtree(args.threads, exes['iqtree'], nexus, args.outgroup)
    exit_successfully(args.rundir, treefile)


if __name__ == '__main__':
    main()
