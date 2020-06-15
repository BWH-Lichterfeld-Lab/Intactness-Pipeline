#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform Hypermutation analysis based on Hypermut v2.0

@author: Ce Gao
"""

import re
import logging

from Bio import SeqIO
from scipy.stats import fisher_exact

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.hypermut')


def get_pos(pattern, string):
    """Get start position for each pattern match in string
    """
    return set(m.start() for m in re.finditer(pattern, string))


def _hypermut(ref_aln, qry_aln):
    """Perform pairwise hypermut analysis between reference and query seqs

    Parameters
    ----------
    ref_str: str
        reference seq
    qry_str: str
        query seq
    """
    # Too many local variables
    # pylint: disable=R0914

    aln = [(i, j) for i, j in zip(ref_aln, qry_aln) if i != '-' or j != '-']

    ref_str = ''.join(i for i, _ in aln)
    qry_str = ''.join(i for _, i in aln)

    indx_down_mut = get_pos('(?=.[AG][^C])', qry_str)
    indx_down_ctl = get_pos('(?=.([CT].|[AG]C))', qry_str)
    indx_from = get_pos('G', ref_str)
    indx_to = get_pos('A', qry_str)

    mut_pot = indx_from.intersection(indx_down_mut)
    ctl_pot = indx_from.intersection(indx_down_ctl)
    mut = mut_pot.intersection(indx_to)
    ctl = ctl_pot.intersection(indx_to)

    cnt_mut = len(mut)
    cnt_mut_pot = len(mut_pot)
    cnt_ctl = len(ctl)
    cnt_ctl_pot = len(ctl_pot)

    _, pval = fisher_exact([[cnt_mut, cnt_mut_pot - cnt_mut],
                            [cnt_ctl, cnt_ctl_pot - cnt_ctl]],
                           alternative='greater')
    try:
        ratio = (cnt_mut / cnt_mut_pot) / (cnt_ctl / cnt_ctl_pot)
    except ZeroDivisionError:
        ratio = float('Inf')

    psc = 'Yes' if pval < 0.05 else 'No'

    return (cnt_mut, cnt_mut_pot, cnt_ctl, cnt_ctl_pot, ratio, pval, psc)


def _hyper_long(ref, qry, seqs):
    """ hyper

    """
    # Perform hypermut analysis
    ret = _hypermut(str(ref.seq), str(qry.seq))

    seqs.call[qry.id]['hypermut'] = ret[-1]
    seqs.info[qry.id]['hypermut'] = ret[:-1]


def _hyper_drop(file_ref, file_drop, seqs):
    """Hypermutation for those short seqs
    """

    # Check the size of file_drop: stop early if no seqs get dropped
    import os.path

    try:
        f_size = os.path.getsize(file_drop)
    except OSError:
        return

    if f_size == 0:
        return

    # Run HypermutShort
    import tempfile
    fh_out = tempfile.NamedTemporaryFile(delete=False)
    from Bio.Blast.Applications import NcbiblastnCommandline
    cline = NcbiblastnCommandline(out=fh_out.name, query=file_drop,
                                  db=file_ref, outfmt=5,
                                  max_hsps=1, task='blastn')

    cline()
    from Bio.Blast import NCBIXML
    blast_records = NCBIXML.parse(open(fh_out.name))

    for blast_record in blast_records:
        qseqid = blast_record.query
        for aln in blast_record.alignments:
            hsps = aln.hsps[0]
            res = _hypermut(hsps.sbjct, hsps.query)

            seqs.call[qseqid]['hypermut'] = res[-1]
            seqs.info[qseqid]['hypermut'] = res[:-1]


def _write_hyper_summary(seqs, file_out):
    with open(file_out, 'wt') as f_o:
        header = "\t".join(['ID', 'Muts', 'PotenMuts', 'Ctls', 'PotenCtls',
                            'Ratio', 'p', 'Hypermut?'])
        print(header, file=f_o)

        line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"

        for qid in seqs.qids:
            print(line.format(qid, *seqs.info[qid]['hypermut'],
                              seqs.call[qid]['hypermut']),
                  file=f_o)


def hypermut(configs, seqs):
    """Hypermut analysis
    """
    # Logging
    logger.info('Running Hypermut analysis')

    file_aln = configs['file_aln']

    # File name of reference genome, e.g., ref.fasta
    file_ref = configs['file_ref']

    # File name of query seqs, e.g., contig.fasta
    file_seqs_drop = configs['file_seqs_drop']

    # Reference
    recs = SeqIO.parse(file_aln, 'fasta')
    for rec in recs:
        if 'HXB2' in rec.id:
            ref = rec
            break

    recs = SeqIO.parse(file_aln, 'fasta')
    for qry in recs:
        if 'HXB2' in qry.id:
            continue
        _hyper_long(ref, qry, seqs)

    # Run hypermut for long contigs
    for qid in seqs.qids:
        if seqs.call[qid]['is_hiv'] == 'No':
            seqs.call[qid]['hypermut'] = 'Pass'
            seqs.info[qid]['hypermut'] = tuple('Pass' for _ in range(6))

    # Run hypermut for short ones
    _hyper_drop(file_ref, file_seqs_drop, seqs)

    # Write report
    _write_hyper_summary(seqs, configs['file_out'])
