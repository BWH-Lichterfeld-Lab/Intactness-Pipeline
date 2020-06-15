#!/usr/bin/env python3
"""
Various alignment utilities

"""

import logging
import sys
from math import inf

from .utils import run_cmd

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.blast')


class BlastHit:
    """Blast Hit.

    Attributes
    ----------
    qlen: int
        query seq length
    strand: str
        subject seq strand
    sstart: int
        start position on subject
    min_sstart: int
        start boundary on subject
    max_send: int
        end boundary on subject
    aligned_pos: set of int
        positions on __query__ that can be aligned to sbjct
    """

    def __init__(self, qlen, min_sstart, max_send):
        self.qlen = qlen
        self.min_sstart = min_sstart
        self.max_send = max_send

        self.strand = 'NA'
        self.sstart = inf

        # dict, mapping sseqid to a _set_ of positions on __query__
        self.aligned_pos = set()

    def update(self, qstart, qend, sstart, send, sstrand):
        '''Update strand, sstart of a query

        If a hit is within [0, min_sstart] or [max_send, end of sbjct], it will
        not be included.

        If two hits have a overlap, the overlapping region will be counted only
        once.
        '''

        # pylint: disable=R0913
        # Too many arguments

        if self.strand == 'NA':
            # The strand has not been set yet
            self.strand = sstrand
        else:
            # This query has multiple hits of the same subject on both strands
            # Or we won't change
            if self.strand != sstrand:
                self.strand = 'mixture'

        if (sstart > self.max_send) or (send < self.min_sstart):
            return

        # if a blast has a overlap within 5'-end region [0, min_sstart],
        # only the region of [min_sstart, send] will be counted.
        # The start_shift is a inferred position in the query sequence by ref
        s_shift = 0 if sstart >= self.min_sstart else self.min_sstart - sstart

        # if a blast has a overlap in 3'-end region of [max_send, >max_ssend],
        # only the non-overlapped region, [sstart, max_send], will be counted
        e_shift = 0 if send <= self.max_send else send - self.max_send

        self.aligned_pos.update(range(qstart + s_shift, qend - e_shift))

        if sstart < self.sstart:
            self.sstart = sstart

    def report_hit(self):
        """Summary hit and report"""
        aln_len = len(self.aligned_pos)
        sstrand = self.strand
        qlen = self.qlen

        sstart = self.sstart
        sstart = sstart if sstart != inf else 'Pass'

        return aln_len, sstrand, qlen, sstart


def _run_blast(file_ldb, file_qry, max_eval, path_out):
    """Run blast

    Parameters:
        file_ldb: path to local database
        file_qry: path to query sequences file
        max_eval: maximum E-value level
        path_out: path to output file

    Output format. DONOT change, will break API
        0 qseqid means Query Seq-id
        1 qlen means Query length
        2 sseqid means Subject Seq-id
        3 slen means Subject length
        4 qstart means Start of alignment in query
        5 qend means End of alignment in query
        6 sstart means Start of alignment in subject
        7 send means End of alignment in subject
        8 evalue means Expect value
        9 bitscore means Bit score
        10 length means Alignment length
        11 pident means Percentage of identical matches
        12 nident means Number of identical matches
        13 gaps means Total number of gap
        14 sstrand means Subject strand
    """
    fmt = ("6 qseqid qlen sseqid slen qstart qend sstart send "
           "evalue bitscore length pident nident gaps sstrand")

    # Build blast program
    cmd = ['blastn',
           '-db', file_ldb,
           '-query', file_qry,
           '-outfmt', fmt,
           '-evalue', max_eval,
           '-out', path_out]

    run_cmd(cmd)


def _parse_blast_output(f_i, min_sstart, max_send):
    """Parse blast output
    """
    # A dict mapping qseqid's to BlastHit objects
    hits = {}

    # Get alignment results
    with open(f_i, 'U') as fh_i:
        for line in fh_i:
            tokens = line.strip().split('\t')

            # Query seq-id
            qseqid = tokens[0]
            qlen = int(tokens[1])
            qstart = int(tokens[4])
            qend = int(tokens[5])
            sstrand = tokens[-1]

            if sstrand == 'plus':
                sstart, send = int(tokens[6]), int(tokens[7])
            else:
                sstart, send = int(tokens[7]), int(tokens[6])

            if qseqid not in hits:
                hits[qseqid] = BlastHit(qlen, min_sstart, max_send)

            hits[qseqid].update(qstart, qend, sstart, send, sstrand)

    return hits


def _call_del_inv(seqs, hits, min_aln_len, min_aln_len_no_primer):
    """Call large deletion and inversion based on blast reulsts"""

    for qid in seqs.qids:
        if qid not in hits:
            seqs.call[qid]['is_hiv'] = 'No'
            seqs.call[qid]['deletion'] = 'Pass'
            seqs.call[qid]['inversion'] = 'Pass'
            seqs.info[qid]['blast'] = (0, seqs.qlen(qid), 'Pass', 'Pass')
        else:
            aln_length = len(hits[qid].aligned_pos)
            aln_strand = hits[qid].strand
            qlen = hits[qid].qlen
            sstart = hits[qid].sstart
            sstart = sstart if sstart != -1 else 'NA'

            seqs.info[qid]['blast'] = aln_length, qlen, sstart, aln_strand

            # Discard queries whose max alignment length are below a threshold
            if aln_length == 0:
                seqs.call[qid]['is_hiv'] = 'No'
                seqs.call[qid]['deletion'] = 'Pass'
                seqs.call[qid]['inversion'] = 'Pass'
            else:
                seqs.call[qid]['is_hiv'] = 'Yes'

                if aln_length >= min_aln_len:
                    seqs.call[qid]['deletion'] = 'No'
                elif seqs.call[qid]['primer'] == 'No' and \
                        aln_length > min_aln_len_no_primer:
                    seqs.call[qid]['deletion'] = 'No'
                else:
                    seqs.call[qid]['deletion'] = 'Yes'

                # Assign selected seq to +/- categories
                if aln_strand == 'plus':
                    seqs.call[qid]['inversion'] = 'No'
                elif aln_strand == 'minus':
                    seqs.call[qid]['inversion'] = 'No'
                # Discard queries mapped both strands, ie. aln[3] == 'mixture'
                else:
                    seqs.call[qid]['inversion'] = 'Yes'


def _write_blast_summary(seqs, path_sum):
    """Write a summary blast results
    """
    # Summarize results
    with open(path_sum, 'w') as f_o:
        header = ['Contig ID',
                  'Sample ID',
                  'Contig Length',
                  'Aligned Length',
                  'Aligned coverage of Contig',
                  'Ref Seq ID',
                  'Aligned Start at Ref',
                  'Ref Strand',
                  'Is HIV?',
                  'Primer?',
                  'Primer Seq'
                  'Large Deletion?',
                  'Internal Inversion?']

        print('\t'.join(header), file=f_o)

        sseqid = seqs.ref_id

        for qid in seqs.qids:
            sample_id = qid.split('_')[0]

            aln_len, qlen, sstart, strand = seqs.info[qid]['blast']
            is_hiv = seqs.call[qid]['is_hiv']
            deletion = seqs.call[qid]['deletion']
            inversion = seqs.call[qid]['inversion']

            line = '{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}'
            line = line.format(qid, sample_id, qlen, aln_len, aln_len/qlen,
                               sseqid, sstart, strand, is_hiv, deletion,
                               inversion)

            print(line, file=f_o)


def _prepare_qseqs(seqs, file_keep, file_keep_plus_ref, file_drop):
    """Prepare query sequences"""
    qids_keep_fwd = set()
    qids_keep_rev = set()
    qids_drop = set()
    for qid in seqs.qids:
        _, _, _, aln_strand = seqs.info[qid]['blast']
        if all([seqs.call[qid]['is_hiv'] == 'Yes',
                seqs.call[qid]['deletion'] == 'No',
                seqs.call[qid]['inversion'] == 'No']):
            if aln_strand == 'plus':
                qids_keep_fwd.add(qid)
            else:
                qids_keep_rev.add(qid)
        else:
            qids_drop.add(qid)

    msg = 'No. of seqs with forward mapping: {}'.format(len(qids_keep_fwd))
    logger.info(msg)
    msg = 'No. of seqs with reverse mapping: {}'.format(len(qids_keep_rev))
    logger.info(msg)
    msg = 'No. of seqs with Non-HIV/Del/Inv: {}'.format(len(qids_drop))
    logger.info(msg)

    if len(qids_keep_fwd) + len(qids_keep_rev) == 0:
        print('No sequence remained. Program was terminated!')
        sys.exit(1)

    # Get reverse sequence right
    seqs.take_rc_for(qids_keep_rev)

    qids_keep = qids_keep_fwd | qids_keep_rev
    qids_drop = set(seqs.qids) - qids_keep
    seqs.write(file_keep, qids_keep)
    seqs.write(file_keep_plus_ref, qids_keep, prepend_ref=True)
    seqs.write(file_drop, qids_drop)


def _prepare_del(seqs, file_del):
    qids_del = []
    qids_rev = []
    for qid in seqs.qids:
        aln_length, qlen, sstart, aln_strand = seqs.info[qid]['blast']
        if seqs.call[qid]['deletion'] == 'Yes' and \
                aln_length >= 7000 and sstart <= 2290:
            qids_del.append(qid)
            if aln_strand != 'plus':
                qids_rev.append(qid)

    seqs.take_rc_for(qids_rev)
    seqs.write(file_del, qids_del)


def blast(configs, seqs):
    """BLAST runner

    Conduct a blast search and parse the results.

    Parameters
    ----------
    configs: dict. Configurations, containing the following fields
        query: query sequence files
        database: local database folder
        output: blast output file
        max_evalue: blast output threshold
        min_alignment_start, max_alignment_end, min_alignment_length:
            alignment params
        region_of_gene: specific gene needs to be checked
    seqs: Sequence obj
        input
    """
    # Logging
    logger.info('Running Blast against local database')

    _run_blast(file_ldb=configs['file_ref'], file_qry=configs['file_qry'],
               max_eval=configs['max_eval'], path_out=configs['file_out'])

    # Parse the blast output
    hits = _parse_blast_output(configs['file_out'],
                               int(configs['min_alignment_start']),
                               int(configs['max_alignment_end']))

    # Call large deletion and internal inversion
    _call_del_inv(seqs, hits,
                  int(configs['min_alignment_length']),
                  int(configs['min_aln_len_no_primer']))

    _write_blast_summary(seqs, configs['file_summary'])

    _prepare_qseqs(seqs,
                   configs['file_seqs_keep'],
                   configs['file_seqs_keep_plus_ref'],
                   configs['file_seqs_drop'])

    _prepare_del(seqs, configs['file_seqs_del'])
