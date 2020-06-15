#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 00:09:49 2017

@author: Rong Chen
@author: Ce Gao

"""

import logging
from collections import Counter

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.summary')

# pylint: disable=R0902

# Disable too many return
# pylint: disable=R0911


def _call(calls):
    """Make final call"""

    final_call = ''
    if calls['is_hiv'] == 'No':
        final_call = 'NonHIV'
        return final_call

    if calls['deletion'] == 'Yes':
        final_call = 'Large Deletion'
        if calls['inversion'] == 'Yes':
            final_call += ' with Internal Inversion'
        elif calls['hypermut'] == 'Yes':
            final_call += ' with Hypermut'
        return final_call

    if calls['inversion'] == 'Yes':
        final_call = 'Internal Inversion'
        return final_call

    if calls['hypermut'] == 'Yes':
        final_call = 'Hypermut'
        return final_call

    if calls['psc'] == 'Yes':
        final_call = 'Premature Stop Codon'
        return final_call

    if calls['defect'] == 'Yes' and calls['primer'] == 'Yes':
        final_call = "5' defect"
        return final_call

    if calls['primer'] == 'No':
        final_call = 'Inferred Intact'
        return final_call

    return 'Intact'


def summary(configs, seqs):
    """Final call"""

    # Disable too many local variables
    # pylint: disable=R0914

    # Logging
    logger.info('Generating summaries')

    # Multiple contig sample
    sample_ids = [qid.split('_')[0] for qid in seqs.qids]
    freq = Counter(sample_ids)
    sids_multi_seq = set(sid for sid, cnt in freq.items() if cnt > 1)

    # Multiple contig sample
    sample_ids = [qid.split('_')[0] for qid in seqs.qids
                  if seqs.call[qid]['is_hiv'] == 'Yes']
    freq = Counter(sample_ids)
    sids_multi_hiv = set(sid for sid, cnt in freq.items() if cnt > 1)

    with open(configs['file_out'], 'w') as fh_o:
        cols = ["Contig ID",
                "Sample ID",
                "Multi-Contig Sample?",
                "Multi-HIV Sample?",
                "Contig Length",
                "Aligned Length",
                "Aligned coverage of Contig",
                "Ref Seq ID",
                "Aligned Start at Ref",
                "Ref Strand",
                "Is HIV?",
                "Primer",
                "Primer Seq",
                "Large Deletion?",
                "Internal Inversion?",
                "Hypermut?",
                "Hypermut pval",
                "PSC?",
                # "PSC Type",
                "gag",
                "pol",
                "env",
                "5' Defect",
                "5' Gaps",
                "5' Inserts",
                "Gag Start Codon Missing?",
                "Gag Start Seq",
                "Final Call",
                "Comments",
                "Contig Sequence"]

        fmts = ["{}",       # Contig ID
                "{}",       # Sample ID
                "{}",       # Multi-Contig Sample?
                "{}",       # Multi-HIV Sample?
                "{}",       # Contig Length
                "{}",       # Aligned Length
                "{:.2f}%",  # Aligned coverage of Contig
                "{}",       # Ref Seq ID
                "{}",       # Aligned Start at Ref
                "{}",       # Ref Strand
                "{}",       # Is HIV?
                "{}",       # Primer
                "{}",       # Primer Seq
                "{}",       # Large Deletion?
                "{}",       # Internal Inversion?
                "{}",       # Hypermut?
                "{}",       # Hypermut pval
                "{}",       # PSC?
                # "{}",       # PSC Type
                "{}",       # gag
                "{}",       # pol
                "{}",       # env
                "{}",       # 5' Defect
                "{}",       # 5' Gaps
                "{}",       # 5' Inserts
                "{}",       # Gag Start Codon
                "{}",       # Gag Start Seq
                "{}",       # Final Call
                "{}",       # Comments
                "{}"]       # Contig sequence

        header = ','.join(cols)
        print(header, file=fh_o)

        fmt_str = ','.join(fmts)
        ref_id = seqs.ref_id
        for qid in seqs.qids:
            sample_id = qid.split('_')[0]
            multi_seq_sample = 'Yes' if sample_id in sids_multi_seq else 'No'
            multi_hiv_sample = 'Yes' if sample_id in sids_multi_hiv else 'No'

            call = seqs.call[qid]
            info = seqs.info[qid]

            aln_len, qlen, sstart, strand = info['blast']
            is_hiv = call['is_hiv']
            primer = call['primer']
            primer_seq = ';'.join(info['primer'])
            deletion = call['deletion']
            inversion = call['inversion']
            hypermut = call['hypermut']
            hypermut_p = info['hypermut'][-1]
            psc = call['psc']
            psc_info = info['psc']
            gag = psc_info['Gag']
            pol = psc_info['Pol']
            env = psc_info['Env']
            defect = call['defect']
            defect_gaps, defect_inserts = info['defect']
            gag_codon = call['gag_codon']
            gag_codon_seq = info['gag_codon']
            comment = seqs.comments[qid]
            qseq = str(seqs.qry[qid].seq)

            line = fmt_str.format(qid,
                                  sample_id,
                                  multi_seq_sample,
                                  multi_hiv_sample,
                                  qlen,
                                  aln_len,
                                  aln_len/qlen*100,
                                  ref_id,
                                  sstart,
                                  strand,
                                  is_hiv,
                                  primer,
                                  primer_seq,
                                  deletion,
                                  inversion,
                                  hypermut,
                                  hypermut_p,
                                  psc,
                                  # psc_type,
                                  gag,
                                  pol,
                                  env,
                                  defect,
                                  defect_gaps,
                                  defect_inserts,
                                  gag_codon,
                                  gag_codon_seq,
                                  _call(call),
                                  comment,
                                  qseq)
            print(line, file=fh_o)
