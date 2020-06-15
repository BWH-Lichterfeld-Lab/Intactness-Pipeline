"""Premature stop codon, frameshift"""

import logging
import os.path
from math import floor
import re

from Bio import SeqIO

from .GeneCutter import submit_GC, process_GC

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.psc')

def psc(configs, seqs):
    """Premature stop codon"""
    # Logging
    logger.info('Checking PSC')

    # Gene Cutter
    psc_list = []
    for qid in seqs.qids:
        call = seqs.call[qid]
        if call['is_hiv'] == 'Yes' and call['deletion'] == 'No' and \
                call['inversion'] == 'No' and call['hypermut'] == 'No':
            psc_list.append(qid)

    if len(psc_list) == 0:
        for qid in seqs.qids:
            seqs.call[qid]['psc'] = 'Pass'
            seqs.info[qid]['psc'] = {'Gag': 'Pass', 'Pol': 'Pass',
                                     'Env': 'Pass'}
        return

    seqs.write('data/seqs/seqs_psc.fasta', psc_list)
    submit_GC(configs['email'])
    process_GC()

    for gene in ['Gag', 'Pol', 'Env']:
        file_aln = configs['path_out'] + \
            '/Gene_Cutter/{}.aa.fasta'.format(gene)

        try:
            f_size = os.path.getsize(file_aln)
        except OSError:
            continue

        if f_size == 0:
            continue

        recs = SeqIO.parse(file_aln, 'fasta')

        # Reference
        ref = next(recs)
        len_ref = len(ref.seq.ungap('-')) - 1

        pos_start = int(configs[gene])

        for qry in recs:
            qry.id = qry.id.rstrip('_' + gene)
            cur_pos = 0  # Current ungapped position on ref
            num_ins = 0  # insertion
            num_del = 0  # deletion
            num_idn = 0  # identical
            num_mis = 0  # mismatch
            num_gap = 0  # gap in the last protein

            for c_ref, c_qry in zip(ref.upper(), qry.upper()):
                if c_ref == '-' and c_qry == '-':
                    continue
                elif c_ref == '-' and c_qry != '-':
                    num_ins += 1
                elif c_ref != '-' and c_qry == '-':
                    num_del += 1
                    cur_pos += 1
                    if cur_pos >= pos_start:
                        num_gap += 1
                elif c_ref == c_qry:
                    num_idn += 1
                    cur_pos += 1
                else:
                    num_mis += 1

            if num_gap >= floor(len_ref * 0.05):
                seqs.call[qry.id]['psc'] = 'Yes'
                msg = "Too many gaps in the last protein of {};".format(gene)
                seqs.comments[qry.id] += msg

            len_qry = len(qry.seq.ungap('-'))
            info = seqs.info[qry.id].setdefault('psc', {})
            fmt_str = "len({});iden({});mis({});ins({});del({})"

            info[gene] = fmt_str.format(len_qry, num_idn, num_mis, num_ins,
                                        num_del)

            if len_ref - len_qry >= 20:
                seqs.comments[qry.id] += "Large Deletion on {};".format(gene)

    with open(configs['path_out'] + '/summary_psc.tsv') as fh_i:
        # Header
        fh_i.readline()

        for line in fh_i:
            qid, gene, psc_type, _ = line.rstrip().split('\t')
            seqs.call[qid]['psc'] = 'Yes'
            seqs.info[qid]['psc'][gene] += ";{}".format(psc_type)

    for qid in seqs.qids:
        if seqs.call[qid]['deletion'] == 'Yes' or \
                seqs.call[qid]['inversion'] == 'Yes' or \
                seqs.call[qid]['hypermut'] == 'Yes' or \
                seqs.call[qid]['is_hiv'] == 'No':
            seqs.call[qid]['psc'] = 'Pass'
            seqs.info[qid]['psc'] = {'Gag': 'Pass', 'Pol': 'Pass',
                                     'Env': 'Pass'}
        elif 'psc' not in seqs.call[qid]:
            seqs.call[qid]['psc'] = 'No'
            if 'psc' not in seqs.info[qid]:
                seqs.info[qid]['psc'] = {'Gag': 'Pass', 'Pol': 'Pass',
                                         'Env': 'Pass'}

    file_cut = configs["path_out"] + "/Gene_Cutter/ALL.AA.PRINT"

    try:
        f_size = os.path.getsize(file_cut)
    except OSError:
        return

    if f_size == 0:
        return

    for qid in seqs.qids:
        if seqs.call[qid]['deletion'] == 'No' and \
                seqs.call[qid]['inversion'] == 'No' and \
                seqs.call[qid]['hypermut'] == 'No' and \
                seqs.call[qid]['is_hiv'] == 'Yes':
            file_out = configs["path_out"] + \
                "/Gene_Cutter/indv_reports/{}.txt".format(qid)
            with open(file_cut) as fh_i, open(file_out, 'w') as fh_o:
                for line in fh_i:
                    pat = "({})|(--)|(HXB2)".format(qid)
                    if re.match(pat, line):
                        print(line, file=fh_o, end='')
