"""
Check whether the start codon is present
"""

import os.path
import logging

from Bio import SeqIO

from .utils import find_gapped_pos

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.gag_codon')


def gag_codon(configs, seqs):
    """Primer check.

    Check whether any of the primer sequences appear in the contig. The primer
    sequences are defined in Primer -> primers options.
    """
    # Logging
    logger.info('Running codon analysis')

    # Output path
    file_aln = configs['file_aln']
    file_out = configs['file_out']
    pos_cdn = int(configs['pos'])

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

        pos = find_gapped_pos(ref, pos=pos_cdn-1, nbp=3)

        qry_str = str(qry.seq)
        qry_codon = ''.join(qry_str[i] for i in pos)

        seqs.info[qry.id]['gag_codon'] = qry_codon
        # Call of missing codon
        seqs.call[qry.id]['gag_codon'] = 'No' if qry_codon == 'ATG' else 'Yes'

    for qid in seqs.qids:
        # Large deletion or internal inversion or non-hiv
        if 'gag_codon' not in seqs.call[qid]:
            seqs.call[qid]['gag_codon'] = 'Pass'
            seqs.info[qid]['gag_codon'] = 'Pass'

    with open(file_out, 'w') as f_o:
        # Header
        print("SeqID\tCodon\tMissing", file=f_o)

        for qid in seqs.qids:
            print("{}\t{}\t{}".format(qid,
                                      seqs.info[qid]['gag_codon'],
                                      seqs.call[qid]['gag_codon']),
                  file=f_o)
