"""5' defects

"""
import logging
import os.path

from Bio import SeqIO

from .utils import find_gapped_pos

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.defect')


def _detect_defect(ref, qry, pos_start, pos_end):

    # Find tag position in reference seq, and start aln pos for the rest.
    pos_start = find_gapped_pos(ref, pos=pos_start - 1)[0]
    pos_end = find_gapped_pos(ref, pos=pos_end - 1)[0]

    # Iterate through records id, start end
    n_del = 0
    n_ins = 0
    for char_ref, char_qry in zip(ref.seq[pos_start:pos_end],
                                  qry.seq[pos_start:pos_end]):
        if (char_ref != '-') and (char_qry == '-'):
            n_del += 1
        elif (char_ref == '-') and (char_qry != '-'):
            n_ins += 1
    return n_del, n_ins


def defect(configs, seqs):
    """Determine if a contig has 5' defect in the tag region.
    """
    # Logging
    logger.info("Checking 5' defects")

    # Path to codon alignment
    file_aln = configs['file_aln']

    try:
        f_size = os.path.getsize(file_aln)
    except OSError:
        f_size = 0

    if f_size > 0:
        # Threshold of gap number to determine 5' defect
        max_gaps = int(configs['max_gaps'])

        pos_start = int(configs['start'])
        pos_end = int(configs['end'])

        aln = SeqIO.parse(file_aln, 'fasta')
        ref = next(aln)
        for qry in aln:
            qid = qry.id.rstrip('_Genome')
            n_del, n_ins = _detect_defect(ref, qry, pos_start, pos_end)

            if n_del >= max_gaps:
                if seqs.call[qid]['primer'] == 'No':
                    seqs.call[qid]['defect'] = 'Possible'
                    seqs.comments[qid] += "Missing primer for 5' defect;"
                else:
                    seqs.call[qid]['defect'] = 'Yes'
            else:
                seqs.call[qid]['defect'] = 'No'

            seqs.info[qid]['defect'] = (n_del, n_ins)

    # Run hypermut for long contigs
    for qid in seqs.qids:
        if 'defect' not in seqs.call[qid]:
            seqs.call[qid]['defect'] = 'Pass'
            seqs.info[qid]['defect'] = ['Pass', 'Pass']

    # Write output
    with open(configs['file_out'], 'w') as fh_o:
        # Header line
        print("Contig ID\tGaps\tInserts\t5' Defect", file=fh_o)

        for qid in seqs.qids:
            print('{}\t{}\t{}\t{}'.format(qid, *seqs.info[qid]['defect'],
                                          seqs.call[qid]['defect']),
                  file=fh_o)
