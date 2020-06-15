"""
Check whether the primer sequences are in the sequenes.
"""

import logging
import re

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.primer')


def primer(configs, seqs):
    """Primer check.

    Check whether any of the primer sequences appear in the contig. The primer
    sequences are defined in Primer -> primers options.
    """
    # Logging
    logger.info('Checking primers')

    # Forward
    primers_fwd = re.split(';\n?', configs['primers'])

    # Reverse complement (translate then reverse the order)
    table = str.maketrans('ATGC', 'TACG')
    primers_rev = [p.translate(table)[::-1] for p in primers_fwd]

    # All the primer sequences
    primers = primers_fwd + primers_rev

    # Query data dict (qseqid -> SeqRecord)
    qry = seqs.qry

    # Search primer sequences in all contigs
    for qseqid in qry:
        primers_found = []
        for p in primers:
            indx = qry[qseqid].seq.find(p)
            if indx != -1:
                primers_found.append("{}({})".format(p, indx+1))

        if primers_found:
            seqs.call[qseqid]['primer'] = "Yes"
        else:
            seqs.call[qseqid]['primer'] = "No"

        seqs.info[qseqid]['primer'] = primers_found

    with open(configs['file_out'], 'w') as fh_o:
        # Header
        print("Contig ID\tSequence\tPrimer", file=fh_o)

        for qid in seqs.qids:
            print("{}\t{}\t{}".format(qid,
                                      ";".join(seqs.info[qid]['primer']),
                                      seqs.call[qid]['primer']),
                  file=fh_o)
