"""Container for sequences"""

import logging
from collections import defaultdict

from Bio import SeqIO
import numpy as np

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.sequence')


def get_file_handle_by_type(file_name, _mode, _type):
    """Open a file according to its mode and return the handle
    """

    _handle = None
    _type = int(_type)
    if _type == 0:
        _handle = open(file_name, _mode)
    elif _type == 1:
        import gzip
        _handle = gzip.open(file_name, _mode)
    elif _type == 2:
        import bz2
        _handle = bz2.BZ2File(file_name, _mode)
    return _handle


def get_seqs(configs, verbose=True):
    """Get sequences according configuration"""
    f_i = configs['file_seq']
    fmt = configs['format']
    compress = configs['compress']

    fh_i = get_file_handle_by_type(f_i, 'r', compress)
    seqs = SeqIO.to_dict(SeqIO.parse(fh_i, fmt))

    if verbose:
        # Print length distribution
        msg = 'There are {} sequences in {}'.format(len(seqs), f_i)
        logger.info(msg)

        # Get Returns maximum, minimum, median, mean, and stdev
        lens = [len(s) for s in seqs.values()]

        len_msg = ['Seq      min:\t{}'.format(min(lens)),
                   'Seq      max:\t{}'.format(max(lens)),
                   'Seq   median:\t{:.2f}'.format(np.median(lens)),
                   'Seq mean/std:\t{:.2f}/{:.2f}'.format(np.mean(lens),
                                                         np.std(lens))]
        for msg in len_msg:
            logger.info(msg)

    return seqs


class Sequences:
    """Container for sequences necessary for the pipeline to run properly.

    Attributes
    ----------
    qry : dict
        Query sequences
    ref : dict
        Reference sequence
    call: dict
        Various call for each contig
    info: dict
        Associated info for each contig

    """

    # disable=too-many-instance-attribute
    # pylint: disable=R0902

    def __init__(self, cfg_qry, cfg_ref):
        logger.info('Importing reference genome and query contigs...')

        self.qry = get_seqs(cfg_qry)
        self.ref = get_seqs(cfg_ref, verbose=False)
        self.call = defaultdict(dict)
        self.info = defaultdict(dict)
        self.comments = defaultdict(str)

    @property
    def ref_id(self):
        """Reference genome id"""
        return list(self.ref.keys())[0]

    @property
    def qids(self):
        """Query seq id"""
        return list(self.qry.keys())

    @property
    def ref_seq(self):
        """Return Seq obj for the reference"""

        return self.ref[self.ref_id].seq

    def qlen(self, qseqid):
        """Get the sequence length of a particular query"""
        return len(self.qry[qseqid].seq)

    def write(self, f_o, qids=None, fmt='fasta', prepend_ref=False):
        """Write sequence

        Parameters
        ----------
        ids: str, IDs of the seqs to be written
        f_o: str, output file name
        fmt: str, output file format
        prepend_ref: logical
            whether to prepend reference geneome to the beginning
        """

        if qids is None:
            qids = self.qry.keys()

        if prepend_ref:
            records = list(self.ref.values())
        else:
            records = []

        records += [self.qry[seq_id] for seq_id in qids]

        SeqIO.write(records, f_o, fmt)

    def take_rc_for(self, qids):
        """Create rc for seqs"""
        for qid in qids:
            self.qry[qid].seq = self.qry[qid].seq.reverse_complement()
