"""Some utilities"""

import sys
import subprocess


def find_gapped_pos(aln, ref_id=None, pos=0, nbp=1):
    """ Find gapped position given a position on the reference genome

    Parameters:
        aln: MultipleSeqAlignment obj
        ref_id: reference genome seq id
        pos: (ungapped) position on the reference genome
        nbp: number of basepairs
    """

    if ref_id is None:
        rec = aln
    else:
        # Find the seqrecord with reference genome id
        rec = next((r for r in aln if r.id == ref_id), None)

        if rec is None:
            sys.exit(1)

    locations = []
    j = 0
    for i, res in enumerate(rec.seq):
        if res != '-':
            if j == pos:
                locations.append(i)
            else:
                if locations:
                    locations.append(i)
                if len(locations) == nbp:
                    break
            j += 1

    return locations


def run_cmd(cmd):
    """Run a external command, raise Exception while failed.

    Parameters:
        cmd: comand line, including options.

    """

    # code = subprocess.call(cmd)
    res = subprocess.run(cmd, stdout=subprocess.PIPE)

    # if code:
    #     raise CommandFailedError(' '.join(cmd))

    return res
