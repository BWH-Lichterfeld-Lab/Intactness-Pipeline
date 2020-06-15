#!/usr/bin/env python3
"""
multiple sequence alignment using muscle
"""

from .utils import run_cmd


def muscle(configs):
    """Multiple sequence alignment using muscle
    """

    file_i = configs['file_seq']
    file_o = configs['file_aln']
    maxiters = configs['maxiters']

    cmd = ['muscle', '-in', file_i, '-out', file_o, '-maxiters', maxiters]
    run_cmd(cmd)
