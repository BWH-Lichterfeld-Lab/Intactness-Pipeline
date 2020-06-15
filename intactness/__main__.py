# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 00:29:42 2017

@author: Ce Gao
@author: Rong Chen
"""
# pylint: disable=C0103
# Disable module level constant name checking

import sys
import logging

# Local module
from .view import View
from .sequence import Sequences
from .configs import configs
from .blast import blast
from .muscle import muscle
from .gag_codon import gag_codon
from .primer import primer
from .hypermut import hypermut
from .psc import psc
from .defect import defect
from .summary import summary


logger = logging.getLogger('pipe')
logger.setLevel(logging.INFO)

import os
mod_path = os.path.dirname(os.path.abspath(__file__))
cfg = configs(os.path.join(mod_path, 'default.cfg'))

seqs = Sequences(cfg['Query'], cfg['Reference'])

primer(cfg['Primer'], seqs)

blast(cfg['BLAST'], seqs)

View(cfg['View']).run()

muscle(cfg['MSA'])

gag_codon(cfg['Codon'], seqs)

hypermut(cfg['Hypermut'], seqs)

psc(cfg['PSC'], seqs)

defect(cfg['Defect'], seqs)

summary(cfg['Summary'], seqs)

print('\nPipeline finished. Good bye!\n')
