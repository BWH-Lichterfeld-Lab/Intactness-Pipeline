# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:23:40 2017

@author: Rong Chen
# """

import os
import gc
import logging

from reportlab.lib import colors
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from PyPDF2 import PdfFileMerger, PdfFileReader

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.view')


class View:
    def __init__(self, configs, min_value=60, scale_fac=1):

        logger.info('Creating alignment views')

        annotation_file = configs['file_HIV_gene']
        blast_file = configs['file_blast']
        self.out = configs['path_out']

        if os.path.exists(self.out):
            msg = '{} already existed'.format(self.out)
            logger.info(msg)
        else:
            os.makedirs(self.out)

        self.start = int(configs['UTR_start'])
        self.min_value = min_value
        self.scale_fac = scale_fac
        self.read_annotation(annotation_file)
        msg = 'Got %d records from %s' % (len(self.annotation),
                                          annotation_file)
        logger.info(msg)
        self.ref_length = max([int(u[1]) for u in self.annotation])
        self.read_blast(blast_file)
        self.group = {}
        for rec in self.record:
            rep_name = rec.split('_')[0]
            if rep_name not in self.group:
                self.group[rep_name] = [rec, ]
            else:
                self.group[rep_name].append(rec)
        msg = 'There are %d query sequences with %d groups'
        msg = msg % (len(self.query_len), len(self.group))
        logger.info(msg)

    def scale_color(self, value):
        """Scale_color """

        if value < self.min_value:
            return self.scale_fac
        return 100.0 - (100.0-self.scale_fac) \
            / (100.0-self.min_value)\
            * (100.0-value)

    def read_annotation(self, filename):
        self.annotation = []
        with open(filename, 'U') as infile:
            for line in infile:
                tokens = line.rstrip().split('\t')
                self.annotation.append(tokens)

    def smart_switch(self, regions, _qlen):
        '''
        if the query has mulitiple blast hits, and half or more hits have
        reversed alignments with reference sequence, The alignments will be
        took reverse complement (rc)
        '''
        # r[-1] is the value of strand.
        no_reverse = [r[-1] for r in regions].count(-1)
        # no need to take rc
        if no_reverse * 2 < len(regions):
            return regions

        # otherwise it will take rc
        reversed_region = []
        for r in regions:
            new_qstart = _qlen-r[1] + 1
            new_qend = _qlen-r[0] + 1
            new_sstart, new_send = r[3], r[2]
            if r[5] == -1:
                new_strand = 1
            else:
                new_strand = -1
            reversed_region.append((new_qstart, new_qend, new_sstart,
                                    new_send, r[4], new_strand))
        return reversed_region

    def read_blast(self, filename):
        '''extracting blast hits'''
        # _record[query]: a list of the blast alignments for a query
        self.record = {}

        # the length for a query
        self.query_len = {}

        with open(filename, 'U') as infile:
            for line in infile:
                tokens = line.rstrip().split('\t')
                _query = tokens[0]  # query ID
                self.query_len[_query] = int(tokens[1])  # query length

                # if strand is plus, strand is set to 1, otherwise set to -1
                if tokens[-1] == 'plus':
                    strand = 1
                else:
                    strand = -1
                # blast alignments, including
                v = (int(tokens[4]), int(tokens[5]),  # qstart, qend,
                     int(tokens[6]), int(tokens[7]),  # sstart, send,
                     float(tokens[11]), strand)       # pident, and strand
                if _query not in self.record:
                    self.record[_query] = [v, ]
                else:
                    self.record[_query].append(v)

        if not self.record:
            print('No IDs in %s' % filename)
            exit(1)
        for _query in self.record:
            self.record[_query] = self.smart_switch(self.record[_query],
                                                    self.query_len[_query])

    def draw_alignment(self, _gdd, query_id, query_length, blast_hits):
        # draw reference
        gdt_features = _gdd.new_track(1, greytrack=False,
                                      start=1,
                                      end=self.ref_length)
        gds_features = gdt_features.new_set()
        max_length = 0
        for annot in self.annotation:
            if annot[-1] == '1':
                strand = 1
            elif annot[-1] == '-1':
                strand = -1
            else:
                strand = None
            start, end = int(annot[0]), int(annot[1])
            if end > max_length:
                max_length = end
            feature = SeqFeature(FeatureLocation(start, end), strand=strand,
                                 id=annot[2])
            gds_features.add_feature(feature, label=True, name=annot[2],
                                     sigil='BOX', label_size=14, label_angle=0,
                                     arrowhead_length=0.1, arrowshaft_height=1)
        gds_features.add_feature(SeqFeature(FeatureLocation(self.start,
                                                            self.start+1)),
                                 label=True, name='PCR start',
                                 label_size=14, color='black')
        _gdd.draw(format='linear', fragments=1, start=1, end=max_length)

        # draw query sequence
        gdt_features = _gdd.new_track(1, greytrack=True, start=0,
                                      end=query_length, name=query_id)
        gds_features = gdt_features.new_set()
        for r in blast_hits:
            feature = SeqFeature(FeatureLocation(r[0], r[1]), strand=r[-1])
            gds_features.add_feature(feature, label=False, sigil='BOX')
        gdt_features.greytrack_fontcolor = colors.black
        gdt_features.greytrack_fontsize = 12
        # draw cross link
        for r in blast_hits:
            alpha = self.scale_color(r[4])
            color = colors.linearlyInterpolatedColor(colors.white,
                                                     colors.firebrick, 0, 100,
                                                     alpha)
            feature = SeqFeature(FeatureLocation(r[0], r[1]), strand=r[-1])
            gds_features.add_feature(feature, label=False, sigil='BOX')
            if max(r[2], r[3]) < self.start:
                continue
            link_xy = CrossLink((_gdd.tracks[1], r[0], r[1]),
                                (_gdd.tracks[2], r[2], r[3]), color)
            _gdd.cross_track_links.append(link_xy)
        _gdd.draw(format='linear', fragments=1)

    def run(self):
        '''save figures to pdf'''
        for g in self.group:
            count = 1
            pdf_file = []
            pdf_title = []

            for r in self.group[g]:
                outfile = '%s_%d.pdf' % (g, count)
                gdd = GenomeDiagram.Diagram('Diagram of Blast Alignment')
                self.draw_alignment(gdd, r, self.query_len[r], self.record[r])
                gdd.write(outfile, output='pdf')
                pdf_file.append(outfile)
                pdf_title.append(r)
                count += 1
            merger = PdfFileMerger()
            count = 1
            for f, t in zip(pdf_file, pdf_title):
                merger.append(PdfFileReader(open(f, 'rb')), bookmark=t,
                              import_bookmarks=True)
                count += 1
            outfile = g + '.pdf'
            if self.out != '':
                outfile = os.path.join(self.out, outfile)
            merger.write(outfile)
            merger.close()
        gc.collect()
        for g in self.group:
            count = 1
            for r in self.group[g]:
                outfile = '%s_%d.pdf' % (g, count)
                # try:
                os.remove(outfile)
                # except:
                pass
                count += 1
