# coding: utf-8
import time
import re
import sys
import logging
import zipfile
import os

import mechanicalsoup

from datetime import datetime
from datetime import datetime
from random import randrange
from collections import defaultdict

# Disable warning
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe.GeneCutter')

URL_BASE = "https://www.hiv.lanl.gov/"

def sleep_btw(early, late):
    time.sleep(randrange(early, late))

def submit_GC(email_address):
    """
    Submit aligned seqs to Gene Cutter and start the job

    """

    # Create browser object
    br = mechanicalsoup.StatefulBrowser()

    # Open website
    url_gene_cutter = URL_BASE + "content/sequence/GENE_CUTTER/cutter.html"
    br.open(url_gene_cutter, verify = False)

    logger.info('Opening Gene Cutter website')
    sleep_btw(0, 5)

    # Upload sequences
    br.select_form('form[method="post"]')
    br['INSERTSTDSEQ'] = 'YES'
    br['UPLOAD'] = 'data/seqs/seqs_psc.fasta'
    response = br.submit_selected()

    logger.info('Uploading sequences to Gene Cutter')
    sleep_btw(0, 5)

    # Fill in contact information and start the job
    # TODO use custom email
    br.select_form('form[method="post"]')
    br['titleFromUser'] = 'PSC ' + datetime.now().strftime("%Y-%m-%d %H:%M")
    br['EMAIL'] = email_address
    br['EMAIL2'] = email_address
    response = br.submit_selected()

    logger.info('Submitting job to Gene Cutter')
    sleep_btw(0, 5)

    # Parse download url
    job_id = re.search(r'<b>/tmp/GENE_CUTTER/(.*)</b>',
            response.content.decode()).group(1)

    logger.info('Gene Cutter Job ID: {}'.format(job_id))
    sleep_btw(0, 5)

    url_download = f'https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/tmp/GENE_CUTTER/{job_id}/genecutter.zip'
    url_all = f'https://www.hiv.lanl.gov/tmp/GENE_CUTTER/{job_id}/all_aa.html'

    # Check the result availability
    while True:
        # TODO Progressively increate wait time
        sleep_btw(60, 61)
        response = br.get(url_download, verify=False)

        if b'Request Rejected' not in response.content:
            with open('data/seqs/genecutter.zip', 'wb') as fh:
                fh.write(response.content)
            fnz = zipfile.ZipFile('data/seqs/genecutter.zip', 'r')
            fnz.extractall('data/seqs/')
            fnz.close()

            os.rename('data/seqs/genecutter', 'data/seqs/Gene_Cutter')
            os.mkdir('data/seqs/Gene_Cutter/indv_reports')

            response = br.get(url_all, verify=False)

            with open('data/seqs/Gene_Cutter/ALL.AA.PRINT', 'w') as out:
                content = re.sub(r'<br>', '\n', response.content.decode())
                content = re.sub(r'<.*?>', '', content)
                out.write(content)

            break


def process_GC():
    # Parse results
    results = defaultdict(dict)

    gene_set = set(['Gag', 'Pol', 'Env'])
    with open('data/seqs/Gene_Cutter/ALL.AA.PRINT') as fn:
        line = fn.readline()
        while line:
            if line.startswith('---------- List of Stop Codons Within Sequences'):
                gene = re.search('^.*\((.*)\).*$', line).group(1)
                if gene in gene_set:
                    line = fn.readline()

                    while True:
                        if line == '\n':
                            line = fn.readline()
                        elif line.startswith('----'):
                            break
                        else:
                            line = line.strip()
                            seq_id, _, pos = line.split(' ')
                            results[gene].setdefault(seq_id, set()).add((int(pos), 'SC', ''))
                            line = fn.readline()
                else:
                    line = fn.readline()

            elif line.startswith('---------- List of Incomplete Codons'):
                gene = re.search('^.*\((.*)\).*$', line).group(1)
                if gene in gene_set:
                    line = fn.readline()

                    while True:
                        if line == '\n':
                            line = fn.readline()
                        elif line.startswith('----'):
                            break
                        else:
                            line = line.strip()
                            seq_id, _, pos, ic_type = line.split(' ')
                            results[gene].setdefault(seq_id, set()).add((int(pos), 'IC', ic_type))
                            line = fn.readline()
                else:
                    line = fn.readline()
            else:
                line = fn.readline()

    if (len(results) == 0):
        with open('data/seqs/summary_psc.tsv', 'w') as fn:
            print('Contig\tRef\tType\tPSC', file=fn)
    else:
        final_results = defaultdict(dict)
        for gene, contigs in results.items():
            for contig, events in contigs.items():
                events = sorted(events)
                idx_used = []
                for idx, event in enumerate(events):
                    if idx <= 1:
                        final_results.setdefault(gene, {}).setdefault(contig, {}).setdefault(event[1], []).append((event[0], event[2]))
                    elif events[idx][1] == 'IC' and \
                        events[idx - 1][1] == 'SC' and \
                        events[idx - 2][1] == 'IC' and \
                        events[idx-1][0] - events[idx-2][0] < 100 and \
                        events[idx][0]   - events[idx-1][0] < 100:
                        del final_results[gene][contig]['IC'][-1]
                        del final_results[gene][contig]['SC'][-1]
                    else:
                        final_results.setdefault(gene, {}).setdefault(contig, {}).setdefault(event[1], []).append((event[0], event[2]))

        with open('data/seqs/summary_psc.tsv', 'w') as fn:
            print('Contig\tRef\tType\tPSC', file=fn)
            # Remove beginning and ending
            for gene, contigs in final_results.items():
                for contig, event_types in contigs.items():
                    for event_type, events in event_types.items():
                        events = [event for event in events if event[0] != 1]
                        if (event_type == 'SC'):
                            events = ";".join(str(event[0]) for event in events)
                        else:
                            events = ";".join(str(event[0])+ ":" + event[1] for event in events)
                        if events == '':
                            continue
                        print("{}\t{}\t{}({})\tYes".format(contig, gene, event_type, events), file=fn)


if __name__ == '__main__':
    submit_GC()
    process_GC()
