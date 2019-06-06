#! /usr/bin/env python

import sys
import re
import itertools
import argparse
from Bio import SeqIO


def convert_jgi_rec(rec):
    match = re.search(r'gene:NEMVEDRAFT_(?P<geneid>\S+)\s', rec.description)
    if match is not None:
        rec.id = match.group('geneid')
    return rec


def add_nv_prefix(rec):
    rec.id = 'nv_' + rec.id
    rec.description = 'nv_' + rec.description
    return rec


def parse_bed_gene_ids(filename):
    with open(filename) as f:
        return [line.strip().split()[3] for line in f]


def add_seq(jgi_dict, nve_dict, acc, id):
    if ',' in id: # one of the JGI genes was joined
        acc += [jgi_dict[k] for k in id.split(',')]
    elif id.startswith('NVE'):
        acc.append(nve_db[id])
    elif id.startswith('v1'):
        acc.append(jgi_db[id])
    elif id.startswith('EMNVEG'):
        pass 
    else:
        raise KeyError, 'Unknown gene id: {}'.format(gene.name)
    return acc


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('JGI_FASTA')
    parser.add_argument('NVE_FASTA')
    parser.add_argument('BED')
    args = parser.parse_args()

    jgi_recs = map(convert_jgi_rec, SeqIO.parse(args.JGI_FASTA,'fasta'))
    jgi_db = {g.id: g for g in jgi_recs}
    nve_db = SeqIO.to_dict(SeqIO.parse(args.NVE_FASTA,'fasta'))

    ids = set(parse_bed_gene_ids(args.BED))

    SeqIO.write(map(add_nv_prefix, 
                    reduce(lambda g,acc: add_seq(jgi_db, nve_db, g, acc), ids, [])),
                sys.stdout, 'fasta')
