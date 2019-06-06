#! /usr/bin/env python

import sys
import argparse
import FeatureIO
from Bio import SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('GENOME')
    parser.add_argument('GENEFILE',nargs='+')
    parser.add_argument('--format',default='bed12',
                        help='format of gene file [bed12]')
    parser.add_argument('--cds', action='store_true',
                        help='retrive cds sequence only')
    # TODO: custom prefixes?

    args = parser.parse_args()

    genome = SeqIO.to_dict(SeqIO.parse(args.GENOME,'fasta'))
    print >> sys.stderr, "got genome"

    seq_method = 'get_cds' if args.cds else 'get_exons'

    for aug in args.GENEFILE:
        for gene in FeatureIO.parse(aug,args.format):
            gene_seq = getattr(gene,seq_method)(str(genome[gene.chrom].seq))
            print '>{}\n{}'.format(gene.name, gene_seq)
