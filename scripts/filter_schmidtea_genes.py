#! /usr/bin/env python

import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
    description='''The schmidtea transcriptome has an extra field for the
    transcript id of a gene which is not included in the gene expression matrix.
    This script checks for genes included in the gene expression matrix
    regardless of their transcript id.''')
parser.add_argument('DGE_FILE')
parser.add_argument('PROTEIN_FILE')
args = parser.parse_args()


dge_genes = { l.split()[0] for l in open(args.DGE_FILE) }
seqs = [sr for sr in SeqIO.parse(open(args.PROTEIN_FILE), 'fasta')
        if '_'.join(sr.name.split('_')[:-1]) in dge_genes ]
SeqIO.write(seqs, sys.stdout, 'fasta')
