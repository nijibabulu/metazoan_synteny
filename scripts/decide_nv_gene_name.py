#! /usr/bin/env python

import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('NEMVE_DUAL_BED')
    parser.add_argument('UMI_TABLE', nargs='+')
    args = parser.parse_args()

    gene_names = set(line.split('\t')[0] 
                     for file in args.UMI_TABLE
                     for line in open(file))

    with open(args.NEMVE_DUAL_BED) as f:
        for line in f:
            chrom, start, end, name, color, strand = line.strip().split('\t')
            names = name.split(',')
            included = [n in gene_names for n in names]
            if sum(included) > 1: 
                raise ValueError, 'Both names are included ' + names
            elif sum(included) == 0 or len(names) == 1:
                table_name = names[0]
            else:
                table_name = names[included.index(True)]
            if not table_name.startswith('EMNVE'):
                print '\t'.join([chrom, start, end, 'nv_'+table_name, color, strand])
