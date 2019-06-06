#! /usr/bin/env python

import os
import operator
import io
import argparse

import Syn

__doc__ = """
Parse a synteny cluster file and output blocks ordered by a cluster and filtering
by certain properties of the clusters.
"""

class GeneInfo(object):
    def __init__(self, gene_id, names):
        self.gene_id = gene_id  # type: str
        self.names = names  # type: list[str]

    def __getitem__(self, item):
        return self.names.__getitem__(item)


class EmapperInfo(GeneInfo):
    @classmethod
    def from_line(cls, line):
        if line.startswith('#'):
            return None
        fields = line.strip().split('\t')
        short_gene_name = fields[4].rstrip('0123456789-')

        # long gene name can be 'NA', in this case there's still mineable
        # information in the OGs, but we're leaving that alone for now.

        long_gene_name = fields[12] if fields[12] != 'NA' else ''
        return EmapperInfo(fields[0], [short_gene_name, long_gene_name])


class TrinotateInfo(GeneInfo):
    @classmethod
    def from_line(cls, line):
        if line.startswith('#'):
            return None
        fields = line.strip().split('\t')
        if fields[2] == '.':
            return None

        unpacked_hits = [h.split('^') for h in fields[2].split('`')]

        # create a tuple of evalue, annotation
        # evalue format is e.g. 'E:1e-3'
        # gene name format is <GENENAME>_<SPECIES>
        hits = [(float(hfs[4][2:]),
                 hfs[0].split('_')[0].rstrip('0123456789-')
                 ) for hfs in unpacked_hits]

        # get the top hit by the lowest evalue
        best = sorted(hits, key=operator.itemgetter(0))[0][1]
        return TrinotateInfo(fields[0], [best])


def parse_name(line):
    fields = line.strip().split('\t')
    return fields[0], fields[12]


class AnnotationMap(object):
    def __init__(self, emapper_info, trinotate_info):
        self.emapper_info = emapper_info
        self.trinotate_info = trinotate_info

    def __getitem__(self, item):
        e_hits = self.emapper_info.get(item, ['', ''])
        t_hit = self.trinotate_info.get(item, [''])
        return e_hits[0] or t_hit[0] or e_hits[1]

    def get(self, item, default):
        return self[item] or default


def get_annotation_map(emapper_out=None, trinotate_out=None):
    emapper_map = {ei.gene_id: ei
                   for ei in map(EmapperInfo.from_line, emapper_out)
                   if ei is not None}

    trinotate_map = {ti.gene_id: ti for ti in
                     map(TrinotateInfo.from_line, trinotate_out)
                     if ti is not None}

    return AnnotationMap(emapper_map, trinotate_map)


class BlockFileCache(object):
    def __init__(self, dir, suffix='.clust'):
        self.cache = {}
        self.dir = dir
        self.suffix = suffix

    def get_name(self, species):
        return '-'.join(sorted(species)) + self.suffix

    def __call__(self, cluster):
        name = self.get_name(list(cluster.species))
        if name not in self.cache:
            self.cache[name] = open(os.path.join(self.dir, name), 'w')
        return self.cache[name]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--num-species', type=int, metavar='N',
                        help='output blocks with at least N species only')
    parser.add_argument('--species-combo-output', metavar='DIR',
                        help='separate file for each combination of species '
                             'and output files into DIR.')
    parser.add_argument('CLUST_FILE')
    parser.add_argument('EMAPPER_OUT', type=argparse.FileType('r'), nargs='?',
                        default=io.StringIO())
    parser.add_argument('TRINOTATE_OUT', type=argparse.FileType('r'), nargs='?',
                        default=io.StringIO())
    parser.add_argument('TWOCOLUMN', nargs='?', help='Optionally give a two '
                        'column input file consisting of genes and annotations '
                        'separated by tabs.')
    args = parser.parse_args()

    with open(args.CLUST_FILE) as f:
        graph = Syn.SynBlockGraph.from_oleg_file(f)

    ann_map = get_annotation_map(args.EMAPPER_OUT, args.TRINOTATE_OUT)
    graph.annotate(ann_map)

    clusters = graph.clusters
    if args.num_species is not None:
        clusters = filter(lambda c: len(c.species) >= args.num_species, clusters)

    if args.species_combo_output:
        cache = BlockFileCache(args.species_combo_output)
        for b in clusters:
            cache(b).write(b.format() + '\n')

        # map(lambda b: cache(b).write(b.format() + '\n'), clusters)

    print '\n'.join(b.format() for b in clusters)
