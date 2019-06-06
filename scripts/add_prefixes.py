#! /usr/bin/env python

import argparse


class Prefixer(object):
    def __init__(self, pairs):
        self.pairs = pairs # type: list[(str, str)]

    def _prefix_geneid(self, ortholog):
        match = [p for p in self.pairs if ortholog.startswith(p[0])]
        if len(match) > 1:
            raise ValueError('Multiple matches for gene {}'.format(ortholog))
        if len(match) < 1:
            raise ValueError('No match for gene {}'.format(ortholog))
        return '{}_{}'.format(match[0][1], ortholog)

    def __call__(self, line):
        raise NotImplementedError()


class OrthofinderPrefixer(Prefixer):
    def __call__(self, line):
        if line.startswith('#'):
            return line
        else:
            fields = line.split()
            return ' '.join([fields[0]] +
                            map(lambda g: self._prefix_geneid(g), fields[1:]))


class EggnogPrefixer(Prefixer):
    def __call__(self, line):
        if line.startswith('#'):
            return line.strip()
        else:
            fields = line.strip().split('\t')
            return '\t'.join([self._prefix_geneid(fields[0])] + fields[1:])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    prefixers = {'eggnog': EggnogPrefixer,
                 'orthofinder': OrthofinderPrefixer}
    parser.add_argument('FORMAT', choices=prefixers.keys())
    parser.add_argument('FILE')
    parser.add_argument('PATTERN_PREFIX', nargs='+',
                        help='takes the format <id_prefix>:<prefix>. add '
                             'prefix to all sequences that begin with '
                             'id_prefix.')
    args = parser.parse_args()

    prefixer = prefixers[args.FORMAT]([p.split(':') for p in args.PATTERN_PREFIX])
    with open(args.FILE) as f:
        print '\n'.join(map(prefixer, f))
