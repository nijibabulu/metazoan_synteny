#! /usr/bin/env python

import itertools

def complement_char(c):
    if c in complement_char.compd:
        return complement_char.compd[c]
    else:
        return c


complement_char.chars = 'acgtumrwsykvhdbnACGTUMRWSYKVHDBN'
complement_char.compl = 'tgcaakywsrmbdhvnTGCAAKYWSRMBDHVN'
complement_char.compd = dict(zip(complement_char.chars, complement_char.compl))


def reverse_complement(seq):
    return ''.join(complement_char(c) for c in seq)[::-1]


class ITreeNode(object):
    def __init__(self, i):
        self.i = i
        self._init_coords(i.start, i.end)

    def _init_coords(self,start,end):
        self.start = start
        self.end = end
        self.max = end
        self.left = None
        self.right = None

    @classmethod
    def from_coords(cls, start, end):
        obj = cls.__new__(cls)
        obj._init_coords(start, end)
        return obj


class ITree(object):
    def __init__(self, nodes=[]):
        self.root = None
        if len(nodes):
            # approximate a balanced tree by adding nodes in order of closeness to the center.
            center = max(n.start for n in nodes) / 2
            keyfunc = lambda n: abs(n.start - center)
            map(lambda n: self.insert(n), sorted(nodes, key=keyfunc))

    def insert(self, i):
        n = ITreeNode(i)
        if self.root is None:
            self.root = n
            return
        ncur = self.root
        while True:
            ncur.max = max(i.end, ncur.max)
            if i.start < ncur.start:
                if ncur.left is None:
                    ncur.left = n
                    break
                else:
                    ncur = ncur.left
            else:
                if ncur.right is None:
                    ncur.right = n
                    break
                else:
                    ncur = ncur.right

    def search(self, i):
        result = []
        if self.root is None:
            return result
        nstack = [self.root]
        while len(nstack):
            ncur = nstack.pop()
            if ncur.start <= i.end and i.start <= ncur.end:
                result.append(ncur.i)
            if ncur.left is not None and ncur.left.max >= i.start:
                nstack.append(ncur.left)
            if ncur.right is not None:
                nstack.append(ncur.right)
        return result


class BaseFeature(object):
    def __init__(self, chrom, start, end, strand='.'):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __len__(self):
        return self.end-self.start


class Gene(BaseFeature):
    # TODO: make a separate "from_blocks" instantiator for bed-style formats
    # separate from gff-style
    def __init__(self, chrom, start, end, name, score, strand, cds_start, cds_end,
                 item_rgb, block_count, block_sizes, block_starts, *args, **kwargs):
        super(Gene, self).__init__(chrom, start, end, strand)
        self.name = name
        self.score = score
        self.cds_start = int(cds_start)
        self.cds_end = int(cds_end)
        self.item_rgb = int(item_rgb)
        self.block_count = int(block_count)
        self.block_sizes = [int(s) for s in block_sizes.split(',') if len(s)]
        self.length = sum(self.block_sizes)
        self.block_starts = [int(s) for s in block_starts.split(',') if len(s)]

        self.__dict__.update(kwargs)

        self.exons = []
        self.cds_exons = []

        for start, size in zip(self.block_starts, self.block_sizes):
            start += self.start
            end = start + size
            sc = sorted([start, end, self.cds_start, self.cds_end])
            self.exons.append((start, end))
            if (sc[0] == start and sc[1] == end) or (
                    sc[0] == self.cds_start and sc[1] == self.cds_end):
                continue
            self.cds_exons.append((max(self.cds_start, start), min(self.cds_end, end)))

        self.length = sum([e[1] - e[0] + 1 for e in self.exons])
        self.cds_length = sum([e[1] - e[0] + 1 for e in self.cds_exons])

        if strand == '+':
            self.fivep = self.start
            self.cds_fivep = self.cds_start
        else:
            self.fivep = self.end
            self.cds_fivep = self.cds_end

    def _get_seq(self, seq, exons):
        gseq = ''.join(seq[c[0] - 1:c[1]]
                       for c in sorted(exons, cmp=lambda a, b: cmp(a[0], b[0])))
        if self.strand == '-':
            gseq = reverse_complement(gseq)
        return gseq

    def get_cds(self, seq):
        return self._get_seq(seq, self.cds_exons)

    def get_exons(self, seq):
        return self._get_seq(seq, self.exons)

    def locus_overlap(self, other):
        if self.start > other.end or other.start > self.end:
            return False
        if self.chrom != other.chrom:
            return False
        if self.strand != other.strand:
            return False
        return True

    def is_isoform(self, other, comparison='exons', check=False):
        if not self.locus_overlap(other):
            return False
        for a, b in itertools.product(getattr(self,comparison),
                                      getattr(other,comparison)):
            if a[0] == b[0] and a[1] == b[1]:
                return True
        return False

    def overlap(self, other, comparison='exons'):
        if not self.locus_overlap(other):
            return False
        for a, b in itertools.product(getattr(self, comparison),
                                      getattr(other, comparison)):
            if a[0] <= b[1] and b[0] <= a[1]:
                return True
        return False

    def overlap_length(self, other, comparison='exons'):
        if not self.locus_overlap(other):
            return False
        length = 0
        for a, b in itertools.product(getattr(self, comparison),
                                      getattr(other, comparison)):
            length += max(0,min(a[1],b[1])-max(a[0],b[0]))
        return length

    def __str__(self):
        return '\t'.join(str(item) for item in
                         [self.chrom, self.start, self.end, self.name, self.score, self.strand,
                          self.cds_start, self.cds_end, self.item_rgb, self.block_count,
                          ','.join(str(s) for s in self.block_sizes),
                          ','.join(str(s) for s in self.block_starts)])


def BedIterator(handle,cls=Gene):
    for line in handle:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) == 1:
            continue
        yield cls(*fields)


def parse_psl_line(matches, misMatches, repMatches, nCount, qNumInsert,
                   qBaseInsert, tNumInsert, tBaseInsert, strand,
                   name, qSize, qStart, qEnd, chrom, tSize, start, end, block_count,
                   block_sizes, qStarts, block_starts,cls=Gene):
    corr_block_starts = ','.join([str(int(s) - int(start))
                                  for s in block_starts.split(',') if len(s)])
    return cls(chrom, start, end, name, matches, strand, 0, 0, 0,
                block_count, block_sizes, corr_block_starts)


def PslIterator(handle,cls=Gene):
    for line in handle:
        fields = line.strip().split()
        yield parse_psl_line(*fields,cls=cls)


def BlatPslIterator(handle,cls=Gene):
    for _ in range(5):  # skip the header
        next(handle)
    for line in handle:
        yield parse_psl_line(line,cls=cls)


def AugustusGtfIterator(handle,cls=Gene):
    while True:
        line = handle.readline()
        if not line:
            return
        if line.startswith("# start gene"):
            break
    while True:
        line = handle.readline()
        if not line:
            return  # stopiteration
        f = line.strip().split('\t')
        if len(f) > 1 and f[2] == 'transcript':
            chrom, start, end, strand, tid = (f[0], int(f[3]), int(f[4]), f[6], f[-1])
            bst, bsz = [[], []]
            cs, ce = [None, None]
            seq = ''
            gene_id = None
            while True:
                line = handle.readline()
                if not line:
                    return  # stopiteration
                # if we have already seen protein sequence or we are initiating
                # the protein sequence add to the current sequence
                if line.startswith('# protein sequence') or len(seq):
                    istart = (line.index('[') + 1) if line.count('[') else 2
                    iend = line.index(']') if line.count(']') else len(line)
                    seq += line[istart:iend].strip()
                else:
                    # TODO: change to a non-block initializer to avoid these
                    # calcluations
                    f = line.strip().split('\t')
                    if f[2] == 'exon':
                        gene_id = f[-1].split()[-1].strip('";')
                        bst.append(int(f[3]) - start)
                        bsz.append(int(f[4]) - int(f[3]))
                    elif f[2] == 'CDS':
                        gene_id = f[-1].split()[-1].strip('";')
                        cs = min(cs, int(f[3])) if cs else int(f[3])
                        ce = max(cs, int(f[4])) if ce else int(f[4])
                if seq is not None and ']' in line:
                    yield cls(
                        chrom, start, end, tid, 0, strand, cs, ce, 0, len(bst),
                        ','.join(str(x) for x in bsz),
                        ','.join(str(x) for x in bst),
                        seq=seq, gene_id=gene_id)
                    break


_readers = {"bed12": BedIterator, "psl": PslIterator, "blatpsl": BlatPslIterator,
            "augustusgtf": AugustusGtfIterator}


def parse(maybe_handle, format, mode='r', cls=Gene, **kwargs):
    # type: (handle_or_string, basestring, basestring, ...) -> list
    # this can be better handled with contextlib.contextmanager
    if isinstance(maybe_handle, basestring):
        fp = open(maybe_handle, mode, **kwargs)
    else:
        fp = maybe_handle

    if format in _readers:
        gen = _readers[format]
        i = gen(fp, cls=cls)

        for gene in i:
            yield gene
    else:
        raise ValueError('Unknown format {}. Should be one of {}'.format(
            format, ','.join(_readers.keys())))


class GeneWriter(object):
    def __init__(self, handle):
        self.handle = handle

    def write_header(self):
        pass

    def write_genes(self, genes):
        for i, gene in enumerate(genes):
            self.count = i + 1
            self.write_gene(gene)

    def write_gene(self, _):
        raise NotImplementedError('{} does not implement write_gene'.format(
            self.__class__))

    def write_footer(self):
        pass

    def write_file(self, genes):
        self.write_header()
        for gene in genes:
            self.write_gene(gene)
        self.write_footer()


class AugustusExonHintWriter(GeneWriter):
    def __init__(self, handle, cds_exons=True, feature_type='exon',
                 source='featureio', augustus_source='E', priority=4):
        super(AugustusExonHintWriter, self).__init__(handle)
        self.exon_attr = 'cds_exons' if cds_exons else 'exons'
        self.feature_type = feature_type
        self.source = source
        self.priority = priority
        self.augustus_source = augustus_source

    def write_gene(self, gene):
        for exon in getattr(gene, self.exon_attr):
            attrs = 'grp={};pri={};src={}'.format(
                gene.name, self.priority, self.augustus_source)
            self.handle.write('\t'.join(str(x) for x in [
                gene.chrom,
                self.source,
                self.feature_type,
                exon[0],
                exon[1],
                '.',
                gene.strand,
                '.',
                attrs]))
            self.handle.write('\n')


class Bed12Writer(GeneWriter):
    def write_gene(self, gene):
        self.handle.write('\t'.join(str(item) for item in
                                    [gene.chrom, gene.start, gene.end, gene.name, gene.score, gene.strand,
                                     gene.cds_start, gene.cds_end, gene.item_rgb, gene.block_count,
                                     ','.join(str(s) for s in gene.block_sizes),
                                     ','.join(str(s) for s in gene.block_starts)]) + '\n')


_writers = {"bed12": Bed12Writer,
            "augustus_exon_hints": AugustusExonHintWriter}


def write(genes, maybe_handle, format, mode='w', **kwargs):
    if isinstance(maybe_handle, basestring):
        fp = open(maybe_handle, mode, **kwargs)
    else:
        fp = maybe_handle

    if format in _writers:
        writer = _writers[format](fp, **kwargs)
        writer.write_file(genes)
    else:
        raise ValueError('Unknown format {}. Should be one of {}'.format(
            format, ','.join(_writers.keys())))


if __name__ == '__main__':
    from Bio import SeqIO

    seq = str(SeqIO.read('unmask_split_nv2i5/Chr1.fa', 'fasta').seq)

    for gene in parse('bookends_rb/Chr1.BE0.1.gff', 'augustusgtf'):
        print gene.name
        print gene.strand
        cds = gene.get_cds(seq)
        print len(cds)
        print gene.cds_length
        print cds
