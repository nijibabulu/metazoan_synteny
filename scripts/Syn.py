#! /usr/bin/env python
import collections
import re
import FeatureIO

class SynBlockGraph(object):
    def __init__(self, blocks):
        """
        Represnt a synteny block graph by a edge dict mapping ids to SynBlocks.
        the clusters themselves are mapped by the cluster_id to a SynCluster.
        """
        self.v = {b.full_id: b for b in blocks}  # type:dict[str, SynBlock]
        self.e = collections.defaultdict(list)  # type:dict[SynBlock,[SynBlock]]

        map(self._connect, blocks)
        cs = (self._make_clusters())

        # let cmap map all blocks to their cluster
        self.cmap = {}  # type: dict[str, SynCluster]
        map(lambda c: map(lambda b: self.cmap.setdefault(b, c), c.blocks), cs)

    def _connect(self, block):
        # type: (SynBlock) -> None
        """Add all edges (connections) to the graph"""
        self.v[block.full_id] = block
        map(lambda cid: self.e[block].append(self.v[cid]), block.conn_ids)
        map(lambda cid: self.e[self.v[cid]].append(block), block.conn_ids)

    def _make_clusters(self):
        """
        perform a breadth first search on all nodes and yield a new cluster
        each time we have a non-empty set
        """
        visited = set()
        for startb in set(self.v.values()):
            if startb in visited:
                continue
            current_cluster = set()
            q = collections.deque([startb])
            while len(q):
                # this could be visited.
                # avoid check by making current_cluster a set
                b = q.popleft()
                visited.add(b)
                adj = [bb for bb in self.e[b] if bb not in visited]
                visited.update(adj)
                q += adj
                current_cluster.add(b)
            yield SynCluster(current_cluster)

    @property
    def clusters(self):
        # type: () -> list[SynCluster]
        return list(set(self.cmap.values()))


    @classmethod
    def from_oleg_file(cls, handle):
        # type: (file) -> SynBlockGraph
        return SynBlockGraph(map(SynBlock.from_oleg_format, handle))

    def annotate(self, annotation_map):
        map(lambda c: c.annotate(annotation_map), self.clusters)


class SynCluster(object):
    def __init__(self, blocks):
        self.id = ','.join(sorted(b.full_id for b in blocks))  # type: str
        self.blocks = blocks  # type: list[SynBlock]

    def add(self, block):
        self.blocks.append(block)

    def annotate(self, annotation_map):
        map(lambda b: b.annotate(annotation_map), self.blocks)

    @property
    def species(self):
        # type: () -> list[str]
        return list(set(b.species for b in self.blocks))

    def format(self):
        return '#{}\n{}'.format(
            self.id, '\n'.join(b.to_oleg_format() for b in self.blocks))


class SynBlock(object):
    def __init__(self, id, species, conn_ids, loc, gene_ids, names=None):
        self.id = id  # type: str
        self.species = species  # type: str
        self.conn_ids = conn_ids  # type: list[str]
        self.loc = loc  # type: FeatureIO.BaseFeature
        self.genes = gene_ids  # type: list[str]
        self.names = [''] * len(gene_ids) if names is None else names
        # type: list[str]

    @classmethod
    def from_oleg_format(cls, line):
        fields = line.strip().split('\t')
        conn_ids = filter(lambda cid: len(cid) > 0, fields[3].split(','))
        genes = fields[9].split(',')
        chrom, _, start, _, end = re.split(r'(:|\.\.)', fields[7])
        location = FeatureIO.BaseFeature(chrom, start, end)
        return SynBlock(fields[0], fields[1], conn_ids, location, genes)

    def annotate(self, annotation_map):
        self.names = [annotation_map.get(g,g) for g in self.genes]

    def to_oleg_format(self):
        fields = [self.id, self.species, len(self.conn_ids),
                  ','.join(self.conn_ids),
                  '{}:{}..{}'.format(self.loc.chrom, self.loc.start,
                                     self.loc.end), len(self.loc),
                  ','.join(self.genes),
                  ','.join(self.names)]
        return '\t'.join(map(str, fields))

    # @property
    # def cluster_id(self):
    #    return ','.join(sorted([self.full_id] + self.conn_ids))

    @property
    def full_id(self):
        return '{}({})'.format(self.id, self.species)
