#! /usr/bin/env python3

import os
import csv
import click
import itertools
from dataclasses import dataclass,field
from typing import Dict, List, Tuple

@dataclass
class BlockPartition:
    name: str
    nonparalogous: int = 0
    paralogous: int = 0

    def __add__(self, other):
        return BlockPartition(self.nonparalogous + other.nonparalogous,
                              self.paralogous + other.paralogous)

    def paralogous_fraction(self):
        return self.paralogous / (self.paralogous + self.nonparalogous)

    def nonparalogous_fraction(self):
        return self.nonparalogous / (self.paralogous + self.nonparalogous)


@dataclass
class OrthoGroup:
    name: str
    genes: List[str]
    orthologs: List[Tuple[str,str]] = field(default_factory=list)
    _paralogs: List[str] = None


    @property
    def paralogs(self):
        if self._paralogs is None:
            all_orthologs = itertools.chain.from_iterable(self.orthologs)
            self._paralogs = set(self.genes) - set(all_orthologs)
        return self._paralogs

@click.command(help="Find any paralogs in a orthofinder-formatted file "
               "(ORTHOGROUPS) within synteny blocks (BLOCKS)")
@click.argument("ORTHOGROUPS")
@click.argument("ORTHOLOGS_DIR")
@click.argument("BLOCKS")
def main(orthogroups, orthologs_dir, blocks):
    with open(orthogroups) as f:
        og_csv = csv.reader(f, delimiter=' ')
        groups = [OrthoGroup(r[0].rstrip(':'), r[1:]) for r in og_csv]
        group_dict = {k: v  for v in groups for k in v.genes}

    for root, dirs, files in os.walk(orthologs_dir):
        for file in files:
            with open(os.path.join(root,file)) as f:
                ortho_csv = csv.reader(f, delimiter=',')
                _ = next(ortho_csv, None)
                for (og,o1,o2) in ortho_csv:
                    # ignore co-orthologs and 
                    if ',' in o1 or ',' in o2:
                        continue
                    group_dict[o1].orthologs.append((o1,o2))

    with open(blocks) as f:
        d: Dict[str, BlockPartition] = dict()
        for l in f:
            if l.startswith("#"):
                continue
            fields = l.strip().split("\t")
            name = f"{fields[0]}({fields[1]})"
            block = fields[6].split(",")
            names = fields[7].split(",")

            bp = BlockPartition(name)
            for i,j in itertools.combinations(range(len(block)), 2):
                ga, gb = block[i],block[j]
                if ga in group_dict and gb in group_dict and \
                   group_dict.get(ga) == group_dict.get(gb) and \
                   ga in group_dict.get(ga).paralogs and \
                   gb in group_dict.get(gb).paralogs:
                    bp.paralogous += 1
                else:
                    bp.nonparalogous += 1
            d.setdefault(block[0][:2], []).append(bp)

    for species in d.keys():
        print("\n".join(f"{species}\t{bp.name}\t{bp.nonparalogous_fraction()}"
                        for bp in d[species]))

if __name__ == "__main__":
    main()
