# Gene Functional Annotation

This describes how we used [EggNog Mapper](https://github.com/eggnogdb/eggnog-mapper) and [Trinotate](https://github.com/Trinotate/Trinotate) to identify the genes in the synteny blocks.  Emapper is the more critical step of the two, and Trinotate can be avoided.

## Emapper 

Generate an emapper annotations for all genes:

```bash
cd processing/05_emapper
for s in aq ml ta nv sm; do 
    emapper.py --override --cpu 1 --database meNOG --report_ortholog  --output_dir out -i ../01_dbs/$s.faa
done
```

Concatenate the results, add prefixes and make readable blocks files:

```bash
cat processing/05_emapper/out/*.annotations > processing/05_emapper/emapper.annotations.txt
python scripts/add_prefixes.py eggnog processing/05_emapper/emapper.annotations.txt Aqu:aq ML:ml v1g:nv NVE:nv Tadh:ta > processing/05_emapper/emapper.annotations.prefixed.txt
```

## Trinotate 

Add species prefixes prior to running trinotate:

```bash
for f in processing/02_annotations/*.cds.fa; do
    perl scripts/addPrefix.pl $f $(basename $f .cds.fa) > processing/02_annotations/$(basename $f .fa).prefix.fa
done
```

Make the obligatory `gene_trans_map`:

```bash
for f in processing/02_annotations/*.prefix.fa; do
    grep '^>' $f | awk '{ print $1 }' | sed -e 's/>\([^[:space:]]\+\)$/\1\t\1/g'
done > gene_trans_map.txt
grep '^>' ../00_transdecoder/dd_Smed_v6/dd_Smed_v6.pcf.contigs.fasta.fnn | awk '{ print $1 }' | sed -e 's/>\([^[:space:]]\+\)$/\1\t\1/g; s/_[[:digit:]]\([[:space:]]\)/\1/' > gene_trans_map.smed.txt
```

Now we create the transdecoder sequence files:

```bash
for s in aq ml nv ta ; do
    python scripts/convert_cds_to_transdecoder.py --out-dir processing/07_trinotate processing/07_trinotate/$s.gene_trans_map.txt processing/02_annotations/$s.cds.prefix.fa processing/06_oma/DB/$s.fa
done
python scripts/convert_cds_to_transdecoder.py --out-dir processing/07_trinotate processing/07_trinotate/gene_trans_map.smed.txt processing/00_transdecoder/dd_Smed_v6/dd_Smed_v6.pcf.contigs.fasta.fnn processing/00_transdecoder/dd_Smed_v6/dd_Smed_v6.pcf.contigs.fasta.faa
```

Copy initial sources over from the trinotate pipeline:

```bash
TRINOTATE_SRC=#CHANGEME#
cp $TRINOTATE_SRC/db/uniprot_sprot.pep  $TRINOTATE_SRC/db/uniprot_sprot.pep.phr  $TRINOTATE_SRC/db/uniprot_sprot.pep.pin  $TRINOTATE_SRC/db/uniprot_sprot.pep.psq $TRINOTATE_SRC/db/Trinotate.sqlite Trinotate/
```

Now we begin the blast and hmmer jobs (note this should probably be split and database locations may be different):

```bash
HMMER_SRC=#CHANGEME#
TRINOTATE_LOC=#CHANGEME#
cat ??.cds.prefix.tdfmt.fa > allcds.fa
cat ??.tdfmt.fa > allprot.fa
hmmscan --domtblout allprot.domtbl $HMMER_SRC/pfam_a.hmm allprot.fa
blastp -query allprot.fa -db $TRINOTATE_LOC/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 > allprot.txt
blastx -query allcds.fa -db $TRINOTATE_LOC/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 > allprot.txt
```

Now we take a pre-created boilerplate Trinotate.sqlite file and add our sequences to them (these were repeated for S. medierranea only into `Trinotate/Trinotate.smed.sqlite`):

```bash
Trinotate Trinotate/Trinotate.sqlite init --gene_trans_map gene_trans_map.txt --transcript_fasta allcds.fa --transdecoder_pep allprot.fa
```

Concatenate the results and load them into the database, spit out the report:

```bash
cat splitprot/*.txt > blastp.out.txt
cat splitcds/*.txt > blastx.out.txt
cat splitprot/*.domtbl > hmmer.out.domtbl
Trinotate Trinotate/Trinotate.sqlite LOAD_swissprot_blastp blastp.out.txt
Trinotate Trinotate/Trinotate.sqlite LOAD_swissprot_blastx blastx.out.txt
Trinotate Trinotate/Trinotate.sqlite LOAD_pfam hmmer.out.domtbl
Trinotate Trinotate/Trinotate.sqlite report > trinotate_report.txt
```


