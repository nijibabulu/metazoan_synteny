
# Data Acquisition

This describes the steps required to get all the data files in shape for the synteny pipeline and single cell analysis in the paper Zimmermann et al. *Ancient animal genome architecture reflects cell type identities*.

## Proteomes

**Assemble and synchronize IDs of protein sequences.**

The genome builds for the the NEE genomes were listed on the linked samples for the GEO submission GSE111068 (referenced in the NEE paper doi: 10.1038/s41559-018-0575-6) downloaded into `data`.

- Amphimedon queenslandica (Aqu1): <https://metazoa.ensembl.org/Amphimedon_queenslandica/Info/Index>
- Mnemiopsis leidyi (ML2): <https://research.nhgri.nih.gov/mnemiopsis/download/download.cgi?dl=genome>
- Trichoplax adhaerens (Triad1): <https://genome.jgi.doe.gov/Triad1/Triad1.home.html>
- Schmidtea mediterranea (`dd_Smed_v6`): <http://planmine.mpi-cbg.de/planmine/report.do?id=2000001#ad-image-0>

The latter has repeired IDs as follows:

```
sed -e 's/^>.*dd_Smes/>dd_Smes/; s/,.*$//' GCA_002600895.1_ASM260089v1_genomic.fna > GCA_002600895.1_ASM260089v1_genomic.dd_Smes_ids.fna
```

Nematostella (doi: 10.1016/j.cell.2018.05.019) was a combination of ENSEMBL metazoa genes (`^EMNVEG.*`), JGI genes downloaded from ENSEMBL Metazoa release 39 and NVE genes as curated described in the paper. The ENSEMBL genes consisted entirely of the metazoan signal recognition particle RNA, one case of Nuclear RNAse P and deprecated identifiers. They were ignored.

Determine which of the two IDs included in his merged dual id bed file were included as a table column (rule of thumb is the JGI ID but this was checked anyway). If the ID was not in either table, we used the NVE id:

```bash
python scripts/decide_nv_gene_name.py \
  data/Nematostella_merged_annotation_DUAL_NAMES.bed \
  data/GSM2942625_MARSseq_UMI_Table_Adult_Nematostella.txt \
  data/GSM2942626_MARSseq_UMI_Table_Larva_Nematostella.txt > \
    processing/02_annotations/nv.bed
```

Proteome sequences were gathered between the two proteomes with `scripts/assemble_nematostella_proteome.py`.

```bash
python scripts/assemble_nematostella_proteome.py \
  data/Nematostella_vectensis.ASM20922v1.pep.all.fa \
  data/nveGenes.good.130208.longCDS.protein.corrected_inframe_stops.fa  \
  processing/02_annotations/nv.bed > \
    processing/01_dbs/nv.faa
```

Additionally, in order to run Trinotate, we needed CDS sequences. This is provided with the other proteomes, but the NVE sequences include UTRs. As we're not concerned with UTRs, I simply extract the CDS sequences:

```bash
python scripts/extract_transcript_sequence.py --cds data/nemVec1.fa data/nveGenes.good.130208.longCDS.bed  > processing/02_annotations/nveGenes.good.130208.longCDS.nt.fa
python scripts/assemble_nematostella_proteome.py \
  data/Nematostella_vectensis.ASM20922v1.cds.all.fa \
  processing/02_annotations/nveGenes.good.130208.longCDS.nt.fa \
  processing/02_annotations/nv.bed > \
    processing/02_annotations/nv.cds.fa
```

Convert the Trichoplax ids:

```bash
sed -e 's/^>jgi|Triad1|\([[:digit:]]\+\).*$/>Tadh_P\1/' data/Triad1_best_proteins.fasta > processing/01_dbs/ta.faa
sed -e 's/^>jgi|Triad1|\([[:digit:]]\+\).*$/>Tadh_P\1/' data/Triad1_best_transcripts.fasta > processing/02_annotations/ta.cds.fa
```

The remaining sequences were copied over:
```bash
cp data/Amphimedon_queenslandica.Aqu1.pep.all.fa processing/01_dbs/aq.faa
cp data/ML2.2.aa processing/01_dbs/ml.faa
```

Schmidtea mediterranea assembly sequences were converted to protein coding sequences via GeneMarkS-T 5.1 (doi `10.1093/nar/gkv227`):

```bash
module load genemarks-t
cd processing/00_transdecoder/dd_Smed_v6
gmst.pl --verbose --fna --faa ../../../data/dd_Smed_v6.pcf.contigs.fasta
```

Only sequences found in the count matrix were used for orthology analysis:

```bash
python scripts/filter_schmidtea_genes.py data/schmidtea_dge.txt processing/00_transdecoder/dd_Smed_v6/dd_Smed_v6.pcf.contigs.fasta.faa > processing/00_transdecoder/dd_Smed_v6/dd_Smed_v6.pcf.contigs.fasta.filtered.faa
```

The original sequence file had 41745 transcripts, the proteins totaled only 26439 and the count matrix had 28065. The intersection of the latter two is 25930.

Once orthoblocks are detected, we will use the values in the expression matrix without the files (we will have to remove them before applying the analyses in SCS).

## Chromosome Files

**Reformat chromosomal annotations for the synteny pipeline**

Convert the Trichoplax IDs to the ones used in the supplemental of the paper:

```bash
sed -e 's/proteinId \([^;]\+\);/proteinId "Tadh_P\1";/' data/Triad1_best_genes.gff > processing/02_annotations/Triad1_best_genes.Tadh_P.ids.gff
```

In order to make files compatible with the microsynteny pipeline pipeline, run the following (note the pipeline depends on the filenames ending in `.chrom`):

```bash
perl scripts/makeMap.pl data/ML2.2.gff3 ml  mRNA ID > processing/02_annotations/ml.chrom # 16548 genes found!
perl scripts/makeMap.pl processing/02_annotations/Triad1_best_genes.Tadh_P.ids.gff ta CDS proteinId > processing/02_annotations/ta.chrom #  11520 genes found!
perl scripts/makeMap.pl data/Amphimedon_queenslandica.Aqu1.39.gff3 aq mRNA transcript_id > processing/02_annotations/aq.chrom # 43615 genes found!
perl scripts/makeMapFromBed.pl processing/02_annotations/nv.bed nv > processing/02_annotations/nv.chrom #  32268 genes found! # 24774 genes found!
perl scripts/makeMapFromBed.pl data/dd_Smed_v6.bed sm > processing/02_annotations/sm.chrom # 23347 genes found!
```

`data/Triad1_best_genes.gff` is a gff version 2 file which did not have meta (mRNA) features. This presumably works fine, and I checked the uniqueness of all names associated with `start_codon` and `stop_codon` features (some genes are partial).

## Gene Names

**Automatically add speaking names to the genes**

### Emapper

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

### Trinotate 

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

## Extract metacells from single cell data

The identity of the metacells are encoded in figures 1c, 1f, 2c and 2f [of the paper](https://dx.doi.org/10.1038/s41559-018-0575-6). In the cases where there are multiple metacells per cell identity, the ranges are not listed in numerical order but instead in the order of the sheet in the corresponding excel workbook (and probably in some internal ordering). In order to extract a map of metacells to cell identies in these cases, I wrote a script to parse the excel sheets and, using hardcoded values based on the figure, create text files:

```bash
python scripts/extract_metacells.py data/SpongeEtcSuppl/41559_2018_575_MOESM5_ESM.xlsx \
    data/SpongeEtcSuppl/41559_2018_575_MOESM6_ESM.xlsx \
    data/SpongeEtcSuppl/41559_2018_575_MOESM7_ESM.xlsx \
    data/SpongeEtcSuppl/41559_2018_575_MOESM8_ESM.xlsx \
    processing/08_metacells/aq_adult.metacells.txt \
    processing/08_metacells/aq_larva.metacells.txt \
    processing/08_metacells/ml.metacells.txt \
    processing/08_metacells/ta.metacells.txt
```

Also, extract marker genes from the supplementary data:

```bash
python scripts/extract_markers.py data/NematostellaSuppl/mmc2.xlsx > processing/10_markers/nv.markers.txt
python scripts/extract_markers.py data/SpongeEtcSuppl/41559_2018_575_MOESM6_ESM.xlsx > processing/10_markers/aq.adult.markers.txt
python scripts/extract_markers.py data/SpongeEtcSuppl/41559_2018_575_MOESM7_ESM.xlsx > processing/10_markers/ml.markers.txt
python scripts/extract_markers.py data/SpongeEtcSuppl/41559_2018_575_MOESM8_ESM.xlsx > processing/10_markers/ta.markers.txt
```

The data for Nematostella taken from the download at (`http://www.wisdom.weizmann.ac.il/~/arnau/Single_cell_datasets/`), and processed as such:

```bash
tail -n +2 data/sebe_pedros_nv_data_2018/Annotation_and_config_files/Adult_metacell_annotation > processing/08_metacells/nv_adult_metacells.txt
tail -n +2 data/sebe_pedros_nv_data_2018/Annotation_and_config_files/Larva_metacell_annotation > processing/08_metacells/nv_larva_metacells.txt
```

This creates files we can then use for the remainder of the synteny pipeline

```bash
for n in $(seq 0 999); do
    mkdir -p blocks/$n
     cd blocks/$n
     perl ../../../../MICROSYNT/prepMicroSynt.pl $(ls ../../../02_annotations/*.chrom | xargs | sed -e 's/ /,/g') 5 ../../random_orthogroups/Orthogroups.prefixed.txt.rand.$n
     cd -
done
```

