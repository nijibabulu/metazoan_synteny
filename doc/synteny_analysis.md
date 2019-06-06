# Synteny Analysis

The synteny analysis is done in 3 basic steps:

1. [Determine orthogroups](#orthofinder) - This is done by taking gene groups after the MCL step of [orthofinder](https://github.com/davidemms/OrthoFinder).
2. [Create clusters of blocks](#orthoblocks) - With the method described in [Simakov et al 2013](https://doi.org/10.1038/nature11696).
3. [Annotate genes in blocks (Optional)](#annotate-synteny-clusters) - Using the gene names found in the [gene_functional_annotation.md](gene_functional_annotation.md) steps. Gene annotation is optional, however it is necessary to reformat the synteny blocks as detailed in this step.

## OrthoFinder

All genes were preprocessed to include their species name:

```bash
sed -ie 's/^>/>aq_/' processing/01_dbs/aq.faa
sed -ie 's/^>/>ml_/' processing/01_dbs/ml.faa
sed -ie 's/^>/>nv_/' processing/01_dbs/nv.faa
sed -ie 's/^>/>ta_/' processing/01_dbs/ta.faa
sed -e 's/^>/>sm_/' processing/00_transdecoder/dd_Smed_v6/dd_Smed_v6.pcf.contigs.fasta.filtered.faa > processing/01_dbs/sm.faa
```

First the blast commands were generated

```
orthofinder -f processing/01_dbs/ -op
```

The remainder of the orthofinder algorithm was run as follows:

```bash
orthofinder -a 4 -t 4 -b processing/01_dbs/Results_Nov01/WorkingDirectory
```

Convert the `Orthogroups.txt` file to the microsynteny pipeline format:

```bash
perl scripts/orthoFinderToOrthogroup.pl processing/03_orthofinder/Results_Nov01/WorkingDirectory/Orthogroups.txt > processing/03_orthofinder/orthogroups.txt
```


## Orthoblocks

Run the preparation of the blocks (note this creates scripts for submission to a slurm-based queueing engine.)

```bash
cd processing/04_blocks
perl ../../scripts/prepMicroSynt.pl $(ls ../02_annotations/??.chrom | xargs | sed -e 's/ /,/g') 5 ../03_orthofinder/orthogroups.txt
```

Generate clusters of blocks:

```bash
cd processing/04_blocks
perl ../../scripts/makeClusters3.pl $(ls ../02_annotations/??.chrom | xargs | sed -e 's/ /,/g') .5.blocks 3 0.3 0.5 > nmax5.clust
```

## Annotate Synteny Clusters

To interpret the clustered synteny blocks and group them by connectivity, as well as add annotation data to a final column you may reformat these files with the following script:

```bash
python scripts/reformat_synteny_blocks.py \
    processing/04_blocks/nmax5.clust \
    processing/05_emapper/emapper.annotations.wsmed.txt \
    processing/07_trinotate/trinotate_report.wsmed.txt > processing/04_blocks/nmax5.clust.refmt
```
For planaria, at this stage, we don't care about the transcript id, since they are not included in the expression matrix. Therefore we remove them manually:

```bash
sed -e 's/\(sm_dd_[^,\s]*\)_[[:digit:]]/\1/g' processing/04_blocks/nmax5.clust > processing/04_blocks/nmax5.noplntxids.clust
sed -e 's/\(sm_dd_[^,\s]*\)_[[:digit:]]/\1/g' processing/04_blocks/nmax5.clust.refmt > processing/04_blocks/nmax5.clust.noplntxids.refmt
sed -e 's/\(sm_dd_[^,\s]*\)_[[:digit:]]/\1/g' processing/03_orthofinder/orthogroups.wschmidtea.txt > processing/03_orthofinder/orthogroups.wschmidtea.noplntxids.txt
sed -e 's/\(sm_dd_[^,\s]*\)_[[:digit:]]/\1/g' processing/03_orthofinder/orthogroups.expanded.geneids.txt > processing/03_orthofinder/orthogroups.expanded.geneids.noplntxids.txt
```

Similarly, we have to remove those from the chrom files for the expression analysis (not synteny):

```bash
sed -e 's/\(sm_dd_[^,\s]*\)_[[:digit:]]/\1/g' processing/02_annotations/sm.chrom > processing/02_annotations/sm.geneid.chrom
```

## Paralog content of Microsynteny blocks

In order to assess the fraction of paralogs in each of the microsynteny blocks, we may use the output of orthofinder implicitly via the orthologs output of OrthoFinder. A script is provided which summarizes this information which requires the original `Orthogroups.txt` file from OrthoFinder, the `Orthologues` output directory, and the reformatted blocks file

```
python3 scripts/blockwise_paralogous_content.py processing/03_orthofinder/Results_Nov01/WorkingDirectory/Orthogroups.txt processing/03_orthofinder/Results_Nov01/WorkingDirectory/Orthologues_Nov01/Orthologues/ processing/04_blocks/nmax5.clust.noplntxids.refmt > processing/04_blocks/paralogous_fractions.txt
```

We can visualize this as in Supplementary Figure 1a, as well as the number of blocks in each species, using the script in `SCS/scripts/plot_paralogous.R`. 
