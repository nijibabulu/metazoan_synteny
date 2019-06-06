# Noncoding elements

The below is mostly a description of how the [AME](http://meme-suite.org/doc/ame.html) analysis was done in detail. Doing this requires fully several hundred hours of CPU time.

Construct masked versions of the genome containing no coding exons. Note these are not gene beds, but based on CDS exons. Note the nematostella transcriptome from Arnau had no exon annotation so I just masked both JGI and NVE annotations.

```
grep CDS Amphimedon_queenslandica.Aqu1.39.gff3 | convert2bed -i gff -o bed > Amphimedon_queenslandica.Aqu1.39.cds_exon.bed
grep CDS ML2.2.gff3 | convert2bed -i gff -o bed > ML2.2.cds_exon.bed
python ../scripts/genes_bed_to_cds_exons.py dd_Smed_v6.bed > dd_Smed_v6.cds_exon.bed
(python ../scripts/genes_bed_to_cds_exons.py nveGenes.good.130208.longCDS.bed ; grep CDS Nematostella_vectensis.ASM20922v1.39.gff3 | convert2bed -i gff -o bed | awk -vOFS='\t' '{print $1,$2,$3}') > nv_jgi_nve_merge.cds_mask.bed
grep CDS Triad1_best_genes.gff | convert2bed -i gff -o bed > Triad1_best_genes.cds_exon.bed

bedtools maskfasta -bed Amphimedon_queenslandica.Aqu1.39.cds_exon.bed -fi  Amphimedon_queenslandica.Aqu1.dna_rm.toplevel.fa  -fo Amphimedon_queenslandica.Aqu1.dna_rm.cds_mask.toplevel.fa
# bedtools maskfasta -bed ML2.2.cds_exon.bed -fi MlScaffold09.nt -fo MlScaffold09.cds_mask.nt
bedtools maskfasta -bed ML2.2.cds_exon.bed -fi Mnemiopsis_leidyi.MneLei_Aug2011.dna_rm.toplevel.fa -fo Mnemiopsis_leidyi.MneLei_Aug2011.dna_rm.toplevel.cds_mask.fa
# bedtools maskfasta -bed dd_Smed_v6.cds_exon.bed -fi dd_Smes_g4.fasta -fo dd_Smes_g4.cds_mask.fasta 
bedtools maskfasta -bed dd_Smed_v6.cds_exon.bed -fi GCA_002600895.1_ASM260089v1_genomic.dd_Smes_ids.fna -fo GCA_002600895.1_ASM260089v1_genomic.dd_Smes_ids.cds_mask.fna
bedtools maskfasta -bed nv_jgi_nve_merge.cds_mask.bed -fi nemVec1.RM.fa -fo nemVec1.RM.cds_mask.fa
# bedtools maskfasta -bed Triad1_best_genes.cds_exon.bed -fi Triad1_masked_genomic_scaffolds.fasta -fo Triad1_masked_genomic_scaffolds.cds_mask.fasta
bedtools maskfasta -bed Triad1_best_genes.cds_exon.bed -fi Trichoplax_adhaerens.ASM15027v1.dna_rm.toplevel.fa -fo Trichoplax_adhaerens.ASM15027v1.dna_rm.toplevel.cds_mask.fa
```

Additionally, hard mask genomes for which there is only lowercase masking:

```
lcmask dd_Smes_g4.cds_mask.fasta > dd_Smes_g4.cds_mask.hardmask.fasta
```

I then obtained the block locations from the Rdata file via the script `SCS/scripts/write_blocks_bed.R` and moved them to the `data/` directory.

These were then used to extract the sequences

```
bedtools getfasta -fi Amphimedon_queenslandica.Aqu1.dna_rm.cds_mask.toplevel.fa -bed aq.observed.bed  -fo aq.observed.fa -name
bedtools getfasta -fi Amphimedon_queenslandica.Aqu1.dna_rm.cds_mask.toplevel.fa -bed aq.sampled.ortholog.bed  -fo aq.sampled.ortholog.fa -name
# bedtools getfasta -fi MlScaffold09.cds_mask.nt -bed ml.observed.bed -fo ml.observed.fa -name
# bedtools getfasta -fi MlScaffold09.cds_mask.nt -bed ml.sampled.ortholog.bed -fo ml.sampled.ortholog.fa -name
bedtools getfasta -fi Mnemiopsis_leidyi.MneLei_Aug2011.dna_rm.toplevel.cds_mask.fa -bed ml.observed.bed -fo ml.observed.fa -name
bedtools getfasta -fi Mnemiopsis_leidyi.MneLei_Aug2011.dna_rm.toplevel.cds_mask.fa -bed ml.sampled.ortholog.bed -fo ml.sampled.ortholog.fa -name
# bedtools getfasta -fi dd_Smes_g4.cds_mask.fasta -bed sm.observed.bed -fo sm.observed.fa -name
# bedtools getfasta -fi dd_Smes_g4.cds_mask.fasta -bed sm.sampled.ortholog.bed -fo sm.sampled.ortholog.fa -name
bedtools getfasta -fi GCA_002600895.1_ASM260089v1_genomic.dd_Smes_ids.cds_mask.fna -bed sm.observed.bed -fo sm.observed.fa -name
bedtools getfasta -fi GCA_002600895.1_ASM260089v1_genomic.dd_Smes_ids.cds_mask.fna -bed sm.sampled.ortholog.bed -fo sm.sampled.ortholog.fa -name
bedtools getfasta -fi Trichoplax_adhaerens.ASM15027v1.dna_rm.toplevel.cds_mask.fa -bed ta.observed.bed -fo ta.observed.fa -name
bedtools getfasta -fi Trichoplax_adhaerens.ASM15027v1.dna_rm.toplevel.cds_mask.fa -bed ta.sampled.ortholog.bed -fo ta.sampled.ortholog.fa -name
bedtools getfasta -fi nemVec1.RM.cds_mask.fa -bed nv.observed.bed -fo nv.observed.fa -name
bedtools getfasta -fi nemVec1.RM.cds_mask.fa -bed nv.sampled.ortholog.bed -fo nv.sampled.ortholog.fa -name

bedtools getfasta -fi Amphimedon_queenslandica.Aqu1.dna_rm.cds_mask.toplevel.fa -bed aq.sampled.all.bed  -fo aq.sampled.all.fa -name
bedtools getfasta -fi nemVec1.RM.cds_mask.fa -bed nv.sampled.all.bed -fo nv.sampled.all.fa -name
bedtools getfasta -fi Trichoplax_adhaerens.ASM15027v1.dna_rm.toplevel.cds_mask.fa -bed ta.sampled.all.bed -fo ta.sampled.all.fa -name
bedtools getfasta -fi dd_Smes_g4.fasta -bed sm.sampled.all.bed -fo sm.sampled.all.fa -name
bedtools getfasta -fi MlScaffold09.cds_mask.nt -bed ml.sampled.all.bed -fo ml.sampled.all.fa -name
```


I then take all the analagous 100 background sequences for each of the synteny blocks as background:

```
mkdir -p {aq,ml,nv,ta,sm}_bg
for sp in aq ml ta nv sm; do 
    python3 ../../scripts/get_background_files.py ../../data/$sp.sampled.ortholog.fa ${sp}_bg
done
```

I then ran `AME` to detect motifs on individual blocks against the entire sampled ortholog blocks as background:

```
AME=$(realpath ../../meme-5.0.4/stage/bin/ame)
JASPAR=$(realpath ../../meme-5.0.4/stage/share/meme-5.0.4/db/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme)
for sp in aq ml ta nv sm; do 
    for f in ${sp}_split/*.fa; do 
        id=$(basename $f .fa)
        bg=../../data/$sp.sampled.ortholog.fa
        echo $AME --o \"${sp}_ame/${id}\" --control \"$bg\" \"$f\" $JASPAR
    done
done
```

Next I ran the group of all blocks, again using the sampled ortholog blocks as background:

```
AME=$(realpath ../../meme-5.0.4/stage/bin/ame)
JASPAR=$(realpath ../../meme-5.0.4/stage/share/meme-5.0.4/db/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme)
for sp in ml ta  sm ; do 
    echo $AME --o ${sp}_ame_all --control ../../data/$sp.sampled.ortholog.fa ../../data/$sp.observed.fa $JASPAR 
done
```

As a control, run on as many sampled blocks from the genome as we found in the observed set. First generate the samples (in data directory):

```
# initial single control:
for sp in aq ml ta  nv sm ; do 
    head -n $(wc -l $sp.observed.fa | awk '{print $1}') $sp.sampled.all.fa > $sp.sampled.all.1.fa
done
```

Now run `ame`:

```
AME=$(realpath ../../meme-5.0.4/stage/bin/ame)
JASPAR=$(realpath ../../meme-5.0.4/stage/share/meme-5.0.4/db/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme)
for sp in aq ml ta  nv sm ; do 
    for f in split_samples/$sp/*.fa; do
        echo $AME --o ${sp}_ame_ctl/$(basename $f .fa) --control ../../data/$sp.sampled.ortholog.fa $f $JASPAR 
    done
done
```

Get a summary:

```
for sp in aq ml ta nv sm; do 
    echo "$sp $(grep 'tp$' ${sp}_ame_all/sequences.tsv | wc -l) $(grep 'tp$' ${sp}_ame_ctl/*/sequences.tsv | wc -l)"; 
done > ame_all_summary.txt

for sp in aq ml ta nv sm; do 
    echo "fg $sp $( tail -n +2 ${sp}_ame_all/sequences.tsv | grep -v '^#' |  awk '{print $2}' | sort -u | wc -l)"
    for f in ${sp}_ame_ctl/$sp*/sequences.tsv; do 
        echo "bg $sp $(tail -n +2 $f | grep -v '^#' | awk '{print $2}' | sort -u | wc -l)"
    done
done > motifcounts_summary.txt
```

Finally, motif analysis was done with the script `SCS/scripts/motif_analysis.R`
