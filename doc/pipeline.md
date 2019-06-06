# Metazoan Synteny Expression Analyses

Aspects of these analyses are being developed into tools elsewhere, but the code is provided here in order to gain insight into the analysis beyond the text. 

This document contains the main steps of the analysis, whereas data acquisition, gene functional annotation and noncoding elements analysis are relegated to separate files. If you wish to jump in, you may start at the "[(Optional) Download Processed Data](#optional-download-processed-data)" to see how the main part of the analysis works.

1. **Data Acquisition** - Detailed in [data_acquisition.md](data_acquisition.md).
1. **Gene Functional Identification** - Annotate the gene names of the proteomes, detailed in [gene_functional_annotation.md](gene_functional_annotation.md)
1. **Synteny Analysis** - Finding conserved microsynteny blocks among the genomes, detailed in [synteny_analysis.md](synteny_analysis.md).
1. [(Optional) Download Processed Data](#optional-download-processed-data) - Download the processed data from the above steps.
1. [Installing R Package Prerequisites](#installing-r-package-prerequisites) - Everything you need to run the R code here.
1. [Load Data](#load-data) - Import the synteny blocks and single cell data into `.Rda` files for use in the analysis pipeline.
1. [Randomize Blocks](#randomize-blocks) - Shuffle and sample blocks from the genome for analysis
1. [Imputation and Correlation Analysis](#imputation-and-correlation-analysis) - All of the figures from the manuscript.
1. **Noncoding elements** - Search for enriched motifs in the synteny blocks compared to random blocks, detailed in [noncoding_elements.md](noncoding_elements.md)

# (Optional) Download Processed Data

A tarball of the necessary files from the above processing that are needed to proceed are available [here](http://fileshare.csb.univie.ac.at/metazoan_synteny_data/extdata.tar.gz). To use it, download it and untar the file into the directory `SCS/inst` directory. If it worked, you should be able to see the file

`SCS/inst/extdata/aq_adult.metacells.txt`

among 21 others.

# Installing R Package Prerequisites

The simplest workflow for this package is to use [devtools](https://github.com/r-lib/devtools), working in the directory of the where the `DESCRIPTION` file is (`SCS/`). In order to obtain the dependencies, you will have to do the following:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges","GenomeInfoDb","WGCNA"))

if (!require("devtools", quietly = TRUE)) 
    install.packages("devtools")

library(devtools)
install_deps()
```

If you wish to do imputation, the following need to be installed (as desired):

```R
install.packages("Rmagic") 
library(devtools)
install_github("gongx030/DrImpute") 
install_github("Vivianstats/scImpute")
```

## Load Data

Most of the following steps can be run from the command line, but since the scripts can take long to run, it is useful to run interactively to be able to trace back and fix issues as they come. 

The `load_data.R` and `preprocess_data.R` scripts loads the data sets and places them in persistent `.Rda` files under `SCS/inst/data_sets`.

You may visualize the intersections as in Supplementary Figure 1B with the script `SCS/scripts/plot_cluster_overlaps.R`.

**Note:** Clusters from the `clusters` object are the clusters used in Supplementary Figure 1. The total count of 270 in the main text refers to the total number of clusters (277) minus the number of clusters which consist only of genes from a single genome. (4 from Nematostella vectensis, 3 from Amphimedon queenslandica).

## Randomization of Blocks

Ransomization takes place in three different ways in our manuscript, by shuffling the gene content among the blocks within a genome ("Shuffled"), by sampling random blocks of genes within the genome, either unconstrained ("Sampled Blocks") or constrained to contain only genes with orthologs in the study set ("Sampled Blocks of Orthologs"). 

These are performed in the script `SCS/scripts/randomize_blocks.R`, in turn stores the data persistently in `.Rda` filesin `SCS/inst/data_sets`.

This script also produces distributions of intervening sizes between syntenic genes within the synteny blocks. To visualize the distribution of these, as well as their fits to the negative binomial distribution a la Supplementary Figure 1C, you may use the `SCS/scripts/plot_fits.R` script.

# Imputation and Correlation Analysis

The various count imputations are driven by the script `SCS/scripts/impute.R`, which can be run via `Rscript` and has a command line interface. These calculations are rather lengthy as they include computing the gene-wise expression correlation. It is necessary to generate at least the correlation for the unimputed counts via this script by specifying the method `none`:

```
Rscript scripts/impute.R --expression=ess.norm --method=none 0 unimputed.norm
```

This will save an `.Rda` file in the `inst/` directory. A bash script in `SCS/scripts/create_imputation_script.sh` generates a list of tasks that includes the settings in the Supplementary Figure 2. Note that, as in the paper, it was infeasible to do lareger imputation computation for *S. mediterranea* due to the large amount of cells.

Despite being precomputed, the files are very large and unruly to handle in parallel, so the relevant values are extracted and block correlations computed with with `SCS/scripts/compute_block_correlations.R`. An example script generating commands that were used for the paper can be found in `SCS/scripts/create_block_correlation_script.sh`.

The various plots were then created, as in Figure 1B and Supplementary Figure 2, in `SCS/scripts/plot_correlations.R`.

Tau and expression bias values were computed in the `SCS/scripts/compute_expression_bias.R`, and are plotted as in Figure 2a and Supplementary Figure 8 in `SCS/scripts/plot_expression_bias.R`.

## Expression Matrix Plots

Supplementary Figures 3-7 represent the full expression data for all synteny blocks blocks that were found. Plotting these is coupled together with the calculation of the background information. Examples (as is done in the paper) are shown in the script `SCS/scripts/plot_expr.R`. The plots will take longer to compute and render, and this can be reduced by changing the number of background samples to 0.

Expression matrices summarized by cell type as in Figure 2B-C is done in `SCS/scripts/plot_expr_summaries.R`. Examples of other plot types are also shown in this script.

## MAGIC t-disparity plots

In order to see if the imputation will resolve the issue of drop out observation, we use magic among others. But MAGIC seems to over-inflate the counts with *A. queenslandica*. In order to demonstrate that the automatic selection of the t-disparity parameter is faulty,  we plot these. The plot information is not exposed in R so we ran this in python in the `imputation` directory:

```
python do_imputation.py ../SCS/inst/extdata/GSM3021561_Amphimedon_queenslandica_adult_UMI_table.txt aq.tplot.png
python do_imputation.py ../SCS/inst/extdata/GSM3021563_Mnemiopsis_leidyi_UMI_table.txt ml.tplot.png
python do_imputation.py ../SCS/inst/extdata/GSM3021564_Trichoplax_adhaerens_UMI_table.txt ta.tplot.png
python do_imputation.py ../SCS/inst/extdata/GSM2942625_MARSseq_UMI_Table_Adult_Nematostella.txt nv.tplot.png
python do_imputation.py ../SCS/inst/extdata/schmidtea_dge.txt sm.tplot.png
```


