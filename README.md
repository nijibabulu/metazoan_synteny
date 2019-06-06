# Metazoan Synteny Expression Analysis

Bob Zimmermann and Oleg Simakov

This repository contains the code used in Zimmermann et al. **Ancient animal genome architecture reflects cell type identities**, which contains elements from Simkov, et al. **Insights into bilaterian evolution from three spiralian genomes**, *Nature* 2013.  Aspects of these scripts are being developed elsewhere in more digestable forms, but the code is provided here in order to gain insight into the analysis beyond the text.

The code is divided up into two major parts: 

* data processing and synteny analysis scripts are contained here in the `scripts/` directory, and
* an R package and accompanying Rscripts which can be run from the command line are in the `SCS/` directory.

A detailed, step-by-step protocol can be found in [doc/pipeline.md](doc/pipeline.md). 

Files to reproduce the results are available at [here](http://fileshare.csb.univie.ac.at/metazoan_synteny_data/extdata.tar.gz). These should be placed and untarred into the `SCS/inst` directory. (See docs for details.)

If any code appears missing or buggy feel free to raise an issue or contact me.
