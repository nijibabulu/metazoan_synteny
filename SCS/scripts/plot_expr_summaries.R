#! /usr/bin/env Rscript
library(tidyverse)
library(ggthemes)
library(cowplot)
library(devtools)
load_all('.')

main <- function() {
  synBlockESs.norm <- load_data_set('synBlockESs.norm')
  genomes <- load_data_set('genomes')

  
  # this plot can be used exclusively for extracting the colorbar
  summaryGridPlot(synBlockESs.norm, genomes, block_names=c('37(aq)', '350(ta)', '38(nv)'),
             normalize="scale", show.colorbar=T,
             output='Fig2.colorbar.pdf')
  summaryGridPlot(synBlockESs.norm, genomes, block_names=c('37(aq)', '350(ta)', '38(nv)'),
             normalize="sum", show.colorbar=F,
             output='Fig2B.sgp.sum.cur.pdf', return_data_frames=T)
  summaryGridPlot(synBlockESs.norm,  genomes, block_names=c('71(aq)', '498(ta)', '70(nv)'),
             normalize="none", show.colorbar=F,
             orientations=c(1,-1,1), output='Fig2C.sgp.norm.cur.pdf', return_data_frames=T)
  
  
  # examples of more standard, seurat-style dot plots
  circlePlot(synBlockESs.norm, genomes, block_names=c('37(aq)', '350(ta)', '38(nv)'),
             normalize="scale", show.colorbar=F)
  circlePlot(synBlockESs.norm,  genomes, block_names=c('71(aq)', '498(ta)', '70(nv)'),
             normalize="scale", show.colorbar=F)
 
  # plot all of the synteny blocks shared among aq, ta and nv
  clusters.reducedA <- load_data_set('clusters.reducedA')
  query(clusters.reducedA, c('aq','ta','nv')) %>%
    getClusters() %>%
    head(n=1) %>%
    map(~circlePlot(synBlockESs.norm, genomes, cluster=., species_order=c('aq','ta','nv'),
                    normalize="scale", show.colorbar=F))
                    #output=file.path('circle', paste0(name(.), '.pdf'))))
  
  # examples of "filled box" circle plots
  circlePlot(synBlockESs.norm, genomes, block_names=c('37(aq)', '350(ta)', '38(nv)'),
             normalize="scale", show.colorbar=F, fill.box=T)
  circlePlot(synBlockESs.norm,  genomes, block_names=c('71(aq)', '498(ta)', '70(nv)'),
             normalize="scale", show.colorbar=F, fill.box=T)
 
}

if(sys.nframe() == 0) {
  main()
}
