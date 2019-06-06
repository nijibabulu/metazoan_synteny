#! /usr/bin/env Rscript
library(tidyverse)
library(ggthemes)
library(cowplot)
library(devtools)
load_all('.')

main <- function() {
  synBlockESSs.norm <- load_data_set('synBlockESs.norm')
  genomes <- load_data_set('genomes')

  circlePlot(synBlockESs.norm, genomes, block_names=c('37(aq)', '350(ta)', '38(nv)'),
             normalize="scale", show.colorbar=T,
             output='rebfigs/Fig2B.cur.colorbar.pdf')
  circlePlot(synBlockESs.norm,  genomes, block_names=c('71(aq)', '498(ta)', '70(nv)'),
             normalize="scale", show.colorbar=F,
             orientations=c(1,-1,1), output='rebfigs/Fig2C.cur.pdf')

  clusters.reducedA <- load_data_set('clusters.reducedA')
  query(clusters.reducedA, c('aq','ta','nv')) %>%
    getClusters() %>%
    map(~circlePlot(synBlockESs.norm, genomes, cluster=., species_order=c('aq','ta','nv'),
                    normalize="scale", show.colorbar=F,
                    output=file.path('rebfigs/circles', paste0(name(.), '.pdf'))))

  query(clusters.reducedA, c('aq','ml','ta','nv')) %>% getClusters() %>%
    map(~circlePlot(synBlockESs.norm, genomes, cluster=., species_order=c('aq','ta','ml','nv'),
                    normalize="scale", show.colorbar=F,
                    output=file.path('rebfigs/circles', paste0(name(.), '.pdf')),
                    out.width=15))

    query(clusters.reducedA, c('aq','ta','nv','sm')) %>% getClusters() %>%
    map(~circlePlot(synBlockESs.norm, genomes, cluster=., species_order=c('aq','ta','nv','sm'),
                    normalize="scale", show.colorbar=F,
                    output=file.path('rebfigs/circles', paste0(name(.), '.pdf')),
                    out.width=15, out.height=10))
}

if(sys.nframe() == 0) {
  main()
}
