#! /usr/bin/env Rscript
library(tidyverse)
library(ggthemes)
library(cowplot)
library(devtools)
load_all('.')

compute_tau <- function(es, gene_wise=T) {
  xi <- if(gene_wise) {
    cellinfo(es) %>%
      crossing(gene=rownames(es)) %>%
      group_by(cellstate,gene) %>%
      summarize(xig=sum(mat(es)[.data$gene,.data$cell,drop=F])/
                    sum(mat(es)[.data$gene,,drop=F])) %>%
      summarize(xi=mean(xig))
  } else {
    cellinfo(es) %>%
      group_by(cellstate) %>%
      summarize(xi=sum(mat(es[,.data$cell,drop=F]))/n())
  }
  xi %>%
    mutate(xi_hat=xi/max(xi)) %>%
    summarize(tau=sum(1-xi_hat)/(n()-1)) %>%
    flatten_dbl()
}

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

compute_expression_bias <- function(es) {
  cellinfo(es) %>%
    crossing(gene=rownames(es)) %>%
    filter(!str_starts(cellstate, 'Unnamed')) %>%
    group_by(cellstate,gene) %>% summarize(xig=mean(mat(es)[gene[1],.data$cell])) %>%
    group_by(gene) %>% mutate(xig_bar=xig/sum(xig)) %>%
    group_by(cellstate) %>%  summarize(xi=gm_mean(xig_bar))
}

main <- function() {
    ess.norm <- load_data_set('ess.norm')
    clusters.reducedA <- load_data_set('clusters.reducedA')

    sps <- load_data_set("sp.names")
    all_blocks <- load_data_set('all_blocks')

    # compute tau values
    blocks_tau <- crossing(species=sps, block_type=names(all_blocks[[1]])) %>%
        add_column(data=pmap(., function(species, block_type) {
            bs <- all_blocks[[species]][[block_type]]
            taus <- map_dbl(bs, ~calcStat.synBlock(., compute_tau, ess.norm, mat.only=F))
            tibble(tau=taus, species=species, name=names(bs))
        })) %>% unnest()

    save_data_set(blocks_tau)
}

if(sys.nframe() == 0L) {
    main()
}
