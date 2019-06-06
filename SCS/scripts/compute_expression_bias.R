#! /usr/bin/env Rscript
library(tidyverse)
library(ggthemes)
library(cowplot)
library(devtools)
load_all('.')

compute_normalized_expression <- function(es, geo_mean=F) {
cellinfo(es) %>%
    crossing(gene=rownames(es)) %>%
    filter(!str_starts(cellstate, 'Unnamed')) %>%
    group_by(cellstate,gene) %>% dplyr::summarize(xig=mean(mat(es)[gene[1],.data$cell], na.rm=T)) %>%
    group_by(gene) %>% mutate(xig_bar=xig/(if (sum(xig) == 0) 1 else sum(xig)))
}

compute_scaled_bias <- function(es, geo_mean=F) {
  cellinfo(es) %>%
    crossing(gene=rownames(es)) %>%
    filter(!str_starts(cellstate, 'Unnamed')) %>%
    group_by(cellstate,gene) %>% dplyr::summarize(xig=mean(mat(es)[gene[1],.data$cell], na.rm=T)) %>%
    group_by(gene) %>% mutate(xig_bar=xig/(if (sum(xig) == 0) 1 else sum(xig))) %>%
    group_by(gene) %>% mutate(xig_bar=scale(xig, center=T)) %>%
    group_by(cellstate) %>% dplyr::summarize(xi=mean(xig_bar))
}

compute_expression_bias <- function(es, geo_mean=F) {
  .mean <- if(geo_mean) gm_mean else mean
  cellinfo(es) %>%
    crossing(gene=rownames(es)) %>%
    filter(!str_starts(cellstate, 'Unnamed')) %>%
    group_by(cellstate,gene) %>% dplyr::summarize(xig=mean(mat(es)[gene[1],.data$cell], na.rm=T)) %>%
    group_by(gene) %>% mutate(xig_bar=xig/(if (sum(xig) == 0) 1 else sum(xig))) %>%
    group_by(cellstate) %>% dplyr::summarize(xi=.mean(xig_bar))
}

compute_tau <- function(es, include_unnamed=F) {
  cells <- cellinfo(es) %>%
    crossing(gene=rownames(es))
  if(!include_unnamed) {
    cells <- cells %>% filter(!str_starts(cellstate, 'Unnamed'))
  }
  cells %>%
    group_by(cellstate,gene) %>% dplyr::summarize(xig=mean(mat(es)[gene[2],.data$cell], na.rm=T)) %>%
    group_by(gene) %>% mutate(xig_hat=1-xig/max(xig)) %>%
    summarize(taug=sum(xig_hat)/(n()-1)) %>%
    summarize(tau=mean(taug)) %>% flatten_dbl()
}

compute_tau_un <- function(es) compute_tau(es, include_unnamed = T)

main <- function() {
    ess.norm <- load_data_set('ess.norm')
    clusters.reducedA <- load_data_set('clusters.reducedA')

    sps <- load_data_set('sp.names')
    all_blocks <- load_data_set('all_blocks')

    # compute tau values
    # ("taus_unnamed" includes the celltypes which had no name, i.e. "Unnamed N")
    blocks_tau <- tibble(species=sps, block_type="observed") %>%
        add_column(data=pmap(., function(species, block_type){
            bs <- all_blocks[[species]][[block_type]]
            taus_unnamed <- map_dbl(bs, ~calcStat.synBlock(., compute_tau_un, ess.norm, mat.only=F))
            taus <- map_dbl(bs, ~calcStat.synBlock(., compute_tau, ess.norm, mat.only=F))
            tibble(tau=taus, tau_unnamed=taus_unnamed, name=names(bs))
        })) %>% unnest()

    save_data_set(blocks_tau, overwrite=T)

    # compute biases for pca
    biases <- tibble(species=sps) %>%
      add_column(data=pmap(.,function(species) {
        tibble(bias=map(all_blocks[[species]]$observed,
                        ~compute_expression_bias(sc_expr(., ess.norm[[species]])) %>%
                          mutate_if(is.factor, as.character)),
               name=names(all_blocks[[species]]$observed))
      })) %>% unnest() %>% unnest()

    save_data_set(biases)

    scaled_biases <- tibble(species=sps) %>%
      add_column(data=pmap(.,function(species) {
        tibble(bias=map(all_blocks[[species]]$observed,
                        ~compute_scaled_bias(sc_expr(., ess.norm[[species]])) %>%
                          mutate_if(is.factor, as.character)),
               name=names(all_blocks[[species]]$observed))
      })) %>% unnest() %>% unnest()

    save_data_set(scaled_biases)
}

if(sys.nframe() == 0L) {
    main()
}
