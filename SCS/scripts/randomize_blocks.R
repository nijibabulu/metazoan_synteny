#! /usr/bin/env Rscript

library(devtools)
library(fitdistrplus)
load_all()


sp.names <- load_data_set("sp.names")
longnames <- load_data_set("longnames")
clusters.reducedA <- load_data_set("clusters.reducedA")
ess <- load_data_set("ess")

load_data <- function(...) {
  res <- sapply(list(...), function(d) {
    load_data_set(d)
    get(d)
  }, simplify=F, USE.NAMES=T)
  names(res) <- list(...)
  res
  
}
#' Gather all block types for a species into a list of lists named by their randomization method (or lack thereof)
#'
#' @param species the species to load
#'
#' @return a list of lists of observed, shuffled, sampled.all and sampled.ortholog blocks (named respectively)
gather_blocks <- function(species) {
  data <- load_data("clusters.reducedA","clusters.shuffledA",
                    "sampled.allA", "sampled.orthologsA")
  list(
    observed=getSpeciesBlocks(data$clusters.reducedA, species),
    shuffled=data$clusters.shuffledA %>% map(~getSpeciesBlocks(., species)) %>% flatten,
    sampled.all=data$sampled.allA[[species]] %>% flatten,
    sampled.ortholog=data$sampled.orthologsA[[species]] %>% flatten
  )
}

# generate Negative binomial fits to the intervening sizes of the blocks.
get.fit <- function(data) {
  fit.nbinom <- fitdist(data, distr='nbinom', method='mle', lower=c(0,0),
                        start=list(size=.Machine$double.eps,prob=.Machine$double.eps))
}


intervening_sizes_data <- sapply(sp.names, function(sp)  intervening_sizes(clusters.reducedA, sp))
is.df <- data.frame(sizes=unlist(intervening_sizes_data), 
                    species=unlist(sapply(names(intervening_sizes_data), 
                                          function(n) rep(longnames[[n]], times=length(is.list[[n]])))))

fit <- sapply(intervening_sizes_data, get.fit, simplify=F)

save_data_set(intervening_sizes_data)
save_data_set(fits)

# perform the actual shuffling and sampling
set.seed(06102018)
clusters.shuffledA <- shuffle.randomize(clusters.reducedA, ess, n.bs=100)
sampled.allA <- sampleAnalagousBlocks(clusters.reducedA, ess, fits, n.samp=100)
sampled.orthologsA <- sampleAnalagousBlocks(clusters.reducedA, ess, fits, genes.subset=unlist(genes.with.orthology), n.samp=100)

save_data_set(clusters.shuffledA, overwrite = T)
save_data_set(sampled.orthologsA,  overwrite = T)
save_data_set(sampled.allA,  overwrite = T)

# save these blocks in a single data structure with a consistent interface
all_blocks <- sp.names %>% map(gather_blocks) %>% set_names(sps)
save_data_set(all_blocks, overwrite = T)