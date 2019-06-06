#! /usr/bin/env Rscript

library(tidyverse)
library(argparse)
library(devtools)
library(reticulate)
library(scDatasets)
library(DrImpute)
library(scImpute)

supported.methods <- function() { c("magic", "drimpute", "scimpute", "scimpute.drop", "unimputed", "unimputed.norm", "unimputed.norm.blocksonly") }

get.script.dir <- function() {
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name = "--file="
    script = sub(file.arg.name, "",  initial.options[grep(file.arg.name,initial.options)])
    dirname(script)
}

preload <- function() {
    load_all(dirname(get.script.dir()))
}

handle.args <- function() {
    parser <- ArgumentParser(description="Plot the results of an imputed matrix.")
    parser$add_argument("-d", "--dir", type="character", default=baseLoc(),
                        metavar="DIR", help=paste("Location of the .rda files [default", baseLoc(),"]"))
    parser$add_argument('-s', "--suffix", default='.impute.cors.rda',
                        help="output to the given suffix [default: .impute.cors.rda]")
    parser$add_argument("-m", "--methods", type="character", default="scimpute,scimpute.drop,drimpute,magic,unimputed,unimputed.norm",
                        help="Load these correlations")
    parser$add_argument("SPECIES", type="character", help="which species to load")

    args <- parser$parse_args()

    args$methods = strsplit(args$methods, ",")[[1]]
    if(any(sapply(args$methods, function(x) !x %in% supported.methods()))) {
      stop(paste("Invalid method listed in", args$methods, ". Supported methods are",
                 paste(supported.methods(), collapse=',')))
    }

    args$param_sets = sapply(args$methods, function(method) param.values(method, args$SPECIES))

    args
}


load.from <- function(file) {
    cat(paste('loading',file,'...\n'))
    e <- new.env()
    load(file, envir=e)
    get(ls(e)[[1]], envir=e)
}

#' save object into an rda file
#'
#' @param obj the object to save
#' @param name the name of the symbol that the object will have
#' @param filename path to a target file
save.as <- function(obj, name, file.name) {
    assign(name, obj)
    save(list=name, file=file.name, compress=F)
}


block.correlation <- function(block, corrs) {
    corrs <- apply(combn(genes(block), 2), 2,
                   function(x) corrs[x[1],x[2]])
    fisherz.inv(mean(fisherz(corrs)))
}


load.data <- function(...) {
    res <- sapply(list(...), function(d) {
               load_data_set(d)
               get(d)
         }, simplify=F, USE.NAMES=T)
    names(res) <- list(...)
    res

}

#' retrieve all the parameter values that are loadable
param.values <- function(method, species) {
   glob <- paste(method, species, "*", "rda", sep=".")
   files <- file.path(baseLoc(), glob) %>%  Sys.glob
   if(length(files) == 0) {
     stop(paste("no files found for species", species, "using method", method))
   }
   files %>%
     basename %>%
     strsplit("[.]") %>%
     map(~.[length(.)-1]) %>%
     unlist %>%
     as.numeric
}

load.imputed.data <- function(method, param, species) {
  base.name <- paste(method, species, param, 'rda', sep='.')
  path <- file.path(baseLoc(), base.name)
  load.from(path)
}

#' Compute the block correlation based on a block and correlation matrix
#'
#' @param block a synBlock object
#' @param cormat a matrix of correlation values
#'
#' @return the correlation value
block.correlation <- function(block, cormat) {
  combs <- block %>% genes %>% combn(2) %>% t() %>% as_tibble(.name_repair = 'minimal')
  corrs <- combs %>% pmap(~if(..1 %in% rownames(cormat) && ..2 %in% rownames(cormat)) cormat[..1, ..2] else 0)
  corrs %>% map(fisherz) %>% flatten_dbl %>% mean %>% fisherz.inv
}

#' Load correlation values from disk, compute the block correlations and return them
#'
#' @param species the species to load
#' @param method the imputation method
#' @param param the parameter of the imputation method
#' @param blocks a list of blocks
#'
#' @return a vector of correlation values
block.correlations <- function(blocks, corrmat) {
  corrs <-
  corrs
}

#' Gather all block types for a species into a list of lists named by their randomization method (or lack thereof)
#'
#' @param species the species to load
#'
#' @return a list of lists of observed, shuffled, sampled.all and sampled.ortholog blocks (named respectively)
gather.blocks <- function(species) {
  data <- load.data("clusters.reducedA","clusters.shuffledA",
                    "sampled.allA", "sampled.orthologsA")
  list(
      observed=getSpeciesBlocks(data$clusters.reducedA, species),
      shuffled=data$clusters.shuffledA %>% map(~getSpeciesBlocks(., species)) %>% flatten,
      sampled.all=data$sampled.allA[[species]] %>% flatten,
      sampled.ortholog=data$sampled.orthologsA[[species]] %>% flatten
    )
}

#' Get a table of methods x blocktypes x parameters of correlation values
#'
#' @param methods a character vector of imputation methods to include. This will be passed to the block.correlations() method
#' @param blocks list of list of block objects named by the block typ (e.g. observed, shuffled)
#' @param species a species name. this is only used for the species type
#'
#' @return a tibble with all block correlations
gather.imputed <- function(methods, blocks, species) {
  meth.tbl <- crossing(
    species=species,
    method=methods,
    block_type=names(blocks)
  )

  param.tbls <- meth.tbl %>%
    select(method) %>%
    pmap(~tibble(param=param.values(..1,species)))


  tbl <- add_column(meth.tbl, param.tbls) %>%
    unnest() %>%
    arrange(method,param)

  cor.tbls <- tbl %>%
      group_by(method,param) %>% nest() %>%
      pmap(function(method, param, data) {
          expr.data <- load.imputed.data(method, param, species)
          cormat <- expr.data$cormats[[species]]
          tbls <- data %>%
              select(block_type) %>%
              pmap(~tibble(
                   name=names(blocks[[.]]),
                   cor=blocks[[.]] %>% map(~block.correlation(., cormat)) %>% flatten_dbl))
          rm(expr.data)
          gc()
          tbls
      })  %>% flatten()

  tbl %>% add_column(cor.tbls) %>% unnest()
}

main <- function() {
    preload()
    args <- handle.args()

    # get all blocks into lists
    blocks <- gather.blocks(args$SPECIES)
    imputed <- gather.imputed(args$methods, blocks, args$SPECIES)

    save(imputed, file=file.path(baseLoc(),paste0(args$SPECIES, args$suffix)))
}

#save(list=args$IDENTIFIER, file=out.path, compress=F)
if(sys.nframe() == 0L) {
    main()
}
