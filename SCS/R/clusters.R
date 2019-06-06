#' @importFrom future plan
#' @importFrom future.apply future_sapply
NULL

#' @section synClust:
#' methods dealing with synteny clusters
NULL


getSpeciesBlocks.synClust <- function(object,species) {
  species.blocks <- object@.Species[[species]]
  if(is.null(species.blocks)) NULL
  else {
    sapply(ls(species.blocks), get, species.blocks)
  }
}
setMethod("getSpeciesBlocks", signature(object="synClust", species="character"),
          getSpeciesBlocks.synClust)

setMethod("species", signature(object="synClust"), function(object) ls(object@.Species))
setMethod("name", signature(object="synClust"), function(object) object@.Name)


as.data.frame.synClust <- function(x, row.names=NULL, optional=FALSE, ...) {
  do.call(rbind, lapply(getBlocks(x), as.data.frame))
}
setMethod("as.data.frame", "synClust", as.data.frame.synClust)

#' Get all the pairs of synteny blocks from a group in the form of a synPair object
#'
#' @param object synClust
#'
#' @return setGrouper object with the names of the species as the keys
getPairs.synClust <- function(object) {
  .getpair <- function(comb) {
    s1 <- object@.Species[[comb[1]]]
    s2 <- object@.Species[[comb[2]]]

    if(length(s1) > 1 || length(s2) > 1) { stop("Cannot make pairs out of a cluster with multiple blocks per species. Call filterBlocks() first") }
    makePair(s1[[ls(s1)[1]]], s2[[ls(s2)[1]]])
  }
  if(length(object@.Species) < 2) { return(new("setGrouper")) }
  species.combs <- combn(ls(object@.Species), 2, simplify=F)
  Reduce(function(sg, comb) add(sg, comb, .getpair(comb)), species.combs, new("setGrouper"))
}
setMethod("getPairs", "synClust", getPairs.synClust)

addBlockToSpecies.synClust <- function(object, species, block) {
  if(!exists(species, where=object@.Species)) {
    object@.Species[[species]] = new.env()
  }
  object@.Species[[species]][[block@.Name]] <- block
  object
}
setMethod("addBlockToSpecies",
          signature(object="synClust", species="character", block="synBlock"),
          addBlockToSpecies.synClust)

getBlocks.synClust <- function(object) {
  keys <- ls(object@.Blocks)
  sapply(keys, get, object@.Blocks)
}
setMethod("getBlocks", signature(object="synClust"), getBlocks.synClust)

filterBlocks.synClust <- function(object, selectfun) {
  .doselect <- function(species) {
    blocklist <- getSpeciesBlocks(object,species)
    if(length(blocklist) > 1) {
      selectfun(blocklist)
    } else { blocklist[[1]] }
  }
  blocks <- sapply(ls(object@.Species), .doselect)
  Reduce(addBlock, blocks, new("synClust", .Name=object@.Name))
}

#' Filter a synClust down to one block per species
#'
#' @param object synClust object to filter
#' @param selectfun function which takes a list of synBlocks and selects the best one
#'
#' @return new synClust with only one cluster per species
setMethod("filterBlocks", signature(object="synClust", selectfun="function"),
          filterBlocks.synClust)

addBlock.synClust <-  function(object, block) {
  #if(exists(block@.Name, where=object@.Blocks)) { stop(paste("Duplicate block", block@.Name)) }
  object@.Blocks[[block@.Name]] <- block
  object <- addBlockToSpecies(object, block@.Species, block)
  object
}
setMethod("addBlock", signature(object="synClust", block="synBlock"),
          addBlock.synClust)


randomize.synClust <- function(object, genomes) {
  randomized <- new("synClust", .Name=paste0(object@.Name, "rand"))
  blocks <- sapply(getBlocks(object), randomize.synBlock, genomes)
  sapply(blocks, function(b) addBlock(randomized, b))
  randomized
}

#' @section synClustSet :
#' methods dealing with sets of synteny clusters
NULL

#' Construct a synClustSet from a set of synteny clusters
#'
#' @param clusters a list of synteny clusters to generate the clust set from
#' @param ranges a list of granges objects containing the genome genes
#'
#' @return a synClustSet object
#'
#' @export
 make.synClustSet <- function(clusters, ranges) {
  scs <- new("synClustSet", .Genomes=ranges)
  Reduce(addClust.synClustSet, clusters, scs)
}

getClusters.synClustSet <- function(object) {
  keys <- ls(object@.Clusters)
  sapply(keys, get, object@.Clusters)
}
setMethod("getClusters", signature(object="synClustSet"), getClusters.synClustSet)

getPairs.synClustSet <- function(object) {
  sets <- sapply(getClusters(object), getPairs)
  grouper <- Reduce(update.setGrouper, sets, new("setGrouper"))
  new("synPairSet", .Pairs=grouper)
}

#' Get a synPairSet out of a synClust set
#'
#' @param object a synClustSet object
#'
#' @return a synPairSet object
setMethod("getPairs", "synClustSet", getPairs.synClustSet)

setMethod("randomize", signature(object="synClustSet"),
          function(object) {
            randomClusters <- sapply(getClusters(object), randomize.synClust, object@.Genomes)
            make.synClustSet(randomClusters, object@.Genomes)
          })

setMethod("addSpecies", signature(object="synClustSet", species="character"),
          function(object,species) { object@.Species <- union(object@.Species, species); return(object)})

addClust.synClustSet <- function(object,clust) {
  object@.Clusters[[clust@.Name]] <- clust
  object <- addSpecies(object, species(clust))
  add.setGrouper(object@.ClustersBySpecies, sets=species(clust), item=clust)
  object
}
setMethod("addClust", signature(object="synClustSet", clust="synClust"), addClust.synClustSet)

addBlock.synClustSet <- function(object, block) {
  if(!exists(block@.ClusterName, where=object@.Clusters)) {
    cluster <- new("synClust", .Name=block@.ClusterName, .Blocks=new.env(), .Species=new.env())
  } else {
    cluster <- object@.Clusters[[block@.ClusterName]]
  }
  addBlockToSpecies(cluster, block@.Species, block)
  object <- addClust(object, cluster)
  object@.Clusters[[block@.ClusterName]] <- addBlock(cluster, block)
  object
}
setMethod("addBlock", signature(object="synClustSet", block="synBlock"), addBlock.synClustSet)


setMethod("clusters", signature(object="synClustSet"),
          function(object) sapply(ls(object@.Clusters), function(key) object@.Clusters[[key]]))

getBlocks.synClustSet <- function(object, flatten=T) {
  clusters <- clusters(object)
  blocks <- sapply(clusters, function(clust) getBlocks(clust), USE.NAMES = F)
  if(flatten & length(blocks) > 1) {
    blocks <- do.call(c, unlist(blocks, recursive=F))
    names(blocks) <- sapply(blocks, function(b) b@.Name)
  }
  blocks
}
setMethod("getBlocks", signature(object="synClustSet"),
           getBlocks.synClustSet)

getSpeciesBlocks.synClustSet <- function(object, species, flatten=T) {
  clusters <- clusters(object)
  blocks <- sapply(clusters, function(clust) getSpeciesBlocks(clust, species), USE.NAMES = F)
  if(all(sapply(blocks, is.null))) return(NULL)
  if(flatten && length(blocks)) {
    blocks <- unlist(blocks, recursive=F)
    if(is.list(blocks)) {
      blocks <- do.call(c, blocks)
    }
    names(blocks) <- sapply(blocks, function(b) b@.Name)
  }
  blocks
}

setMethod("getSpeciesBlocks", signature(object="synClustSet", species="character"),
          getSpeciesBlocks.synClustSet)

query.synClustSet <- function(object, species, species.count=NULL, strict=F) {
  clusters <- getByGroup(object@.ClustersBySpecies, species)
  if(!strict) {
    others <- setdiff(names(object@.Genomes), species)
    other.sets <- sapply(1:length(others), function(m)
      sapply(combn(others,m, simplify = F), function(ss)
        c(species, ss), simplify=F))
    other.clusters <-  sapply(unlist(other.sets, recursive=F), function(ss)
      getByGroup(object@.ClustersBySpecies, ss), simplify=F)
    clusters <- c(clusters, unlist(other.clusters))
  }
  if(!is.null(species.count)) {
    clusters <- Filter(function(cluster) length(cluster@.Species) == species.count, clusters)
  }
  result <- make.synClustSet(clusters, object@.Genomes)
  result
}


#' Query a synteny cluster set by included species and number of species
#'
#' @param object a synClustSet object
#' @param species a character set with species that must be included in the result
#' @param species.count integer count of the exact number of species required. If missing then return all clusters.
#'
#' @return a synClustSet with the requested species and maximum species counts
#'
#' @export
setMethod("query", signature("synClustSet", "character"), query.synClustSet)

setMethod("speciesGroups", signature(object="synClustSet"), function(object) getGroups(object@.ClustersBySpecies))

setMethod("intersection.data.frame", "synClustSet", function(object) intersection.data.frame(object@.ClustersBySpecies))

setMethod("shuffle.randomize", "synClustSet", function(object, es, n.bs=1000, n.cores=1, corr.fun=function(m) cor(m, method='spearman')) {
    # precompute the correlation of all genes against each other and effectively sample the correlation
    # matrix according to the lengths of the blocks

    blocksBySpecies <- sapply(unlist(species(object)),function (sp) getSpeciesBlocks(object, sp))

    # convenience vectors to indicate the size of the blocks for splitting later
    blockSplits <- sapply(blocksBySpecies, function(bs) rep(seq(length(bs)), times=sapply(bs, length)))

    # names of the genes to sample
    genesBySpecies <- sapply(blocksBySpecies, function(bs) unlist(sapply(bs, genes), use.names = F))

    expr.mats <- sapply(names(genesBySpecies), function(sp) mat(es[[sp]])[genesBySpecies[[sp]],])
    corr.mats <- sapply(expr.mats, function(m) corr.fun(t(as.matrix(m))))

    if(n.cores > 1) {
        plan(multiprocess, workers=n.cores)
        .apply = future_sapply
    }
    else {
        .apply = sapply
    }

    res <- .apply(unlist(species(object)),
                  function(sp) {
                    corr.vals <- do.call(rbind, sapply(1:n.bs, function(x) {
                      # shuffle the matrix by resampling all the genes
                      gene.order <- sample(genesBySpecies[[sp]], length(genesBySpecies[[sp]]))

                      # compute the blocks by splitting the new gene order by the block lengths
                      blocks <- split(gene.order, blockSplits[[sp]])

                      # compute the fisherz-transformed mean correlation of each block
                      sapply(blocks, function(b) {
                        corr.vals <- apply(combn(b, 2), 2, function(p) corr.mats[[sp]][p[1],p[2]])
                        fisherz.inv(mean(fisherz(corr.vals)))
                      })

                    }, simplify = F))
                    corr.vals
                  }
    )
})

filterGenesInBlocks.synClust <- function(object, filter.function) {
  blocks <- sapply(getBlocks(object), filterGenesInBlock, filter.function)
  Reduce(addBlock, blocks, new("synClust", .Name=object@.Name))
}
setMethod("filterGenesInBlocks", "synClust", filterGenesInBlocks.synClust)


filterGenesInBlocks.synClustSet <- function(object, filter.function) {
  filtered.clusters <- sapply(getClusters(object), filterGenesInBlocks, filter.function)
  make.synClustSet(filtered.clusters, object@.Genomes)
}
#' Filter genes in the blocks by an arbitrary fuction
#'
#' @param object a synClustSet
#' @param filter.function a function with the signature function(gene.id, synBlock) -> Boolean
#'
#' @return new synClustSet with modified blocks
#'
#' The callback function is called for each of the blocks and only genes which have a TRUE value will be retained in the blocks.
#'
#' @export
setMethod("filterGenesInBlocks", "synClustSet", filterGenesInBlocks.synClustSet)

#' Filter genes in the blocks by an arbitrary fuction
#'
#' @param object a synClustSet
#' @param ess an scExpressionSet object, columns as cells rows as genes
#' @param min.umis minimum number of total umis to consider
#'
#' @return new synClustSet with modified blocks
#'
#' Filter all genes in blocks for which the total number of UMIs in detected for the gene is below the
#' given threshold. If there is no data for the number of UMIs, the gene is also dropped.
#'
#' @export
filterGenesByTotalUMIs.synClustSet <- function(object, ess, min.tot.umis) {
  valid.gene <- function(gene.id, block) {
    sp.expr <- mat(ess[[species(block)]])
    if(is.null(sp.expr)) {
      return(TRUE) # if we don't have an expression set leave it in.
    }
    if(! gene.id %in% rownames(sp.expr)) {
      return(FALSE)
    }
    else {
      sum(sp.expr[gene.id,]) >= min.tot.umis
    }
  }

  filterGenesInBlocks(object, valid.gene)
}


filterBlocksByStat.synClustSet <- function(object, stat, max.stat=T) {
  filter.function <- function(bs) {
    stat.values <- sapply(bs, function(b) b@.Stats[,stat])
    bs[order(stat.values, decreasing = max.stat)][[1]]
  }
  filterBlocks(object, filter.function)
}

#' Filter a synClustSet down to one block per species based on a calculated statistic
#'
#' @param object synClust object to filter
#' @param stat precomputed stastistic to filter by
#' @param max.stat whether to take the block whose stat is maximum among the set (set false for the minimum)
#'
#' @return new synClustSet with only one cluster per species
setMethod("filterBlocksByStat", signature(object="synClustSet", stat="character"),
          filterBlocksByStat.synClustSet)


#' Filter a synClustSet down to one block per species
#'
#' @param object synClust object to filter
#' @param selectfun function which takes a list of synBlocks and selects the best one
#'
#' @return new synClustSet with only one cluster per species
filterBlocks.synClustSet <- function(object, selectfun) {
  clusters <- sapply(getClusters(object), function(c) filterBlocks(c, selectfun))
  make.synClustSet(clusters, object@.Genomes)
}
setMethod("filterBlocks", signature(object="synClustSet", selectfun="function"),
          filterBlocks.synClustSet)

#' Add a statistic to all blocks in a synClustSet object
#'
#' @param object a synClustSetObject
#' @param expr a list of scExpressionSet objects named by species
#' @param func a callback function which calculates the statistic. It will receive
#'              the expression data associated with the genes in each block only. If
#'              mat.only is set to TRUE, the scExpressionSet will be unwrapped to the
#'              matrix component only. This is ideal for statistics which do not need
#'              cell info, or primitive functions such as `sum`.
#' @param varName a name which will be associated with the statistic.
#' @param mat.only pass the matrix itself only to the `func` (See `func` for details.)
#'
#' @return a amended  synClustSet object containing the corresponding statistics
addStat.synClustSet <- function(object, expr, func, varName, mat.only=T) {
  #blocks <- getBlocks(object)
  blocks <- sapply(getBlocks(object),
                   function (b) addStat.synBlock(b, expr, func, varName, mat.only))
  Reduce(addBlock, blocks, new("synClustSet", .Genomes=genomes(object)))

  #blocks <- sapply(getBlocks(object),
                   #function (b) setStat(b, varName, .calcStat(b, func, expr)))
}

as.data.frame.synClustSet <- function(x, row.names=NULL, optional=FALSE, ...) {
  do.call(rbind, lapply(getClusters(x), as.data.frame))
}
setMethod("as.data.frame", "synClustSet", as.data.frame.synClustSet)

cleanBlocks.synClust <- function(object, min.genes=NULL) {
  blocks <- getBlocks(object)
  if(!is.null(min.genes)) {
   # print(sapply(blocks, function(b) length(genes(b))))
    blocks <- Filter(function(b) length(genes(b)) >= min.genes, blocks)
  }
  Reduce(addBlock, blocks, new("synClust", .Name=object@.Name))
}

setMethod("filterClusters", "synClustSet",
          function(object, min.genes=NULL, min.species=NULL) {
            clusters <- getClusters(object)
            if(!is.null(min.genes)) {
              clusters <- sapply(clusters, cleanBlocks.synClust, min.genes=min.genes)
            }
            if(!is.null(min.species)) {
              clusters <- Filter(function(c) length(c@.Species) >= min.species, clusters)
            }
            make.synClustSet(clusters, object@.Genomes)
})

setMethod("genomes", "synClustSet", function(object) object@.Genomes)

setMethod("species", "synClustSet", function(object) unlist(object@.Species))

setMethod("show", signature(object="synClustSet"),
          function(object) {
            cat("synClustSet with ", length(object@.Clusters), "clusters\n")
            cat("Species: ", unlist(object@.Species), "\n")
            cat("Clusters by Included Species: \n")
            show(object@.ClustersBySpecies)
          })


#' Get a list of  included species in the synClustSet
#'
#' @param object a synClustSet
#'
#' @return a list of characters indicating species names
#' @export
setMethod("species", "synClustSet", function(object) object@.Species)

#' randomize the genes in a synClustSet by shuffling
#'
#' @param object a synClustSet to be randomized
#' @param seed set the random seed
#'
#' @return new synClustSet with all synteny relationships broken
setMethod("randomize", "synClustSet", function(object, seed=180918) {
  set.seed(seed)

  # flatten all the genes of each species into lists
  .genesInSpecies <- function(sp) unlist(sapply(getSpeciesBlocks(object, sp), genes), use.names = F)
  genesBySpecies <- sapply(unlist(species(object)), .genesInSpecies)

  # shuffle all genes in a random order
  randGenesBySpecies <- sapply(names(genesBySpecies),
                               function(sp) sample(genesBySpecies[[sp]], length(genesBySpecies[[sp]])))

  # get the lengths of all the blocks
  blockLengths <- sapply(unlist(species(object)),
                         function(sp) sapply(getSpeciesBlocks(object, sp), length))

  # split the genes into groups by the lengths of the blocks
  .splitGenes <- function(genes, block.sizes) split(genes, rep(seq_along(block.sizes), times=block.sizes))
  randBlockGenes <- sapply(names(randGenesBySpecies),
                           function(sp) .splitGenes(randGenesBySpecies[[sp]], blockLengths[[sp]]))

  clusterNames <- sapply(unlist(species(object)),
                         function(sp) sapply(getSpeciesBlocks(object, sp), clusterName))

  # make synteny blocks out of the randomized gene names
  .makeSynBlock <- function(species, genes, name, clusterName) {
    synBlock(.Name=sprintf('%srand', name),
             .ClusterName=clusterName,
             .Species=species,
             .Conns=list("None(Random)"),
             .Location=paste0("None(Random)"),
             .Length=length(genes),
             .GeneIds=list(genes))
  }
  .synBlockMaker <- function(sp, spGenes, spClusterNames) {
    sapply(seq_along(spGenes),
           function(i) .makeSynBlock(sp, spGenes[[i]], i, spClusterNames[[i]]) )
  }
  blocks <- sapply(names(randBlockGenes),
                   function(sp) .synBlockMaker(sp, randBlockGenes[[sp]], clusterNames[[sp]]))

  # create clusters based on the names of the clusters in the original blocks
  .updateClusts <- function(clusts, block) {
    cname <- clusterName(block)
    if(exists(cname, envir=clusts)) {
      clust <- get(cname, envir=clusts)
    } else {
      clust <- new("synClust", .Name=cname, .Blocks=new.env(), .Species=new.env())
    }
    clust <- addBlock(clust, block)
    assign(cname, clust, envir=clusts)
    clusts
  }
  clustEnvs <- Reduce(.updateClusts, unlist(blocks), new.env())
  clusts <- sapply(ls(clustEnvs), get, clustEnvs)

  # finally, make a synClustSet
  make.synClustSet(clusts, genomes(object))
})

#' Get all intervening sizes of the synteny blocks in a synClustSet.
#'
#' @param object a synClustSet object
#' @param species which species (present in the synClustSet) retreive from
#'
#' @export
setMethod("intervening_sizes", "synClustSet",  function(object, species) {
  unlist(sapply(getSpeciesBlocks(object, species),
                function(b) intervening_sizes(b, genomes(object)[[species]])))
})

#' Sample an analogous set of blocks from the genomes of a synClustSet object
#'
#' @param object a synClustSet object
#' @param expr a species-named list of expression matrices from which to take expression statistics
#' @param fits a species-named list of fit parameters for the intervening genes model. currently supports only the negative binomial distribution.
#' @param genes.subset use sampleBlocks to get random blocks in the genome, subsetted by genes.subset and
# the row contents of the expression matrix expr.
#' @param bs.stat summary statistic for all bootstrapped samples of the genome (defaults to median)
#' @param n.samp number of bootstraps to apply
#' @param seed random seed
#' @param verbose show progress messages
#'
#' @return a list of matrices containing mean pairwise correlation values for each of the sampled blocks (in order analagous to the given observed blocks)
#'
#' @export
sampleAnalagousBlocks.synClustSet <- function(object, expr, fits, genes.subset=NULL, n.samp=1000, seed=03102018, verbose=T) {
  sapply(names(expr), function(sp) {
    if(verbose) cat('\n', sp, '\n')
    sampleBlocks(genomes(object)[[sp]], 4, getSpeciesBlocks(object, sp),
                 intervening.params=as.list(fits[[sp]]$estimate),
                 prefilter = rownames(expr[[sp]]),
                 genes.subset = genes.subset,
                 seed = seed, n=n.samp, verbose=verbose)
  }, simplify=F)
}
setMethod("sampleAnalagousBlocks", "synClustSet", sampleAnalagousBlocks.synClustSet)
