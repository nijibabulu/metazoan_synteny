#! /usr/bin/env Rscript

library(devtools)
load_all('.')


# filter the matrix against genes with 0 UMIs and cells with < 200 UMIs
sc <- load_data_set("sc")
exprA.filt <- sapply(sc, filt.mat)
save_data_set(exprA.filt,  overwrite = T)

# create expression set objects (including cellinfo metadata)
sp.names <- load_data_set('sp.names')
cellinfo <- load_data_set('cellinfo')
ess <- sapply(sp.names, function(sp) scExpressionSet(.Mat=exprA.filt[[sp]], .Cellinfo = cellinfo[[sp]]))
save_data_set(ess)

# normalize all the matrices
ess.norm <- sapply(ess, exprNormalize)
save_data_set(ess.norm,  overwrite=T)

# parse synteny clusters
genomes <- load_data_set('genomes')
clusters <- parse_cluster_file(file.path(extPath(), 'nmax5.clust.noplntxids.refmt'), genomes, remove_prefix=T)

# add counts
clusters <- addStat.synClustSet(clusters, ess, sum, 'umi.count')
save_data_set(clusters,  overwrite=T)

# Note: These are the clusters used in Extended Data Figure 1. The total count of 270 in 
# the main text refers to the total number of clusters (277) minus the number of clusters 
# which consist only of genes from a single genome. (4 from Nematostella vectensis, 3 from 
# Amphimedon queenslandica)

# filter the synteny clusters according to their UMI counts
clusters.cleanblocksA <- filterGenesByTotalUMIs.synClustSet(clusters, ess, 0)
save_data_set(clusters.cleanblocksA,  overwrite=T)

# remove synteny clusters with no species remaining and/or fewer than 3 genes
clusters.filtA <- filterClusters(clusters.cleanblocksA, min.species=1, min.genes=3)
clusters.filtA <- addStat.synClustSet(clusters.filtA, ess, sum, 'umi.count')
save_data_set(clusters.filtA,  overwrite=T)

# in cases where a cluster contains multiple blocks per species, pick the one block which has the most UMIs
clusters.reducedA <- filterBlocksByStat(clusters.filtA, 'umi.count')
save_data_set(clusters.reducedA, overwrite=T)

# create synteny block expression sets
synBlockESs.norm <- sapply(names(ess.norm), function(sp) makeSynBlockExprSet(clusters.reducedA, sp, ess.norm[[sp]]))
save_data_set(synBlockESs.norm, overwrite=T)
