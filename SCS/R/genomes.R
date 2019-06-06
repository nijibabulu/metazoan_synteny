#' @importFrom fastmatch %fin%
#' @importFrom GenomicRanges seqnames 
#' @importFrom IRanges start end
#' @importFrom GenomeInfoDb sortSeqlevels dropSeqlevels
NULL

#' Generate all blocks in a genome, excluding any which may overlap with a given set of blocks
#'
#' @param genome a GRanges object representing all genes in the genome
#' @param size the size of the blocks to be generated
#' @param nmax maximum number of intervening genes in the blocks
#' @param exclude list of synBlocks to exclude from the generated blocks. All windows containing these genes will be excluded as well.
#'
#' @return a list of gene names representing all qualifying blocks
allBlocks <- function(genome, size, nmax, exclude=c()) {
  .generateBlocks <- function(idxm, genes) {
    all.genes <- genes[min(idxm):max(idxm)]
    if(any(sapply(exclude.genes, function(g) g %in% all.genes))) { return(NULL) }
    apply(idxm, 2, function(i) genes[i])
  }
  .getBlocks <- function(genes, win.size) {
    if(length(genes) < win.size) { return() }
    # create a list of matrices of indices that contains all combinations of genes
    # within the window and length of the genes on this chromosome
    idxms <- lapply(seq(length(genes)-win.size), function(i) combn(i:(i+win.size), size))
    # convert this into a matrix of gene names
    blm <- do.call(cbind,lapply(idxms, .generateBlocks, names(genes)))
    # if the genes are all null (because they are already within the excluded blocks), t
    # then we can't do anything with this matrix
    if(is.null(blm)) return(NULL)
    # convert this into a list of gene names
    split(blm, rep(1:ncol(blm), each=nrow(blm)))
  }
  sorted <- sortSeqlevels(genome)
  sorted <- sort(sorted)
  sorted.by.chr <- split(sorted, seqnames(sorted))
  win.size <- size*nmax

  exclude.genes <- unlist(sapply(exclude, genes))

  do.call(c, lapply(sorted.by.chr, .getBlocks, win.size))
}

#' Sample blocks from a genome.
#'
#' @param genome a GRanges object containing all genes in a genome
#' @param nmax maximum number of intervening seqeunces
#' @param nonrandom.blocks a set of blocks to avoid overlapping in generating the random sample sapce
#' @param genes.subset if not NULL, only sample blocks which have the given genes (useful if you would, for example, want only genes which have an ortholog to another gene)
#' @param prefilter remove some genes from the genome before sampling (e.g. genes with no expression data)
#' @param intervening.rfunc random variable generator function for intervening genes (default is negative binomial)
#' @param intervening.params parameters for the random function
#' @param max.tries maximum tries to generate a block. if overstepped, the script is stopped.
#' @param seed random seed
#' @param n number of samples to generate
#'
#' @return a list of synBlock objects. If n > 1, a list of lists of synBlock objects
sampleBlocks <- function(genome, nmax, nonrandom.blocks, prefilter=NULL, genes.subset=NULL, max.tries=1000000000,
                         intervening.rfunc=rnbinom, intervening.params=list(size=1, prob=1), seed=030918, n=1,
                         check.visited=F, verbose=T) {
  set.seed(seed)

  filtered <- if(is.null(prefilter)) genome else genome[names(genome) %in% as.character(prefilter)]

  sorted <- sortSeqlevels(filtered)
  sorted <- sort(sorted, ignore.strand=T)
  sorted.by.chr <- as.list(split(sorted, seqnames(sorted)))

  shortest.block <-  min(sapply(nonrandom.blocks, length))
  unusable.seqnames <- names(which(sapply(sorted.by.chr, length) < shortest.block))

  # remove chromosomes which do not have enough genes on them to contain a block
  filtered <- dropSeqlevels(sorted, unusable.seqnames, pruning.mode = 'coarse')
  filtered <- sortSeqlevels(filtered)
  filtered <- sort(filtered, ignore.strand=T)
  filtered.by.chr <- as.list(split(filtered, seqnames(filtered)))

  # get all ranges
  if(check.visited) {
    visited <- suppressWarnings(do.call(c, sapply(unlist(nonrandom.blocks, use.names = F), grange)))
    seqlevels(visited) <- seqlevels(genome)
  }

  genome.size <- length(filtered)
  chrsizes <- sapply(filtered.by.chr, length)
  # index of the chrom implied by a random number from 1:genome.size
  chr.idxs <- rep(1:length(chrsizes), chrsizes)
  # the global start index of each scaffold
  chr.starts <- as.integer(cumsum(chrsizes))

  get.next <- function(block) {
    size <- length(block)
    for(i in seq(max.tries)) {
      idx <- sample.int(n=genome.size, size=1)
      chr.idx <- chr.idxs[idx]
      # the position chromosome starts at is the cumulative sum of all the chromosomes prior to the current one
      # find the index of this current chromosome to start on
      start.on.chr <- 1+idx- max(0,chr.starts[chr.idx-1])
      cur.chr <- filtered.by.chr[[chr.idx]]

      # get all indices of the block by sampling the poisson distribution as many times as the size of the block
      .next.size <- function() { n <- do.call(intervening.rfunc, args=c(n=1,intervening.params)) ; if(n <= nmax)  n  else .next.size() }
      intervening.sizes <- replicate(size-1, .next.size())
      idxs <- start.on.chr + c(0, cumsum(intervening.sizes + 1))

      # if this goes off the boundaries of the chromosome then give up
      if(idxs[size] > chrsizes[chr.idx]) { next }
      rchrom <- names(chrsizes)[chr.idx]
      rstart <- start(cur.chr)[start.on.chr]
      rend <- end(cur.chr)[idxs[size]]
      gr <- GRanges(rchrom, ranges=IRanges(start=rstart, end=rend))

      if(check.visited) {
        if(length(findOverlaps(gr, visited))) {  next }

        # add the interval to the overlaps
        visited <<- suppressWarnings(c(visited, gr))
      }

      gene.ids <- names(cur.chr)[idxs]

      # check if all the genes are present in the genes.subset, if given
      if(!is.null(genes.subset) && !all(sapply(gene.ids, function(g) g %fin% genes.subset))) { next }

      # return the randomized block
      return(synBlock(.Name=sprintf('%srand', block@.Name),
               .ClusterName=block@.ClusterName,
               .Species=block@.Species,
               .Conns=list("None(Random)"),
               .Location=paste0(as.character(rchrom), ':', rstart, '..', rend),
               .Length=size,
               .GeneIds=list(gene.ids)))
    }
    warning(paste("Could not find a block for", name(block)), "after", max.tries, "tries.")
  }

  if(n > 1) {
    if(verbose) pb <- txtProgressBar(min=1, max=n, style=3)
    res <- sapply(1:n, function(x) {
      if(verbose) setTxtProgressBar(pb, x)
      sapply(nonrandom.blocks, get.next)
    }, simplify=F)
  } else {
    res <- sapply(nonrandom.blocks, get.next)
  }

  res
}
