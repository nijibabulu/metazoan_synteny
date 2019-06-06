# NOTE This is necessary for the Matrix methods to correctly dispatch
#' @importFrom purrr map
#' @importFrom dplyr distinct mutate group_by
#' @importFrom Matrix t colSums rowSums
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom ggrepel geom_text_repel
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid theme_cowplot get_legend
#' @importFrom data.table fread
#' @importClassesFrom Matrix Matrix dgCMatrix
NULL

#' Retrieve the cell information from a single cell expression set object
#'
#' @param object a scExpressionSet object
#'
#' @export
setMethod("cellinfo", "scExpressionSet", function(object) object@.Cellinfo)

#' Retrieve the expression matrix of a single cell expression set
#'
#' @param object a scExpressionSet object
#'
#' @export
setMethod("mat", "scExpressionSet", function(object) object@.Mat)

#' Subset an expressionSet
#'
#' @param x an scExpressionset object
#' @param i rows (genes) to include in the subset
#' @param j columns (cells) to include in the subset
#' @param drop passed to `[`
#'
#' @export
setMethod("[", "scExpressionSet", function(x, i, j, ..., drop = FALSE) {
  if(missing(drop)) drop <- FALSE
  if(missing(i) && missing(j)) stop("call to `[` on scExpressionSet with no parameters")
  if(!missing(j)) {
    x@.Cellinfo <- x@.Cellinfo[j,, ..., drop = drop]
    x@.Mat <- x@.Mat[,j, ..., drop=drop]
  }
  if(!missing(i)) {
    x@.Mat <- x@.Mat[i,, ..., drop=drop]
  }
  x
})

setMethod("show", "scExpressionSet", function(object) {
  cat("scExpressionSet object\n")
  dims <- dim(object@.Mat)
  cat(sprintf(' %d x %d sparse Matrix of class "%s"\n', dims[1], dims[2], class(object@.Mat)))
  cat(sprintf(' cellinfo with slots: %s\n', paste(colnames(object@.Cellinfo), collapse=', ')))
  invisible(sapply(colnames(object@.Cellinfo),
                   function(n) cat(sprintf('  unique %s: %d\n', n, length(unique(object@.Cellinfo[,n]))))))
})

#' @export
setMethod("dim", "scExpressionSet", function(x) dim(x@.Mat))

#' @export
setMethod("dimnames", "scExpressionSet", function(x) dimnames(x@.Mat))

#setMethod("rownames", "scExpressionSet", function(x) rownames(x@.Mat))
#setMethod("colnames", "scExpressionSet", function(x) colnames(x@.Mat))
setMethod("rbind2", signature(x="scExpressionSet",y="scExpressionSet"),
          function(x,y) {
            if(!identical(x@.Cellinfo, y@.Cellinfo)) {
              stop("cannot rbind on two scExpressionSets with differing cellinfos")
            }
            scExpressionSet(.Mat=rbind(mat(x),mat(y)), .Cellinfo=x@.Cellinfo)
          }
)

#' Read the GSM files from GEO archives
#'
#' @param filename GSM formatted file (basic matrix tab-separated matrix)
#' @param rowname.in.header set to T if there is a header item which correspondes to the row names (e.g. geneId)
#'
#' @return matrix of integers
#'
#' @export
read.gsm <- function(filename, rowname.in.header=F) {
  header <- strsplit(readLines(filename, 1), '\t')[[1]]
  if(rowname.in.header) header <- header[2:length(header)]
  tab <- fread(filename, skip = 1, header=F,
               colClasses=c('character', rep('integer', length(header))),
               col.names=c('geneid', header))
  # this is hugely inefficient but necessary (so far) because Matrix does not accept a data.frame
  m <- as.matrix(tab[,-1])
  rownames(m) <- tab$geneid
  Matrix(m, sparse=T)
}


# in the analysis Arnau does, he filters cells with < 200 UMIs and rows with less than 2 reads.
# clusters <- parse_cluster_file(file.path(extPath, 'nmax5.clust.reclust'))
filt.mat <- function(m, min.cell.umis=200, min.gene.umis=0) {
  m[rowSums(m) >= min.gene.umis, colSums(m) >= min.cell.umis]
}

#' Sample random rows in an expression set and summarize them
#'
#' @param object an scExpressionSet object
#' @param n.genes number of genes from the expression set to sample
#' @param n.samp number of times to sample the expression set
#' @param aggregate an aggregation function to perform on each sample (by default, sum the columns with Matrix::colSums)
#' @param transform transformation to apply to the resulting aggregated data (e.g. take a log)
#'
#' @return a dgCmatrix with the sampled statistics
#'
#' @export
sample.scExpressionSet <- function(object, n.genes, n.samp, aggregate=Matrix::colSums, transform=exprLog.dgCMatrix) {
  sapply(
    1:n.samp,
    function(x) {
      sampled.mat <- mat(object)[sample.int(nrow(object), n.genes),]
      sampled.mat.transsums <- transform(Matrix::colSums(sampled.mat))
      res <-setNames(sampled.mat.transsums, seq_along(colnames(object)))
      res
    })
}


#' Filter an scExpression set by various criteria
#'
#' @param object an scExpressionSet object
#' @param min.metacell.size the minimum number of cells associated to a metacell in order to retain it
#' @param min.cellstate.size the minimum number of cells associated to a cellstate in order to retain it
#' @param cellstate.name.exclude a pattern string for names of cellstates to exclude
#'
#' @return a filtere scExpressionSet object
#'
#' @export
filter.scExpressionSet <- function(object, min.metacell.size=0, min.cellstate.size=0, cellstate.name.exclude=NULL) {
  .filterSize <- function(object, colname, min.val) object %>% add_count(!!as.symbol(colname)) %>%
    filter(n >= min.val) %>% dplyr::select(-one_of('n')) %>% as.data.frame %>%
    as.data.frame(row.names=as.character(.$cell))
  object@.Cellinfo <- cellinfo(object) %>%
    .filterSize('cellstate',min.cellstate.size) %>%
    .filterSize('metacell',min.metacell.size)
  if(!is.null(cellstate.name.exclude)) {
    object@.Cellinfo <- object@.Cellinfo[grep(cellstate.name.exclude, cellinfo(object)$cellstate, invert = T),]
    object@.Cellinfo <- droplevels(object@.Cellinfo)
  }
  object@.Mat <- object@.Mat[,rownames(object@.Cellinfo)]
  object
}

#' Aggregate an expression set by phenotype metadata
#'
#' @param object an scExprtessionSet object
#' @param level the level to aggregate by. must be a column in the cellinfo slot (e.g. cellstate or metacell)
#' @param aggregate function for aggregation (e.g. sum)
#'
#' @return an abbreviated scExpressionSet
#'
#' @export
aggregate.scExpressionSet <- function(object, level, aggregate=sum, min.cells=0) {
  if(!level %in% colnames(cellinfo(object))) {
    stop(sprintf("level %s does not exist in given expression set.", level))
  }


  sorted <- object[,order(cellinfo(object)[,level])]
  groups <- split(cellinfo(sorted), cellinfo(sorted)[,level])
  if(min.cells > 0) {
    groups <- groups[lapply(groups, nrow) > min.cells]
  }
  mat <- do.call(cbind, sapply(groups, function(group)
    apply(cbind(matrix(0,nrow(object)),
                mat(object)[,rownames(group)]),
          1,aggregate),
    simplify = F
  )
  )

  # create a new cell info object with the information replaced
  level.szs <- sapply(colnames(cellinfo(object)),
                      function(c) length(unique(cellinfo(object)[,c])))
  level.hier <- names(level.szs)[order(level.szs)]

  ci <- cellinfo(object) %>% dplyr::select(one_of(level.hier[1:which(level.hier == level)])) %>% distinct()
  missing <- colnames(cellinfo(object))[!colnames(cellinfo(object)) %in% colnames(ci)]
  sapply(missing, function(m) {ci[,m] <<- ci[,level]})

  rownames(ci) <- ci[,level]

  scExpressionSet(.Mat=as(mat, 'dgCMatrix'), .Cellinfo=ci)
}

plotMat.synBlockExprSet <- function(
  object,
  es=NULL,
  include.unnamed=F,
  seed=240918,
  n.samp=1000,
  colors=c("white","orange","red","purple","black"), #"white",
  na.color="gray",
  plot.colors=F,
  text.cex =0.6,
  cell.label.size=3,
  cell.label.bottom.margin=-7,
  cell.label.upper.limit=0.5,
  scale=F,
  cell.label.proportion=0.1,
  cell.segement.pointer.color='black',
  min.cellstate.cells=0,
  min.metacell.cells=0,
  summary.limits=c(7,9),
  show.cloud=T,
  summary.breaks=1,
  show.block.names=F,
  show.block.ids=F,
  show.gene.names=F,
  show.points=T,
  show.summary=T,
  background.samples=0,
  cell.column='cell',
  metacell.column='metacell',
  minor.breaks.column=metacell.column
) {
  stopifnot(class(object) == 'synBlockExprSet')
  set.seed(seed)

  if(background.samples > 0 && is.null(es)) {
    stop("need an expression matrix for background samples.")
  }
  if(show.block.names && show.block.ids) {
    stop("cannot show both the block names and ids.")
  }

  blockExprs <- getBlockExprs.synBlockExprSet(object)

  plot.es <- do.call(rbind, sapply(blockExprs, getExpr.synBlockExpr))

  plot.es <- filter.scExpressionSet(plot.es,
                                   min.metacell.size = min.metacell.cells,
                                   min.cellstate.size =  min.cellstate.cells)

  if(!include.unnamed) {
    plot.es <- filter.scExpressionSet(plot.es, cellstate.name.exclude = 'Unnamed')
  }

  if(scale) {
    #m = apply(mat(plot.es), 1, mean, na.rm=T)
    #s = apply(mat(plot.es), 1, sd, na.rm=T)
    #mat.disp <- (mat(plot.es) - m) / s
    mat.disp <- scale(mat(plot.es), center=T)
    mat.disp[as.matrix(mat(plot.es)) == 0] = NA
    lim <- min(-min(mat.disp, na.rm=T), max(mat.disp, na.rm=T))
    mat.disp[mat.disp < -lim] = -lim
    mat.disp[mat.disp > lim] = lim
  } else {
    mat.disp <- exprLog.dgCMatrix(mat(plot.es))
  }

  .getBreaks <- function(labels) cumsum(rle(as.character(labels))$lengths)
  .getBreakCenters <- function(breaks) {
    sapply(seq_along(breaks), function(i) (breaks[i]+max(0,breaks[i-1]))/2)
  }

  cellBreaks <- .getBreaks(cellinfo(plot.es)[,minor.breaks.column])
  cellstateBreaks <- .getBreaks(cellinfo(plot.es)$cellstate)
  cellstateLabelPositions <- .getBreakCenters(cellstateBreaks)
  cellstateLabels <- ifelse(grepl('Unnamed', levels(cellinfo(plot.es)$cellstate)),'', levels(cellinfo(plot.es)$cellstate))

  synBreaks <- cumsum(sapply(blockExprs, length.synBlockExpr))
  synLabels <- blockExprs %>% map(getBlock.synBlockExpr) %>% map(geneSummary)
  synNames <- blockExprs %>% map(getBlock.synBlockExpr) %>% map(name)
  synLabelBreaks <- rownames(mat.disp)[round(.getBreakCenters(synBreaks))]

  geneLabels <- unlist(sapply(blockExprs, function(b) geneNames(getBlock.synBlockExpr(b))))

  empty.theme <- theme(axis.title = element_blank(), axis.text=element_blank(),
                       axis.ticks=element_blank(), axis.line=element_blank(),
                       legend.position='none', plot.margin = unit(c(0,7,0,0),'pt'))
  if(scale) {
    #lim <- 2*min(-min(mat.disp, na.rm=T), max(mat.disp, na.rm=T))
    scale.fill <- scale_fill_gradientn(colours=colors, limits=c(-lim,lim), na.value=na.color)
    #scale.fill <- scale_fill_gradientn(colours=colors)
  } else {
    ulim <- ceiling(max(mat.disp, na.rm=T))
    scale.fill <- scale_fill_gradientn(colours=colors, limits=c(0,ulim),breaks=c(0,ulim),
                                       values=c(0,seq(.Machine$double.eps,1,length.out=length(colors)-1)),
                                       na.value=na.color)
  }
  mat.m <- merge(melt(as.matrix(mat.disp)), cellinfo(plot.es), by.x='Var2', by.y=cell.column, all.y=F)
  mat.m$Var2 <- factor(as.character(mat.m$Var2), cellinfo(plot.es)[,cell.column] )

  # plot the cell names
  p.points <- ggplot(data.frame(x=cellstateLabelPositions, label=cellstateLabels),
                     aes(x=x,y=0,label=label)) +
    geom_point(size=0) + scale_x_continuous(expand=c(0,0), limits = c(0,ncol(plot.es))) +
    geom_text_repel(nudge_y=cell.label.upper.limit/2, size=cell.label.size, segment.size=0.25,
                    segment.colour = cell.segement.pointer.color) +
    lims(y=c(0,cell.label.upper.limit)) +
    theme(axis.title = element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          legend.position='none', plot.margin = unit(c(0,7,cell.label.bottom.margin,0),'pt'))
  #theme(axis.title = element_blank(), axis.text.y=element_blank(),
  #legend.position='none', plot.margin = unit(c(0,7,0,0),'pt'))
  #nudge_y = 0.05, direction = "x", segment.size = 0.2
  p <- ggplot(mat.m, aes(x=Var2, y=Var1, color=cellstate)) +
    geom_raster(aes(fill=value)) +
    scale.fill +
    geom_hline(yintercept=synBreaks[1:length(synBreaks)-1]+0.5,size=0.2) +
    geom_vline(xintercept=cellBreaks+0.5,color='grey',size=0.1) +
    geom_vline(xintercept=cellstateBreaks+0.5,color='black',size=0.5) +
    theme_bw() +
    empty.theme +
    theme(axis.text.x = element_text(size=9, angle=45, hjust=0),
          panel.border=element_rect(color="black",fill=NA,size=1),
          legend.position = 'top', legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,-10,-20,-10), legend.justification='center') +
    guides(fill=guide_colorbar(title='Normalized log expression',
                               title.position='top',
                               title.theme = element_text(size=10, hjust=0.5),
                               ticks=F, barheight = 0.5, barwidth=8,
                               frame.colour='black',frame.linewidth = 1,
                               label.theme = element_text(size=9)
    ), color=guide_legend())


  p.legend <- get_legend(p)
  p <- p + theme(legend.position='none')

  if(show.gene.names || show.block.names || show.block.ids) {
    p <- p + theme(axis.text.y=element_text()) + if(show.gene.names) {
      scale_y_discrete(breaks=synLabelBreaks, labels=unname(synLabels))
    } else if(show.block.ids) {
      scale_y_discrete(breaks=synLabelBreaks, labels=unname(synNames))
    } else {
      scale_y_discrete(labels=unname(geneLabels))
    }
  }


  # problem with this may have been mostly that the normalization was not rational before. there was a per-column normalization (plus per-gene) that was done *after* the matrix was subsetted. this would also be incompatible with the sampling.
  if(show.summary) {
    matsum.m <- mat(plot.es) %>% colSums %>% exprLog.dgCMatrix %>% melt
    matsum.m$cell <- seq_along(rownames(matsum.m))

    p.summ <- ggplot(matsum.m, aes(x=cell,y=value)) +
      #stat_summary_bin(fun.y = "mean", geom = "point", bins=300, size=0.5) +
      scale_fill_gradientn(colours=c("white","#313695"), na.value=na.color) +
      scale_y_continuous(expand=c(0,0), limits=summary.limits, breaks=seq(summary.limits[1],summary.limits[2],summary.breaks)) +
      scale_x_continuous(expand=c(0,0)) +
      labs(y='Total Expression') +
      theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(),
            axis.text.y=element_text(size=8), axis.title.y=element_text(size=8),
            # axis.ticks.y=element_blank(),
            axis.title.x=element_blank(), #
            panel.border=element_blank(), legend.position='none')
    if(show.cloud) {
      p.summ <- p.summ + stat_density2d(aes(fill = ..density..^2), geom = "tile", contour = FALSE, n = 200)
    }
    if(show.points) {
      p.summ <- p.summ + geom_point(alpha=0.5, size=0.0025)
    }

    if(background.samples > 0) {
      set.seed(seed)
      es.filt <- es[,colnames(plot.es)] # expr[rowSums(expr) > min(rowSums(mat)),]

      matsamples.m <- melt(sample.scExpressionSet(es.filt, n.genes=nrow(plot.es), n.samp=n.samp))
      p.summ <- p.summ +
        geom_boxplot(data=matsamples.m, aes(x=Var1, y=value, group=Var1), fill='darkgray', linetype=0,  outlier.shape=NA)
    }
    data.grid <- plot_grid(ggdraw(p.legend),p.points,p,p.summ,ncol=1,align='v', axis='lr',rel_heights = c(0.15,cell.label.proportion,0.55,0.2))
  }
  else {
    data.grid <- plot_grid(ggdraw(p.legend),p.points,p,ncol=1,align='v', axis='lr',rel_heights = c(0.15,cell.label.proportion,0.75))
  }

  data.grid
}

### TODO: get this working with new refactoring!

plotMat.synClusterSet <- function(
  object,
  species,
  exprSets,
  min.cellstate.cells=0,
  scale=F
) {

  # retrieve the block expressions for teh corresponding species
  blockExprs <- sapply(species, function(sp) getBlockExprs.synBlockExprSet(exprSets[[sp]]))

  # subset those expressions by the ones that are available in *all* clusters
  clusters <- query(object,species)
  blockSubsets <- sapply(species, function(sp) {
    sapply(getClusters(clusters), function(c) {
      sapply(names(getSpeciesBlocks(c, sp)), function(n)
        if(n %in% names(blockExprs[[sp]])) blockExprs[[sp]][[n]] else NA)
    })
  }, simplify=F)

  if(any(sapply(blockSubsets, length) == 0)) {
    warning(paste(c("no plot for species pair", species), collapse=" "))
    return(ggplot(data.frame()) + geom_blank())
  }

  # construct scExpressionSets for each species containing the filtered blockExpr objects

  .safe.rbind <- function(...) { l <- list(...); if(length(l) > 1) do.call(rbind, l) else l[[1]] }

  ess <- sapply(species, function(sp)
    do.call(.safe.rbind, sapply(blockSubsets[[sp]], getExpr.synBlockExpr, simplify = F)))

  if(length(min.cellstate.cells) == 1) {
    min.cellstate.cells <- setNames(rep(min.cellstate.cells, length(species)), species)
  } else {
    min.cellstate.cells <- setNames(min.cellstate.cells, species)
  }

  plots <- sapply(species, function(sp)
    plotMat.scExpressionSet(ess[[sp]], blockSubsets[[sp]], show.summary=F, show.block.ids = T,
                            min.cellstate.cells=min.cellstate.cells[[sp]], scale=scale),
    simplify=F)

  plot_grid(plotlist=plots,nrow=1)
}



#' Reorder a set of synteny blocks by their correlation
#'
#' @param object a synBlockExpressionSet
#' @param method one of 'pearson' for pearson correlation or 'bicor' for biweight midcorrelation a la WGCNA
#'
#' @return a reordered synteny expression set
#'
#' @export
reorderByCorr.synBlockExprSet <- function(object, dist.method='pearson', clust.method='ward.D2') {
  methods <- list(
    'maxpearson'=maxpearson.synBlockExprSet,
    'pearson'=pearson.synBlockExprSet,
    'bicor'=bicor.synBlockExprSet)
  if(!dist.method %in% names(methods)) {
    stop(paste("Method for reorderByCorr.synBlockExprSet must be one of", paste(names(methods))))
  }

  cmat <- methods[[dist.method]](object)
  hc <- hclust(as.dist(1-cmat), method=clust.method)
  be <- object@.BlockExprs[hc$order]

  # since hclust is an S3 class, we can either make a fake hclust class, import something else,
  # or just use a representational workaround. I chose the latter.
  new("synBlockExprSet", .Species=object@.Species, .BlockExprs = be, .Corr=cmat, .HClust=list(hc))
}





makeSynBlockExpr <- function(block, es) {
  stopifnot(class(block) == 'synBlock')

  new("synBlockExpr", .Block=block, .Expr=sc_expr(block, es))
}


#' Create a synBlockExpr object from a set of synteny clusters on a single species.
#'
#' @param scs a synClustSet object
#' @param species a character vector that has the species name
#' @param es an scExpressionSet object
#'
#' @return a synBlockExprSet object
#'
#' @export
makeSynBlockExprSet <- function(scs, species, es, min.genes=2) {
  stopifnot(class(scs) == 'synClustSet')
  stopifnot(class(species) == 'character')

  # only take blocks with at least 3 genes in the matrix present
  blocks <- Filter(function(b) sum(genes(b) %in% rownames(es)) > min.genes,
                   getSpeciesBlocks(scs, species))

  blockExprs <- sapply(blocks, function(b) makeSynBlockExpr(b, es))


  new("synBlockExprSet",
      .Species=species,
      .BlockExprs=blockExprs)
}


setMethod(show, signature(object="synBlockExpr"),
          function(object) cat("synBlockExpr with", nrow(object@.Expr), "genes"))

setMethod(show, signature(object="synBlockExprSet"),
          function(object) cat("synBlockExprSet for species", object@.Species,
                               "with", length(object@.BlockExprs), "blocks"))

.checkReordered.synBlockExprSet <- function(object) {
  if(length(object@.Corr) == 0 || length(object@.HClust) == 0) {
    stop("synBlockExprSet is not reordered. Call reorderByCorr.synBlockExprSet first.")
  }
}
getBlockExprs.synBlockExprSet <- function(object) object@.BlockExprs
getHClust.synBlockExprSet <- function(object) { .checkReordered.synBlockExprSet(object);  object@.HClust[[1]] }
getCorrMat.synBlockExprSet <- function(object) {.checkReordered.synBlockExprSet(object) ;  object@.Corr }

getExpr.synBlockExpr <- function(object) object@.Expr
getBlock.synBlockExpr <- function(object) object@.Block
length.synBlockExpr <- function(object) nrow(object@.Expr)

getBlocks.synBlockExprSet <- function(object) {
  stopifnot(class(object) == 'synBlockExprSet')
  sapply(object@.BlockExprs, getBlock.synBlockExpr)
}

blockLengths.synBlockExprSet <- function(object) {
  stopifnot(class(object) == 'synBlockExprSet')
  sapply(getBlocks.synBlockExprSet(object), length)
}

geneSummaries.synBlockExprSet <- function(object) {
  stopifnot(class(object) == 'synBlockExprSet')
  sapply(getBlocks.synBlockExprSet(object), geneSummary)
}

corr.synBlockExpr <- function(object, other, corr.func, stat='mean') {
  stopifnot(class(object) == 'synBlockExpr')
  stopifnot(class(other) == 'synBlockExpr')
  stopifnot(stat %in% c('mean','max'))

  corrs <- c()
  for(i in seq(nrow(object@.Expr))) {
    for(j in seq(nrow(other@.Expr))) {
      v <- suppressWarnings(corr.func(mat(object@.Expr)[i,], mat(other@.Expr)[j,]))
      corrs <- c(corrs, v)
    }
  }

  if(stat == 'mean') {
    fisherz.inv(mean(fisherz(corrs)))
  } else if(stat == 'max') {
    max(corrs)
  }
}

bicor.synBlockExpr <- function(object, other) corr.synBlockExpr(object, other, bicor)
pearson.synBlockExpr <- function(object, other) corr.synBlockExpr(object, other, cor)
maxpearson.synBlockExpr <- function(object, other) corr.synBlockExpr(object, other, cor, stat='max')


corr.synBlockExprSet <- function(object, corr.func) {
  stopifnot(class(object) == 'synBlockExprSet')

  n <- length(object@.BlockExprs)
  cmat <- diag(1, n, n)

  for(i in seq(n-1)) {
    for(j in seq(i+1, n)) {
      v <- corr.func(object@.BlockExprs[[i]], object@.BlockExprs[[j]])
      if(is.nan(v) || is.na(v)) { v <- 0 }

      cmat[i,j] <- v
      cmat[j,i] <- v
    }
  }

  cmat
}

bicor.synBlockExprSet <- function(object) corr.synBlockExprSet(object, bicor.synBlockExpr)
pearson.synBlockExprSet <- function(object) corr.synBlockExprSet(object, pearson.synBlockExpr)
maxpearson.synBlockExprSet <- function(object) corr.synBlockExprSet(object, maxpearson.synBlockExpr)


#' Normalize by umis/cell
#'
#' @param mat matrix of umi counts with genes in rows and cells in columns
#' @param umi_scale the initial scaling factor after normalizing for UMIs.
#'
#' @return UMI-normalized matrix
exprNormalize.dgCMatrix <- function(mat, umi_scale=1000) {
  # normalize by the total umis per cell
  totu = pmax(1,colSums(mat))
  mat = t(t(mat)/totu)*umi_scale

  mat
}

setMethod("exprNormalize", "scExpressionSet", function(object, ...)
  scExpressionSet(.Mat=exprNormalize.dgCMatrix(mat(object), ...), .Cellinfo=object@.Cellinfo))

#' Log transform a matrix for display
#'
#' @param mat matrix of umi counts with genes in rows and cells in columns
#' @param scale_factor factor to multiply the matrix by
#'
#' @return log2 transformed matrix
exprLog.dgCMatrix <- function(mat, scale.factor=7) log2(1+scale.factor*mat)

#' Transform a matrix for display
#'
#' @param mat matrix of umi counts with genes in rows and cells in columns
#' @param smooth_n if a neighborhood smooth is desired; set this to > 1. This only makes sense if the cells are ordered by similarity.
#' @param umi_scale the initial scaling factor after normalizing for UMIs.
#'
#' @export
exprTransform.dgCMatrix <- function(object, smooth_n = 1, umi_scale=1000) {
  mat <- exprNormalize.dgCMatrix(mat, umi_scale)

  # convert to log scale
  lus_1 = exprLog.dgCMatrix(mat)
  #lus_1 = mat

  # move all per-cell expression values to their value above the per-gene median
  #lus = apply(lus_1 - apply(lus_1, 1, median), 2, function(x) pmax(x, 0))
  lus = lus_1

  # if smoothing, smooth per-gene by the requested neighborhood size
  if(smooth_n > 1) {
    apply(lus, 1, function(x) rollmean(x, smooth_n, fill=0))
  } else {
    lus
  }
}
setMethod("exprTransform", "dgCMatrix", exprTransform.dgCMatrix)
setMethod("exprTransform", "scExpressionSet", function(object) {
  object@.Mat <- exprTransform.dgCMatrix(mat(object))
  object
}
)

#' Plot summarizing expression per celltype or metacell
#'
#' You must provide either block_names or a cluster to plot.
#'
#' @param object a list of synBlockExprSet objects
#' @param cluster a cluster to plot, mutally exclusive with block_names
#' @param block_names names of the blocks to plot
#' @param level level to summarize. Should be one of the metadata columns in cell info
#' @param species_order when the cluster is specified, ensure that the order of the species will be as provided
#' @param normalize normalize genewise either not at all ("none"), by the maximum genewise value ("max") or
#'                  by the total genewise expression ("sum")
#' @param cols colors to plot
#' @param rm.unnamed remove cell types containing the word "Unnamed"
#' @param output a filename to output to
#' @param out.width the width of the output plot
#' @param out.height the height of the output plot
#' @param show.colorbar plot the color bar
#'
#' @return a ggplot instance
summaryGridPlot.synBlockExprSet <- function(object,  genomes, block_names=NULL, cluster=NULL,
                                       orientations=NULL, level='cellstate', species_order=NULL,
                                       normalize=c("none","max","sum","scale"),
                                       rm.unnamed=T, output=NULL, show.colorbar=T,
                                       out.width=12, out.height=4, bar.width=5,
                                       cols=if(normalize == "scale") c("red","orange","white","purple","purple4")
                                                                else c("white","orange","red","purple","purple4"),
                                       return_data_frames=F
                                       ) {
  if(!xor(is.null(block_names) , is.null(cluster))) {
    stop("provide either a block_names or a cluster parameter to summaryGridPlot.synBlockExprSet")
  }

  normalize <- match.arg(normalize)

  allblocks <- object %>% map(getBlockExprs.synBlockExprSet) %>% unlist
  names(allblocks) <- allblocks %>% map(getBlock.synBlockExpr) %>% map(name)
  if(!is.null(cluster)) {
    if(!is.null(species_order)) {
      blocks <- cluster %>% getBlocks()
      blocks <- blocks[map(blocks, species) %in% species_order]
      block_names <- names(blocks)[order(match(map(blocks, species), species_order))]
    } else {
      block_names <- cluster %>% getBlocks() %>% names()
    }
  }
  blocks <- allblocks[block_names]
  if(is.null(orientations)) {
    orientations <- rep(1, length(block_names))
  }
  stopifnot(length(block_names) == length(orientations))

  dfs <- sapply(1:length(blocks), function(i) {
    sbe <- blocks[[i]]
    e <- getExpr.synBlockExpr(sbe)
    b <- getBlock.synBlockExpr(sbe)
    ci <- if(rm.unnamed == T)
      droplevels(cellinfo(e)[grep('Unnamed', cellinfo(e)$cellstate, invert=T),])
    else
      cellinfo(e)
    df <- do.call(rbind, sapply(split(ci, ci[,level]), function(ct) {
      mat <- as.matrix(mat(e))[,as.character(ct[,'cell'])]

      safe.mean <- function(x) { if(sum(x > 0) > 0) mean(x, na.rm=T) else 0 }
      data.frame(avg.expr=apply(mat,1,safe.mean),
                 pct.expr=apply(mat,1,function(x) sum(x > 0)/length(x)),
                 CellType=ct[,level][[1]], Gene=geneNames(b), GeneId=genes(b)) #%>%
    }, simplify = F))
    if(normalize != "none") {
      if(normalize == "scale") {
        df <- df %>%
          group_by(Gene) %>%
          dplyr::mutate(.,avg.expr=scale(avg.expr))
      } else {
        z <- ifelse(normalize == "sum", sum, max)
        df <- df %>%
          group_by(Gene) %>%
          dplyr::mutate(.,avg.expr=avg.expr/z(avg.expr))
      }
    }
    df$species <- species(getBlock.synBlockExpr(sbe))
    df
  }, simplify=F)

  .round <- ifelse(normalize == "none", ceiling, identity)
  ulim <- max(sapply(dfs, function(df) .round(max(df$avg.expr, na.rm=T))))
  llim <- min(sapply(dfs, function(df) .round(min(df$avg.expr, na.rm=T))))

  plots <- sapply(1:length(dfs), function(i) {
    df <- dfs[[i]]
    sbe <- blocks[[i]]
    b <- getBlock.synBlockExpr(sbe)
    sp <- df$species[1]
    g <- if(orientations[i] == -1) rev(names(genomes[[sp]])) else names(genomes[[sp]])
    df$GeneId <- factor(df$GeneId, g[g %in% levels(df$GeneId)])
    #levels(df$GeneId) <- genes(b)
    df$CellType <- factor(df$CellType, rev(levels(df$CellType)))

    map.geneid <- function(gid)
      unlist(setNames(geneNames(b), unlist(genes(b)))[gid])
    pl <- min(df$pct.expr)
    ph <- max(df$pct.expr)
    if(normalize == "scale") {
      lim <- max(-llim, ulim)
      range <- c(-lim, lim)
      breaks <- c(-lim,0,lim)
    } else {
      range <- c(0, ulim)
      breaks <- range
    }
    color_bar <- if(!show.colorbar) 'none' else guide_colorbar(
        title=ifelse(normalize == "none", 'Mean Log Expr.',
                     ifelse(normalize == "scale", "Scaled Expression", "Gene Norm. Expr.")),
        title.position='top', title.theme = element_text(size=10, hjust=0.5),
        title.hjust = 0.5, ticks=F, barheight = 0.5, barwidth=bar.width, frame.colour='black',
        frame.linewidth = 1, direction='horizontal', label.theme = element_text(size=9),
        order = 1, label = (normalize == "none" || normalize == "scale"))
    ggplot(df, aes(y=CellType, x=GeneId, size=pct.expr, fill=avg.expr)) +
      geom_raster() +
      scale_fill_gradientn(colours=cols, limits=range, breaks=breaks,
                                           labels= function(x) round(x,2), guide="none") +
      scale_color_gradientn(colours=cols, limits=range, breaks=breaks,
                                           labels= function(x) round(x,2), guide="none") +
      scale_size_continuous(breaks=round(c(ceiling(100*pl)/100,(pl+ph)/2,floor(100*ph)/100), digits=2)) +
      theme_cowplot() +
      theme(axis.text.x=element_text(angle=90),
            axis.title = element_blank(), legend.position = 'top',
            legend.justification = 0.5,
            plot.caption = element_text(hjust=.5, face='bold')) +
      scale_x_discrete(labels=map.geneid) +
      labs(caption=block_names[i]) +
      guides(color=color_bar,
             fill='none',
             size=guide_legend(title='Fract. Cells',
                               title.position='top',
                               title.theme = element_text(size=10, hjust=0.5),
                               label.theme = element_text(size=9), order=2,
                               label.position='bottom')
      )
  }, simplify=F)

  p <- plot_grid(plotlist=plots, nrow=1)
  if(!is.null(output)) {
    ggsave(p, file=output, width=out.width, height=out.height)
  }
  if(return_data_frames) {
    dfs
  } else {
    p
  }
}
setMethod("summaryGridPlot", "list", summaryGridPlot.synBlockExprSet)

#' Plot a circle plot summarizing expression per celltype or metacell
#'
#' You must provide either block_names or a cluster to plot.
#'
#' @param object a list of synBlockExprSet objects
#' @param cluster a cluster to plot, mutally exclusive with block_names
#' @param block_names names of the blocks to plot
#' @param level level to summarize. Should be one of the metadata columns in cell info
#' @param species_order when the cluster is specified, ensure that the order of the species will be as provided
#' @param normalize normalize genewise either not at all ("none"), by the maximum genewise value ("max") or
#'                  by the total genewise expression ("sum")
#' @param cols colors to plot
#' @param rm.unnamed remove cell types containing the word "Unnamed"
#' @param output a filename to output to
#' @param out.width the width of the output plot
#' @param out.height the height of the output plot
#' @param show.colorbar plot the color bar
#' @param fill.box plot a full raster with an outline of the circle in the center (experimental)
#'
#' @return a ggplot instance
circlePlot.synBlockExprSet <- function(object,  genomes, block_names=NULL, cluster=NULL,
                                       orientations=NULL, level='cellstate', species_order=NULL,
                                       normalize=c("none","max","sum","scale"),
                                       rm.unnamed=T, output=NULL, show.colorbar=T, fill.box=F,
                                       out.width=12, out.height=4, bar.width=5,
                                       cols=if(normalize == "scale") c("red","orange","white","purple","purple4")
                                                                else c("white","orange","red","purple","purple4")
                                       ) {
  if(!xor(is.null(block_names) , is.null(cluster))) {
    stop("provide either a block_names or a cluster parameter to circlePlot.synBlockExprSet")
  }

  normalize <- match.arg(normalize)

  allblocks <- object %>% map(getBlockExprs.synBlockExprSet) %>% unlist
  names(allblocks) <- allblocks %>% map(getBlock.synBlockExpr) %>% map(name)
  if(!is.null(cluster)) {
    if(!is.null(species_order)) {
      blocks <- cluster %>% getBlocks()
      blocks <- blocks[map(blocks, species) %in% species_order]
      block_names <- names(blocks)[order(match(map(blocks, species), species_order))]
    } else {
      block_names <- cluster %>% getBlocks() %>% names()
    }
  }
  blocks <- allblocks[block_names]
  if(is.null(orientations)) {
    orientations <- rep(1, length(block_names))
  }
  stopifnot(length(block_names) == length(orientations))

  dfs <- sapply(1:length(blocks), function(i) {
    sbe <- blocks[[i]]
    e <- getExpr.synBlockExpr(sbe)
    b <- getBlock.synBlockExpr(sbe)
    ci <- if(rm.unnamed == T)
      droplevels(cellinfo(e)[grep('Unnamed', cellinfo(e)$cellstate, invert=T),])
    else
      cellinfo(e)
    df <- do.call(rbind, sapply(split(ci, ci[,level]), function(ct) {
      mat <- as.matrix(mat(e))[,as.character(ct[,'cell'])]

      safe.mean <- function(x) { if(sum(x > 0) > 0) mean(x[x > 0], na.rm=T) else 0 }
      data.frame(avg.expr=apply(mat,1,safe.mean),
                 pct.expr=apply(mat,1,function(x) sum(x > 0)/length(x)),
                 CellType=ct[,level][[1]], Gene=geneNames(b), GeneId=genes(b)) #%>%
    }, simplify = F))
    if(normalize != "none") {
      if(normalize == "scale") {
        df <- df %>%
          group_by(Gene) %>%
          dplyr::mutate(.,avg.expr=scale(avg.expr, center = T))
      } else {
        z <- ifelse(normalize == "sum", sum, max)
        df <- df %>%
          group_by(Gene) %>%
          dplyr::mutate(.,avg.expr=avg.expr/z(avg.expr))
      }
    }
    df$species <- species(getBlock.synBlockExpr(sbe))
    df
  }, simplify=F)

  .round <- ifelse(normalize == "none", ceiling, identity)
  ulim <- max(sapply(dfs, function(df) .round(max(df$avg.expr, na.rm=T))))
  llim <- min(sapply(dfs, function(df) .round(min(df$avg.expr, na.rm=T))))

  plots <- sapply(1:length(dfs), function(i) {
    df <- dfs[[i]]
    sbe <- blocks[[i]]
    b <- getBlock.synBlockExpr(sbe)
    sp <- df$species[1]
    g <- if(orientations[i] == -1) rev(names(genomes[[sp]])) else names(genomes[[sp]])
    df$GeneId <- factor(df$GeneId, g[g %in% levels(df$GeneId)])
    #levels(df$GeneId) <- genes(b)
    df$CellType <- factor(df$CellType, rev(levels(df$CellType)))

    map.geneid <- function(gid)
      unlist(setNames(geneNames(b), unlist(genes(b)))[gid])
    pl <- min(df$pct.expr)
    ph <- max(df$pct.expr)
    if(normalize == "scale") {
      lim <- max(-llim, ulim)
      range <- c(-lim, lim)
      breaks <- c(-lim,0,lim)
    } else {
      range <- c(0, ulim)
      breaks <- range
    }
    color_bar <- if(!show.colorbar) 'none' else guide_colorbar(
        title=ifelse(normalize == "none", 'Mean Log Expr.',
                     ifelse(normalize == "scale", "Scaled Expression", "Gene Norm. Expr.")),
        title.position='top', title.theme = element_text(size=10, hjust=0.5),
        title.hjust = 0.5, ticks=F, barheight = 0.5, barwidth=bar.width, frame.colour='black',
        frame.linewidth = 1, direction='horizontal', label.theme = element_text(size=9),
        order = 1, label = (normalize == "none" || normalize == "scale"))
    dots <- if(fill.box) { 
      ggplot(df, aes(y=CellType, x=GeneId, size=pct.expr, fill=avg.expr)) +
        geom_raster() +
        geom_point(color='black', stroke=1) +
        geom_point(aes(size=.9*pct.expr,fill=avg.expr, color=avg.expr)) 
    } else {
      ggplot(df, aes(y=CellType, x=GeneId, size=pct.expr, fill=avg.expr)) +
        geom_point(aes(size=pct.expr,fill=avg.expr, color=avg.expr))  
    }
    
    dots +
      scale_fill_gradientn(colours=cols, limits=range, breaks=breaks,
                                           labels= function(x) round(x,2), guide="none") +
      scale_color_gradientn(colours=cols, limits=range, breaks=breaks,
                                           labels= function(x) round(x,2), guide="none") +
      scale_size_continuous(breaks=round(c(ceiling(100*pl)/100,(pl+ph)/2,floor(100*ph)/100), digits=2)) +
      theme_cowplot() +
      theme(axis.text.x=element_text(angle=90),
            axis.title = element_blank(), legend.position = 'top',
            legend.justification = 0.5,
            plot.caption = element_text(hjust=.5, face='bold')) +
      scale_x_discrete(labels=map.geneid) +
      labs(caption=block_names[i]) +
      guides(color=color_bar,
             fill='none',
             size=guide_legend(title='Fract. Cells',
                               title.position='top',
                               title.theme = element_text(size=10, hjust=0.5),
                               label.theme = element_text(size=9), order=2,
                               label.position='bottom')
      )
  }, simplify=F)

  p <- plot_grid(plotlist=plots, nrow=1)
  if(!is.null(output)) {
    ggsave(p, file=output, width=out.width, height=out.height)
  }
  else {
    p
  }
}

setMethod("circlePlot", "list", circlePlot.synBlockExprSet)
