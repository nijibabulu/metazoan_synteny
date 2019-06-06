#' @importFrom fastmatch %fin%
NULL

#' @section synBlock:
#' methods dealing with synteny blocks
NULL

#' Get the number of intervening genes between each gene in a synteny block
#'
#' @param object a synBlock objet
#' @param genome the underlying GRanges object
#'
#' @export
intervening_sizes.synBlock <- function(object, genome) {
  ranges <- genome[genes(object)]
  ranges <- sort(ranges, ignore.strand=TRUE)
  sapply(1:(length(ranges)-1),
         function(i) {
           which(names(genome) == names(ranges)[i+1]) -
             which(names(genome) == names(ranges)[i]) - 1
         }
  )
}

setMethod("intervening_sizes", "synBlock", intervening_sizes.synBlock)

#' Find a random synteny block
#'
#' @param genome A GRanges object with genes
#' @param length A number of genes to retrieve
#' @param restrict.genes If not NULL, then only look for blocks containing this character vector of genes.
#'
#' @return a GRanges object with the number of genes of interest
#'
#' @export
random_block <- function(genome, size, restrict.genes=NULL, max.tries=200) {
  for(i in seq(max.tries)) {
    random_gene <- genome[sample.int(length(genome),1)]
    block <- .fill_block(random_gene, genome, size-1)
    if(!is.null(block)) {
      #message(paste('made randomized block after', i, 'tries'))
      return(block)
    }
  }

  stop(paste('Could not find a block of size', size, ' in genome ', genome$genome[1]))
}

.fill_block <- function(block, genome, size) {
  if(size == 0) {
    return(block)
  } else {
    gene <- tail(block, 1)
    next_index <- follow(gene, genome)
    if(is.na(next_index)) {
      return(NULL)
    }
    next_gene <- genome[next_index]
    return(.fill_block(c(block, next_gene), genome, size-1))
  }
}

setMethod("name", "synBlock", function(object) object@.Name)

setMethod("clusterName", "synBlock", function(object) object@.ClusterName)

setMethod("species", "synBlock", function(object) object@.Species)

scaffold.synBlock <- function(object)  strsplit(object@.Location, ':')[[1]][1]

startend.synBlock <- function(object) {
  as.integer(strsplit(strsplit(object@.Location, ':')[[1]][2], '..', fixed=T)[[1]])
}

start.synBlock <- function(object) startend.synBlock(object)[1]
end.synBlock <- function(object) startend.synBlock(object)[2]

setMethod("grange", "synBlock", function(object) {
  GRanges(seqnames=scaffold.synBlock(object),
          ranges=IRanges(start=start.synBlock(object),
                         end=end.synBlock(object)))
})

setMethod("show", "synBlock",
          function(object) { cat(object@.Name, "from cluster", object@.ClusterName,
                                 " (", paste(object@.Species, collapse=" "), ") ",
                                 object@.Location, " ", paste(object@.GeneNames, collapse=" "))})

setMethod("genes", signature(object="synBlock"), function(object) unlist(object@.GeneIds))

setMethod("length", "synBlock", function(x) length(genes(x)))

setMethod("geneSummary", signature(object="synBlock"),
          function(object) {
            sprintf("%s (%d genes)", paste(unique(unlist(object@.GeneNames)), collapse=','), length(object@.GeneNames[[1]]))
          })


calcStat.synBlock <- function(object, statfun, expr, mat.only=T) {
  expr.s <- expr[[object@.Species]]
  if(is.null(expr.s)) { as.numeric(NA) }
  else {
    expr.subset <- sc_expr(object, expr.s)
    if(is.null(dim(expr.subset)) || dim(expr.subset)[1] == 0) {
        as.numeric(NA) 
    } else {
        if(mat.only) {
            expr.subset <- mat(expr.subset)
        }
        statfun(expr.subset)
    }
  }
}


addStat.synBlock <- function(object, expr, func, varName, mat.only=T) {
  setStat(object, varName, calcStat.synBlock(object, func, expr, mat.only))
}


setMethod("setStat", signature(object="synBlock", name="character", value="numeric"),
          function(object, name, value) { object@.Stats[,name] = value; return(object) })

setMethod("getStat", signature(object="synBlock", name="character"),
          function(object, name=NULL) {
            if(is.null(name)) {
              object@.Stats
            } else {
              object@.Stats[,name]
            }
          })

as.data.frame.synBlock <- function(x, row.names = NULL, optional = FALSE, ...) {
  cbind(
    data.frame(name=x@.Name,
               species=x@.Species,
               conn.species=paste(sort(unlist(x@.ConnSpecies)), collapse=','),
               location=x@.Location,
               length=x@.Length,
               gene.ids=paste(unlist(x@.GeneIds), collapse=','),
               gene.names=paste(unlist(x@.GeneNames), collapse=',')),
    x@.Stats[, -which(colnames(x@.Stats) == 'name')]
  )
}
setMethod("as.data.frame", "synBlock", as.data.frame.synBlock)

setMethod("geneNames", signature(object="synBlock"), function(object) unlist(object@.GeneNames))
setMethod("geneNames<-", 'synBlock', function(object, names) {
  if(length(names) != length(object@.GeneNames)) {
    stop(paste(length(names), " != ", length(object@.GeneNames),
               " gene names! Must replace gene names with an equal-lengthed vector."))
  }
  object@.GeneNames <- names
})

randomize.synBlock <- function(object, genomes) {
  genome <- genomes[[object@.Species]]
  gr <- random_block(genome, length(genes(object)))
  start <- min(gr@ranges@start)
  end <- max(gr@ranges@start + gr@ranges@width)
  location <- paste0(as.character(gr@seqnames[1]), ':', start, '..', end)

  synBlock(.Name=sprintf('%srand', object@.Name),
           .ClusterName=object@.ClusterName,
           .Species=object@.Species,
           .Conns=list("None(Random)"),
           .Location=location,
           .Length=end-start,
           .GeneIds=list(names(gr)))
}

sc_expr.matrixLike <- function(object, expressionSet) {
  g <- genes(object)
  to_use <- g[g %in% rownames(expressionSet)]
  expressionSet[to_use,]
}
#' Get UMI counts of all genes in a given synBlock
#'
#' @param object synteny block object
#' @param expressionSet matrix containing all sequences
#'
#' @return matrix of UMI counts
#'
#' @export
setMethod("sc_expr", signature(object="synBlock", expressionSet="dgCMatrix"), sc_expr.matrixLike)
setMethod("sc_expr", signature(object="synBlock", expressionSet="matrix"), sc_expr.matrixLike)
setMethod("sc_expr", signature(object="synBlock", expressionSet="scExpressionSet"),sc_expr.matrixLike)

setMethod("makePair", signature(object="synBlock", other="synBlock"),
          function(object,other) {
            .getstats <- function(stats, suffix) {
              stats <- stats[, -which(colnames(object@.Stats) == 'name')]
              colnames(stats) <- paste0(colnames(stats), suffix)
              stats
            }

            new("synPair",
                .Name=paste(c(object@.Name,other@.Name), collapse=","),
                .Blocks=list(object,other),
                .Species=as.list(c(object@.Species,other@.Species)),
                .GeneIds=paste0(object@.Species,":",paste(object@.GeneIds[[1]], collapse=","),'; ',
                                other@.Species,":",paste(other@.GeneIds[[1]], collapse=",")),
                .GeneNames=paste0(object@.Species,":",paste(object@.GeneNames[[1]], collapse=","),'; ',
                                  other@.Species,":",paste(other@.GeneNames[[1]], collapse=",")),
                .Stats=cbind(.getstats(object@.Stats,".x"), .getstats(other@.Stats,".y"),
                             species.x=object@.Species, species.y=other@.Species)
            )
          })

filterGenesInBlock.synBlock <- function(object, filter.function) {
  gene.inds <- which(sapply(genes(object), filter.function, block=object))
  synBlock(.Name=object@.Name,
           .ClusterName=object@.ClusterName,
           .Species=object@.Species,
           .ConnSpecies=object@.ConnSpecies,
           .Conns=object@.Conns,
           .Location=object@.Location,
           .Length=object@.Length,
           .GeneIds=list(object@.GeneIds[[1]][gene.inds]),
           .GeneNames=list(object@.GeneNames[[1]][gene.inds]),
           .Stats = object@.Stats)
}
setMethod("filterGenesInBlock", "synBlock", filterGenesInBlock.synBlock)

