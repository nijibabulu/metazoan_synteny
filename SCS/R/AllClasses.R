
#' A hash table like object for intersection sets
#' @rdname setGrouper
#' @export
setClass("setGrouper", slots=list(.SetMap="environment"))
setMethod("initialize", "setGrouper", function(.Object, ...) {
  .Object <- callNextMethod()
  .Object@.SetMap <- new.env(hash= T, parent=emptyenv())
  .Object
})

#' A simple single cell expression set data structure that handles cell metadata
#' @rdname scExpresssionSet
#' @export
scExpressionSet <- setClass("scExpressionSet",
                            slots=list(.Mat = "dgCMatrix",
                                       .Cellinfo = "data.frame"))
setMethod("initialize", "scExpressionSet", function(.Object, ...) {
  .Object <- callNextMethod()
  cells <- intersect(rownames(.Object@.Cellinfo), colnames(.Object@.Mat))
  .Object@.Mat <- .Object@.Mat[,cells]
  .Object@.Cellinfo <- .Object@.Cellinfo[cells,]
  .Object
})


#' A pair of related synteny blocks
#' @rdname synPair
#' @export
synPair <- setClass("synPair",
                     slots=list(.Name="character",
                                .Blocks="list",
                                .Species="list",
                                .GeneIds="character",
                                .GeneNames="character",
                                .Stats="data.frame"))

#' A set of synPair objects organized by a setGrouper
#' @rdname synPairSet
#' @export
synPair <- setClass("synPairSet",
                    slots=list(.Pairs="setGrouper"))



#' A synteny block data structure
#' @rdname synBlock
#' @export
synBlock <- setClass("synBlock",
                     slots=list(.Name="character",
                                .ClusterName="character",
                                .Species="character",
                                .Conns="list",
                                .ConnSpecies="list",
                                .Location="character",
                                .Length="numeric",
                                .GeneIds="list",
                                .GeneNames="list",
                                .Stats="data.frame"))

setMethod("initialize", "synBlock", function(.Object, ...) {
  .Object <- callNextMethod()
  .Object@.Stats <- data.frame(name=.Object@.Name)
  .Object
})

#' A cluster of synteny blocks related by orthology.
#' @rdname synClust
#' @export
setClass("synClust",
         slots=list(.Name="character",
                    .Blocks="environment",
                    .Species="environment"))
setMethod("initialize", "synClust",
          function(.Object, ...) {
            .Object <- callNextMethod()
            .Object@.Blocks <- new.env(hash= T, parent=emptyenv() )
            .Object@.Species <- new.env(hash= T, parent=emptyenv() )
            .Object
          })

#' A set of synteny clusters
#' @rdname  synClustSet
#' @export
setClass("synClustSet",
         slots=list(.Clusters="environment",
                    .Species="list",
                    .ClustersBySpecies="setGrouper",
                    .Genomes="list",
                    .Stats="data.frame"))
setMethod("initialize", "synClustSet", function(.Object, ...) {
  .Object <- callNextMethod()
  .Object@.Clusters <- new.env(hash= T, parent=emptyenv() )
  .Object@.ClustersBySpecies <- new('setGrouper')
  .Object
})

#' A synteny block data structure decorated with expression data
#' @rdname synBlockExpr
#' @export
synBlockExpr <- setClass("synBlockExpr",
                         slots=list(.Block="synBlock",
                                    .Expr="scExpressionSet"))

#' A set of synteny expression sets
#' @rdname synBlockExprSet
#' @export
synBlockExprSet <- setClass("synBlockExprSet",
                            slots=list(.Species="character",
                                       .BlockExprs="list",
                                       .Corr="matrix",
                                       .HClust="list"))


#.Blocks="list",
#.BlockLengths="integer",
#.BlockNames="character",
#.GroupedGenes="list",
#.Expr="dgCMatrix"))

