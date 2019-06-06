#' @section synPair:
#' methods dealing with pairs of synteny blocks
NULL

as.data.frame.synPair <- function(x, row.names = NULL, optional = FALSE, ...) {
  cbind(
    data.frame(name=x@.Name,
               blocks=paste(sapply(x@.Blocks, function(b) b@.Name), collapse=','),
               species=paste(x@.Species, collapse=','),
               gene.ids=x@.GeneIds,
               gene.names=x@.GeneNames,
               gene.summary=geneSummary(x)),
    x@.Stats
  )
}
setMethod("as.data.frame", "synPair", as.data.frame.synPair)

setMethod("geneSummary", "synPair",
          function(object) {
            sprintf("%s: %s\n%s: %s",
                    object@.Species[[1]], geneSummary(object@.Blocks[[1]]),
                    object@.Species[[2]], geneSummary(object@.Blocks[[2]]))
          })

#' @section synPairSet:
#' methods dealing with sets of synteny pairs
NULL

as.data.frame.synPairSet <- function(x, row.names = NULL, optional = FALSE, ...) {
  #first flatten the list
  pairs <- do.call(c, unlist(values.setGrouper(x@.Pairs)))
  do.call(rbind, lapply(pairs, as.data.frame))
}
setMethod("as.data.frame", "synPairSet", as.data.frame.synPairSet)

setMethod("getPairs", "synPairSet",
          function(object) unlist(values.setGrouper(object@.Pairs)))
