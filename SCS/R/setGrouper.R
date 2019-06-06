.make_setGrouper_key <- function(sets) {
  paste(sort(sets), collapse='&')
}

add.setGrouper <- function(object, sets, item) {
  key = .make_setGrouper_key(sets)
  if(!exists(key, where=object@.SetMap)) {
    object@.SetMap[[key]] <- list()
  }
  object@.SetMap[[key]] <- c(object@.SetMap[[key]],item)
  object
}

#' @export
setMethod("add", signature(object="setGrouper", sets="character", item="ANY"), add.setGrouper)


values.setGrouper <- function(x, ...) {sapply(ls(x@.SetMap), function(k) x@.SetMap[[k]])}
setMethod("values", "setGrouper", values.setGrouper)

getByGroup.setGrouper <- function(object, sets) {
  key = .make_setGrouper_key(sets)
  if(nchar(key) == 0) list()
  else object@.SetMap[[key]]
}
#' @export
setMethod("getByGroup", "setGrouper", getByGroup.setGrouper)

#' @export
setMethod("show", signature(object="setGrouper"),
          function(object) {
            cat("setGrouper with ", length(object@.SetMap), " subsets\n")
            invisible(sapply(ls(object@.SetMap),
                             function(key) cat(key, " ", length(object@.SetMap[[key]]), "\n")))
          })


#' @export
setMethod("getGroups", signature(object="setGrouper"), function(object) strsplit(ls(object@.SetMap), '&'))

update.setGrouper <- function(object,other)  {
  keys <- union(ls(object@.SetMap), ls(other@.SetMap))
  .combine <- function(o, key) { o@.SetMap[[key]] <- c(object@.SetMap[[key]], other@.SetMap[[key]]); o }
  Reduce(.combine, keys, new("setGrouper"))
}
#' @export
setMethod("update", "setGrouper", update.setGrouper)


#' Make a data indicating the size of each intersection
#'
#' These data frames are suitable for e.g. UpSetR
#'
#' @param object a setGrouper object
#'
#' @export
intersection.data.frame.setGrouper <- function(object) {
  all.groups <- Reduce(union, getGroups(object))
  all.intersections <- setNames(expand.grid(rep(list(0:1), length(all.groups))), all.groups)
  as.data.frame(do.call(rbind, apply(all.intersections, 1, function(r) {
    group <- getByGroup(object, colnames(all.intersections)[which(r == 1)])
    nmemb <- length(group)
    if(nmemb == 0) NULL
    else do.call(rbind,lapply(seq(nmemb), function(n) r))
  })))
}
setMethod("intersection.data.frame", "setGrouper", function(object) intersection.data.frame.setGrouper(object))

