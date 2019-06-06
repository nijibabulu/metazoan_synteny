#! /usr/bin/env Rscript

library(argparse)
library(devtools)

supported.methods <- function() { c("magic", "drimpute", "scimpute", "scimputedrop", "none") }

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
    default.python="/apps/python3/3.7.0/bin/python3"
    default.method=supported.methods()[1]

    sp.names <- load.data('sp.names')[[1]]
    default.sp.names <- paste(sp.names, collapse=',')

    parser <- ArgumentParser(description="Compute the imputed matrix and correlations of the normalized count matrices in the SCS project. Note currently this only computes based on the ess.norm data.")
    parser$add_argument("-m", "--method", metavar="METHOD",
                        choices=supported.methods(), default=default.method,
                        help=paste("Method to use for imputation, one of",
                                   paste(supported.methods(), collapse=', '),
                                   ". Default:", default.method))
    parser$add_argument("-t", "--test", action="store_true",
                        help="Run a test on a small subset of the data")
    parser$add_argument("-d", "--out-dir", type="character", default=baseLoc(),
                        metavar="DIR", help=paste("Output to DIR [default",
                                               baseLoc(),"]"))
    parser$add_argument("-s", "--species", type="character",
                        default=default.sp.names,
                        help="comma-separated list of species for to compute")
    parser$add_argument("-p", "--with-python", type="character", default=default.python,
                        help=paste("which python to use [default", default.python,"]"))
    parser$add_argument("-e", "--expression", type="character", default="ess.norm",
                        help="which expression data to start with [default: ess.norm]")
    parser$add_argument("-b", "--only-blocks", action="store_true",
                        help="Only compute correlation for the blocks of the data from the all_blocks dataset.")
    parser$add_argument("T", type="double", help="value for diffusion parameter")
    parser$add_argument("IDENTIFIER", type="character",
                        help="name of the data set; will be used in the file name and identifier")

    args <- parser$parse_args()
    
    tryCatch( {
      if(args$method == "magic") {
         library(reticulate)
         use_python(args$with_python)
         library(Rmagic)
      } else if(args$method == "drimpute") {
         library(DrImpute)
      } else if(args$method == "scimpute" || args$method == "scimputedrop") {
         library(scImpute)
      }
    },
    error = function(e) {
      stop(paste(e, "Please install the library associated with the requested impute method,", args$method))
    })
    
   
    args$out.path = file.path(args$out_dir, paste0(args$IDENTIFIER, '.rda'))
    if(!try(file.create(args$out.path)))
        stop()

    args.sp.names <- strsplit(args$species, ',')[[1]]
    if(length(args.sp.names) < 1) {
        stop(paste0("no species found in '",args$species,"'"))
    }

    if(any(!args.sp.names %in% sp.names)) {
        stop(paste0("one or more are not available in '",args$species,"'"))
    }
    args$sps <- args.sp.names

    args
}

load.data <- function(...) {
    res <- sapply(list(...), function(d) {
               load_data_set(d)
               get(d)
         }, simplify=F, USE.NAMES=T)
    names(res) <- list(...)
    res

}

magic.impute <- function(mat, t) {
    magic.mat <- sqrt(t(mat))
    magic.res <- magic(magic.mat, genes="all_genes", t=t, verbose=F)
    t(2**magic.res$result)
}

drimpute.impute <- function(mat, k) {
    drimpute.mat <- preprocessSC(mat, min.expressed.gene = 0)
    drimpute.mat.log <- log(drimpute.mat + 1)
    drimpute.imp <- DrImpute(drimpute.mat.log, k=k:(k+5))
    exp(drimpute.imp)
}


scimpute.impute <- function(mat, k) {
    # since there is a solid threshold (0.5), based on the expected counts per million UMIs
    # we have to scale up our normalization, which was done in counts per thousand.
    # otherwise no imputation will take place.
    scimpute.mat <- log10(10^3*mat + 1.01)
    scimpute.res <- scImpute:::imputation_model8(count=scimpute.mat, labeled=FALSE, point=log10(1.01),
                                                 drop_thre=0.5, Kcluster=k,
                                                 out_dir=paste0(tempdir(), '/'),
                                                 ncores=8)
    10^scimpute.res$count_imp - 1.01
}

scimpute.drop.impute <- function(mat, d, k=10) {
  # since there is a solid threshold (0.5), based on the expected counts per million UMIs
  # we have to scale up our normalization, which was done in counts per thousand.
  # otherwise no imputation will take place.
  scimpute.mat <- log10(10^3*mat + 1.01)
  scimpute.res <- scImpute:::imputation_model8(count=scimpute.mat, labeled=FALSE, point=log10(1.01),
                                               drop_thre=d, Kcluster=k,
                                               out_dir=paste0(tempdir(), '/'),
                                               ncores=8)

  (10^scimpute.res$count_imp - 1.01)/10^3
}

impute.counts <- function(ess, sps, cellinfo, param, test=F, method=supported.methods()) {
    method <- match.arg(method)
    if (method == "none") {
        res <- ess[sps]
        if(test) {
            sapply(names(res), function(sp) res[[sp]][,1:300])
        }
        else {
            res
        }
    }
    else {
        sapply(sps, function(sp) {
          message(sp)
          impute.mat <- as.matrix(mat(ess[[sp]]))
          if(test)
              impute.mat <- impute.mat[,1:300]
          impute.res <- if(method == "magic") magic.impute(impute.mat, t=param)
            else if(method == "drimpute") drimpute.impute(impute.mat, k=param)
            else if(method == "scimpute") scimpute.impute(impute.mat, k=param)
            else if(method == "scimputedrop") scimpute.drop.impute(impute.mat, d=param)
          sparse.mat <- as(as.matrix(impute.res), 'dgCMatrix')
          scExpressionSet(.Mat=sparse.mat, .Cellinfo = cellinfo[[sp]])
        })
    }
}

shuffle.counts <- function(clusters, ess, n.bs=100) {
    set.seed(06102018)
    shuffle.randomize(clusters, ess, n.bs=n.bs)
}

gene.correlations.restricted <- function(ess, sps) {
  all_blocks <- load_data_set('all_blocks')
  sapply(sps, function(sp) {
    message(sp)
    blocks <- unlist(all_blocks[[sp]])

    # get pairs of genes for which we are interested in the correlation
    pairs <- sapply(blocks, function(block) {
        combs <- combn(genes(block), 2)
        # convert the gene names to matrix indices for speed of construction
        idxs <- match(combs, rownames(ess[[sp]]))
        # split into a list
        idx.list <- split(idxs, col(combs))
        # sort the indices so the matrix is symmetric
        lapply(idx.list, sort)
    })
    upairs <- do.call(rbind,unique(unlist(pairs, recursive=F)))

    # add the lower right corner of the matrix so that the correlation matrix
    # is the same size
    upairs <- rbind(upairs, c(nrow(ess[[sp]]), nrow(ess[[sp]])))

    texprmat <- t(as.matrix(mat(ess[[sp]])))

    # append the index matrix with a 3rd column of correlations
    corrs <- apply(upairs, 1, function(pair) {
      cor(texprmat[,pair[1]], texprmat[,pair[2]], method='spearman')
    })
    mat <- sparseMatrix(i=upairs[,1], j=upairs[,2], x=corrs, symmetric=T)
    dimnames(mat) <- list(rownames(ess[[sp]]), rownames(ess[[sp]]))
    mat
  })
}

gene.correlations <- function(ess, sps) {
    sapply(sps, function(sp) {
               message(sp)
               print(dim(t(mat(ess[[sp]]))))
               #m <- as.matrix(t(mat(ess[[sp]])))
               #m.r <- apply(m, 2L, rank, na.last='keep')
               #m.r <- m.r-rowMeans(m.r)
               #m.r = m.r / sqrt(rowSums(m.r^2))
               #cr <- tcrossprod(m.r)
               cor(t(as.matrix(mat(ess[[sp]]))), method='spearman')
        }, simplify=F)
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

# > system.time(x <- tpcor(mat[1:10000,]))
# 242.352   4.314 244.427
# > system.time(x <- cor(t(mat[1:10000,])))
# 307.262   2.164 309.456
# > system.type(y <- WGCNA::cor(t(mat[1:10000,]), method='spearman'))
# 374.629   3.751 378.416
# > system.time(x <- cor(t(mat[1:10000,]), method='spearman'))
# 307.787   4.183 311.999

main <- function() {
    preload()
    args <- handle.args()
    ess <- load.data(args$expression)[[1]]
    data <- load.data("longnames", "clusters.reducedA", "cellinfo")

    message("imputing...")
    ess.imputed <- impute.counts(ess, args$sps, data$cellinfo,
                                 args$T, method=args$method, test=args$test)
    # message("shuffling...")
    # shuf.counts <- shuffle.counts(data$clusters.reducedA, ess.magic)
    message("computing correlation...")
    if(args$only_blocks) {
      cormats <- gene.correlations.restricted(ess.imputed, args$sps)
    } else {
      cormats <- gene.correlations(ess.imputed, args$sps)
    }
    res <- list("ess.imputed"=ess.imputed,
                 # "shuf.counts"=shuf.counts,
                 "cormats"=cormats)
    save.as(res, args$IDENTIFIER, args$out.path)
}

#save(list=args$IDENTIFIER, file=out.path, compress=F)
if(sys.nframe() == 0L) {
    main()
}
