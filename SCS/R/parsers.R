#' @importFrom GenomicRanges GRanges 
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom Matrix Matrix
NULL

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


#' remove the prefix from a vector of gene names (e.g. nv_v1g199235 becomes v1g199235)
#'
#' @param gene_names a character vector of gene names
#'
#' @return a character vector of gene names with the prefix removed
#'
#' @export
remove_gene_prefixes <- function(gene_names) {
  split_names <- strsplit(as.character(gene_names), '_')
  prefix_free <- lapply(split_names, function(c) c[-1])
  sapply(prefix_free, function(c) paste(c, collapse='_'))
}

parse_block <- function(blockline, cluster, remove_prefix=F) {
  fields <- strsplit(blockline, '\t')[[1]]
  gene_ids <- as.list(strsplit(fields[7], ',')[[1]])
  if(remove_prefix) {
    gene_ids <- remove_gene_prefixes(gene_ids)
  }
  gene_names <- as.list(strsplit(fields[8], ',')[[1]])

  # manipulate the connection string to extract the species names
  split.conns <- strsplit(strsplit(fields[4], "(", fixed=T)[[1]], ")", fixed=T)
  conn.species <- unique(c(fields[2],
                           lapply(split.conns[2:length(split.conns)],
                                function(c) c[1])))

  synBlock(.Name=sprintf('%s(%s)', fields[1], fields[2]),
           .ClusterName=cluster,
           .Species=fields[2],
           .ConnSpecies=conn.species,
           .Conns=as.list(strsplit(fields[4], ',')[[1]]),
           .Location=fields[5],
           .Length=as.numeric(fields[6]),
           .GeneIds=list(gene_ids),
           .GeneNames=list(gene_names))
}

parse_syn_cluster <- function(lines, remove_prefix=F) {
  name = substring(lines[[1]],2)
  clust <- new("synClust", .Name=name, .Blocks=new.env(), .Species=new.env())
  blocks <- sapply(lines[2:length(lines)], function(v) parse_block(v, name, remove_prefix = remove_prefix))

  sapply(blocks, function(b) addBlock.synClust(clust, b))

  clust
}

#' Parse a cluster file into a synClustSet
#'
#' @param filename reordered file with blocks for each cluster
#' @param genomes list of chrom files with gene positions (can be created with sapply(chromfiles, parse_chrom_file))
#' @param remove_prefix remove the genome prefix from the gene name (e.g. nv_v1g199235 becomes v1g199235)
#'
#' @return synClustSet
#'
#' @export
parse_cluster_file <- function(clustfile, genomes, remove_prefix=F) {
  lines <- readLines(clustfile)
  header_lines <- grep('^#', lines)
  if(length(header_lines) == 0) return(make.synClustSet(c(), genomes))
  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
  cluster_groups <- splitAt(lines, header_lines)
  clusters <- sapply(cluster_groups, function(g) parse_syn_cluster(g, remove_prefix=remove_prefix))
  names(genomes) <- sapply(genomes, function(gr) gr[1]$genome)
  make.synClustSet(clusters, genomes)
}

#' Parse a chrom file into a GRanges object.
#'
#' @param filename name of the file
#' @param remove_prefix remove the genome prefix from the gene name ()
#'
#' @return GRanges object
#'
#' @export
parse_chrom_file <- function(filename, remove_prefix=F) {
  tab <- read.table(filename, sep='\t',
                    col.names = c('genome', 'geneName', 'chrom', 'strand', 'start', 'end'))
  if(remove_prefix) {

    tab$geneName <- remove_gene_prefixes(tab$geneName)
  }
  ranges <- IRanges(start=tab$start, end=tab$end, names=tab$geneName)

  # prepare for synteny analysis: sort sequence levels and sort by gene position, ignoring strand
  gs <- GRanges(seqnames=tab$chrom, ranges=ranges, strand=tab$strand, genome=tab$genome)
  gs <- sortSeqlevels(gs)
  gs <- sort(gs, ignore.strand=T)
  gs
}

#' Parse a orthogroups file into lists of genes.
#'
#' The file is assumed to have the format
#' NAME<tab>GROUP_SIZE<tab>GENE_ID
#' and gene ids are prefixed with the species ID and suffixed by species specific gene ID
#'
#'
#' @param filename name of the file
#' @param orthologs.only include only genes which have orthologs (in other species)
#'
#' @return a list of genes, separated by species
parse_orthogroups_genes <- function(filename, orthologs.only=F) {
  .add_gene <- function(l, geneid) {
    species <- strsplit(geneid,'_')[[1]][1]
    gene <- sub('[^_]*_','',geneid)
    l[[species]] <- c(l[[species]], gene)
    l
  }
  .modify_list <- function(l,cvec) {
    # XXX this is slow
    species <- unique(lapply(strsplit(cvec[3:length(cvec)], '_'), function(x) x[1]))
    if(orthologs.only && length(species) < 2) {
      return(l)
    }
    Reduce(.add_gene, cvec[3:length(cvec)], l)
  }
  lines <- readLines(filename)
  Reduce(.modify_list, strsplit(lines, '\t'), list())
}
