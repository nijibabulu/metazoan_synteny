#! /usr/bin/env Rscript
library(tidyverse)
library(magrittr)
library(ggthemes)
library(cowplot)
library(ggsignif)
library(devtools)
load_all('.')

load_file <- function(name) {
  e <- new.env()
  load(name, envir = e)
  if(length(e) == 1) {
    return(e[[ls(e)[[1]]]])
  } else {
    return(as.list(e))
  }
}


map_signif_level <- c(`****` = 1e-04, `***` = 0.001,  `**` = 0.01, `*` = 0.05, ns = 1)

#' Plot a grid of boxplots showing significance levels of block correlations among sets of blocks under different imputation regimes and parameters therein
#'
#' @param tbl a tibble with the columns species, method, block_type, param,
#'            name, and cor (the latter being the value of interest).
#' @param key the block type to compare other block types to (observed, as opposed to randomized)
#' @param comparisons a list of comparisons to show significance levels
#' @param bracket_y positions for each of of the comparisons to display the significance brakcets
#'
#' @return a ggplot.
make_grid_plot <- function(tbl,
                      key="observed",
                      comparisons=list("shuffled", "sampled.ortholog", "sampled.all") %>% map(~c(key, .)),
                      ylims=c(0,0.4),
                      bracket_y=NULL) {
  if(is.null(bracket_y)) {
    h=ylims[2]-ylims[1]
    bracket_y=c(.75,.825,.9)*h + ylims[1]
  }

  ggplot(tbl, aes_string(x='block_type', y='cor', fill='block_type')) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_manual(values=ptol_pal()(4)[c(1,3,2,4)],guide=guide_legend(nrow=2)) +
    facet_grid(species~method+param) +
    theme_cowplot() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    #do.call(geom_signif, argz) +
    geom_signif(comparisons=comparisons,
                test="wilcox.test", test.args=list(paired=F,exact=F),na.rm=T,
                map_signif_level = map_signif_level,
                color="black", tip_length = 0, size=.1,
                y_position=bracket_y, vjust=.6, data=NULL) +
    #stat_compare_means(aes(group=species),
    #                   comparisons=comparisons, na.rm=T, label='p.signif',
    #                   label.y=c(.36,.33,.3), tip.length = 0.01, vjust=.5) +
    scale_y_continuous(name="Mean Pairwise Correlation", limits=ylims) +
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1,0,0,0), units='cm'),
          legend.position = 'none',
          #legend.justification = 'center',
          strip.text = element_text(margin=margin(5,0,5,0,'pt')),
          axis.ticks.x = element_blank())
}

#' Plot the distribution of UMI values among blocks as compared to the rest of the matrix
#'
#' @param counts named list of counts matrices
#' @param blocks synClustSet object of synteny blocks whose count list should be highlighted
#' @param species the species from which to draw the counts
plot_umi_distribution <- function(counts, blocks, species) {
  synGenes <- blocks %>% getSpeciesBlocks(species) %>% map(genes) %>% flatten_chr()
  count_dists <- names(counts) %>%
    map(function(name) {
      mat(counts[[name]]) %>% as.matrix() %>%
          as.data.table(keep.rownames='gene') %>%
          melt(id.vars=c("gene"), variable.name="cell", value.name="count") %>%
        .[count > 0,] %>%
        .[,syntenic := gene %in% synGenes] %>%
        .[,matrix := name]
    }) %>% rbindlist()
  ggplot(count_dists, aes(x=count,fill=syntenic)) +
    geom_histogram(alpha=.5,bins=100) +
    scale_fill_manual(values=c("#999999", "#ef8a62")) +
    scale_x_log10() +
    facet_wrap(.~matrix)
}

#' Plot a boxplot showing significance levels of block correlations among sets of blocks
#'
#' @param tbl a tibble with the columns species, method, block_type, param,
#'            name, and cor (the latter being the value of interest).
#' @param key the block type to compare other block types to (observed, as opposed to randomized)
#' @param comparisons a list of comparisons to show significance levels
#' @param bracket_y positions for each of of the comparisons to display the significance brakcets
#'
#' @return a ggplot.
make_plot <- function(tbl,
                      key="observed",
                      comparisons=list("shuffled", "sampled.ortholog", "sampled.all") %>% map(~c(key, .)),
                      bracket_y=NULL,
                      ylims=c(0,0.4)) {

  if(is.null(bracket_y)) {
    h=ylims[2]-ylims[1]
    bracket_y=c(.9,.825,.75)*h + ylims[1]
  }

  size.summary <- tbl %>% filter(block_type==key) %>% group_by(species) %>% summarize(label=paste0("n=",n()))
  ggplot(tbl, aes_string(x='block_type', y='cor', fill='block_type')) +
    geom_boxplot(outlier.shape=NA) + #outlier.size = .1) +
    scale_fill_manual(values=ptol_pal()(4)[c(1,3,2,4)],guide=guide_legend(nrow=2)) +
    facet_grid(~species) +
    theme_cowplot() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    geom_signif(comparisons=comparisons,
                test="wilcox.test", test.args=list(paired=F,exact=F),na.rm=T,
                map_signif_level = map_signif_level,
                color="black", tip_length = 0.01, size=.5, #vjust=.6,
                y_position=bracket_y,  data=NULL) +
    scale_y_continuous(name="Mean Pairwise Correlation", limits=c(0,0.4)) +
    theme(legend.title = element_blank(),
          plot.margin = unit(c(1,0,0,0), units='cm'),
          legend.position = 'bottom',
          legend.justification = 'center',
          strip.text = element_text(margin=margin(5,0,5,0,'pt')),
          axis.ticks.x = element_blank()) +
    geom_text(data=size.summary, aes(x=2.5,y=0.4,hjust=0.5,label=label), inherit.aes=F)
}

#' plot the matrix of synblocks
#'
#' @param scs a synClusterSet object
#' @param counts an scExpressionSet object
#' @param outname a path to save to
#' @param species a species name
#' @param ... arguments to pass to plotMat.synBlockExprSet
#'
#' @return nothing
plot_blocks <- function(scs, counts, outname=NULL, species='aq', ...) {
  sbes <- makeSynBlockExprSet(scs, species='aq', es=counts)
  plotMat.synBlockExprSet(sbes, show.summary = F, ...)
  if(!is.null(outname)) {
    ggsave(outname)
  }
}

#' compute a block correlation assuming the matrix consists only of the genes of interest.
#'
#' @param mat a matrix of genes by cells
#'
#' @return the maximum of 0 and the block correlation
block_correlation <- function(mat, use="everything") {
  corrs <- mat %>% as.matrix() %>% t() %>% cor(method='spearman', use=use) %>% fisherz()
  mean(corrs[row(corrs) != col(corrs)], na.rm=T) %>% fisherz.inv() %>% max(0)
}

#' plot the cellstate-wise correlation of a set of clusters
#'
#' @param scs a synClusterSet object with clusters of interest
#' @param counts a scExpressionSet object
#' @param species species of interest
#' @param outfile file to save to
#'
#' @return nothing
plot_cellwise_corrs <- function(scs, counts, species='aq',  outfile=NULL, scale=F, use="everything") {
  sbes <- makeSynBlockExprSet(scs, species=species, es=counts)
  exprs <- sbes %>% getBlockExprs.synBlockExprSet() %>% map(getExpr.synBlockExpr)
  if(scale == T) {
    exprs %<>% map(~new("scExpressionSet",
                        .Mat=scale(mat(.)) %>% as("dgCMatrix"),
                        .Cellinfo=cellinfo(.)))
    use = 'pairwise.complete.obs'
    print(use)
  }
  cell_groups <- counts %>% cellinfo() %>% select(cellstate, cell) %>% group_by(cellstate) %>% nest()
  cellwise_correlations <- exprs %>% map(
    function(expr)
      pmap(cell_groups,  ~mat(..3)[,..2$cell] %>% block_correlation(use=use), expr)
    ) %>% as_tibble() %>% add_column(cellstate=cell_groups$cellstate) %>% unnest()
  cellwise_correlation_t <- cellwise_correlations %>%
    gather(block, corr, -cellstate) %>%
    mutate(block=fct_relevel(block, names(exprs)))

  ggplot(cellwise_correlation_t, aes(y=block,x=cellstate,fill=corr)) +
    geom_raster() +
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_fill_gradient(low="white", high="red")

  ggsave(outfile)
}

diverging.scale <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

main <- function() {
  # import all imputation data into a single tibble
  impute.cors <- Sys.glob('inst/*.cor.fullimputed.rda') %>% map(load_file) %>% reduce(bind_rows) %>%
    mutate(species=fct_relevel(species, 'aq', 'ml', 'ta','nv', 'sm'))

  # plot the original plot
  impute.cors %>% filter(method=='unimputed') %>% make_plot()
  ggsave('fig1a.pdf')

  # plot the magic correlation values with only even parameter settings
  impute.cors.subset <- impute.cors %>% filter((method=='magic' & param %in% c(5,6,9,17)) |
                           (method=='drimpute' & param %in% c(5,10,15))|
                           (method=='scimpute' & param %in% c(5,10,15)))
  impute.cors.subset[impute.cors.subset$method=='drimpute',] %<>% mutate(param=str_glue("{param}-{param+5}"))
  impute.cors.subset %>% mutate(
    param=fct_reorder(param, str_split(impute.cors.subset$param, '-') %>%
                        map_chr(~.[1]) %>%
                        as.integer()),
    method=fct_relevel(method, c('magic','drimpute','scimpute'))) %>%
    make_grid_plot(ylims=c(0,1.5))
  ggsave('all_summary_corrplotmatrix.pdf', width=12, height=10)
  impute.cors %>% filter(method=='magic' & param %% 2 == 0) %>% make_grid_plot(ylims=c(0,1.5))
  ggsave('magic_corrplotmatrix.pdf', width=12, height=10)


  # plot the dr impute values with odd parameter settings (5-15)
  impute.cors %>% filter(method=='drimpute') %>% mutate(param=str_glue("{param}-{param+5}")) %>% make_grid_plot(ylims=c(0,1.5))
  ggsave('drimpute_corrplotmatrix.pdf', width=12, height=10)

  # plot the scImpute values with odd parameter settings (5-15)
  impute.cors %>% filter(method=='scimpute') %>% make_grid_plot(ylims=c(0,1.5))
  ggsave('scimpute_corrplotmatrix.pdf', width=12, height=10)
}

if(sys.nframe() == 0L) {
    main()
}
