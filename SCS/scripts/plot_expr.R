library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(devtools)
load_all()



main <- function() {
  genes.with.orthology <- load_data_set('genes.with.orthology')
  ess.norm <- load_data_set('ess.norm')
  sp.names <- load_data_set('sp.names')
  ortho.ess.norm <- map2(ess.norm[sp.names], genes.with.orthology[sp.names],
                         function(es, genes) es[which(rownames(es) %in% genes),])
  synBlockESs.norm <- load_data_set('synBlockESs.norm')

  clusters.reducedA <- load_data_set('clusters.reducedA')

  pAQ <- plotMat.synBlockExprSet(synBlockESs.norm$aq, background.samples = 00,
                                 es=ortho.ess.norm$aq, summary.limits = c(5,10),
                                 show.points=T, show.cloud=F, add.spacers = T,
                                 spacer.size=1, proportional.matrix = T,
                                 y.label.sizes = 7,show.block.ids = T,
                                 output='expr.aq.stretch.1.pdf')
  pML <- plotMat.synBlockExprSet(synBlockESs.norm$ml, background.samples = 100,
                                 es=ortho.ess.norm$ml, summary.limits = c(0,10),
                                 show.points=T, include.unnamed=T, show.cloud=F, add.spacers=T,
                                 spacer.size=1, proportional.matrix = T,
                                 y.label.sizes = 7,show.block.ids = T,
                                 output='expr.ml.stretch.1.pdf')
  pTA <- plotMat.synBlockExprSet(synBlockESs.norm$ta, background.samples = 100,
                                 es=ortho.ess.norm$ta, summary.limits = c(4,10),
                                 show.points=T, show.cloud=F, add.spacers = T,
                                 spacer.size=1, proportional.matrix = T,
                                 y.label.sizes = 7,show.block.ids = T,
                                 output='expr.ta.stretch.1.pdf')
  pNV <- plotMat.synBlockExprSet(synBlockESs.norm$nv, background.samples = 100,
                                 es=ortho.ess.norm$nv, summary.limits = c(0,10),
                                 show.points=T, show.cloud=F, add.spacers = T,
                                 spacer.size=1, proportional.matrix = T,
                                 y.label.sizes = 7,show.block.ids = T,
                                 output='expr.nv.stretch.1.pdf')
  pSM <- plotMat.synBlockExprSet(synBlockESs.norm$sm,  background.samples = 100,
                                 es=ortho.ess.norm$sm, summary.limits = c(0,9), min.cellstate.cells = 100,
                                 cell.label.size=2.5, cell.label.proportion=0.2, cell.label.upper.limit =5,
                                 cell.segement.pointer.color='gray',
                                 show.points=T, show.cloud=F, add.spacers = T,
                                 spacer.size=1, proportional.matrix = T,
                                 y.label.sizes = 7,show.block.ids = T,
                                 output='expr.sm.stretch.1.pdf')
  common <- query(clusters.reducedA, c('aq','nv','ta'), strict=F)
  ess.norm <- load_data_set('ess.norm')
  commonSBES <- makeSynBlockExprSet(common, 'aq', ess.norm$aq)
  plotMat.synBlockExprSet(commonSBES, show.summary=F, proportional.matrix=T, output='expr.common.pdf')
  plotMat.synBlockExprSet(commonSBES, show.summary=F, proportional.matrix=T, output='expr.common.labeled.pdf', show.block.ids = T)
}

if(sys.nframe() == 0) {
  main()
}
