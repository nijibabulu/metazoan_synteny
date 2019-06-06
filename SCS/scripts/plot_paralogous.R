#! /usr/bin/env Rscript

library(tidyverse)
library(ggthemes)
library(cowplot)
library(argparse)
library(devtools)
load_all('.')

handle_args <- function() {
  parser <- ArgumentParser(description=paste("Plot the paralogous fraction of each block",
                                             "as a fraction of all blocks in the species.",
                                             "The initialization scripts loading the data",
                                             "are expected to have been run."))
  parser$add_argument("PARALOGS_FILE", 
                      help=paste("a tab-delimited file represeting one block per line,",
                                 "with the first column being the category (e.g. species),",
                                 "the second column the name of the synteny block",
                                 "and the second representing the paralogous fraction of",
                                 "that block."))
  parser$add_argument("OUTPUT_FILE", help="filename to output the image to")
  parser$parse_args()
}

main <- function() {
  args <- handle_args()
  longnames <- load_data_set('longnames')
  clusters <- load_data_set('clusters')
  block_names <- names(getBlocks(clusters))
  paralog.data <- read.table(args$PARALOGS_FILE, 
                             col.names=c('Species','Name','Fraction')) %>%
    filter(Name %in% block_names)

  paralog.data$Species <- fct_relevel(paralog.data$Species, c('aq','ml','ta','nv','sm'))
  levels(paralog.data$Species) <- unlist(longnames[levels(paralog.data$Species)])

  categorical.paralogs <- paralog.data %>% mutate(Category=case_when(
    Fraction==0 ~ "Paralogs Only",
    Fraction==1 ~ "Orthologs Only",
    TRUE ~ "Mixed"
  )) %>% mutate(Category=fct_relevel(Category,  "Orthologs Only", "Mixed","Paralogs Only")) %>%
    as_tibble()
  totals <- categorical.paralogs %>% group_by(Species) %>% summarise(BlockCount=n())
  ggplot(categorical.paralogs, aes(x=Species, fill=Category)) + geom_bar(position='stack') +
    theme_cowplot() +
    scale_fill_colorblind(guide=guide_legend(direction = "vertical")) +
    geom_text(data=totals, aes(label=BlockCount, x=Species,y=BlockCount), nudge_y=8, inherit.aes=F)  +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          legend.position='bottom',
          legend.justification = 'center',
          legend.title = element_blank()
          ) + ylab('Count')
  ggsave(args$OUTPUT_FILE, height=5,width=3)
}

if(sys.nframe() == 0) {
  main()
}
