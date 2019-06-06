#! /usr/bin/env Rscript

library(devtools)
library(dplyr)
load_all('.')

# constant values used throughout the analysis 
longnames <- list('aq'="A. queenslandica", 
                  "ta"="T. adhaerens",
                  "ml"="M. leidyi", 
                  "nv"="N. vectensis", 
                  "sm"="S. mediterranea" 
)
sp.names <- names(longnames)
save_data_set(sp.names,  overwrite=T)
save_data_set(longnames,  overwrite=T)

# raw single cell data
sc <- list(
  'nv'= read.gsm(file.path(extPath(),'GSM2942625_MARSseq_UMI_Table_Adult_Nematostella.txt')),
  'aq'= read.gsm(file.path(extPath(),'GSM3021561_Amphimedon_queenslandica_adult_UMI_table.txt')),
  'ml'= read.gsm(file.path(extPath(),'GSM3021563_Mnemiopsis_leidyi_UMI_table.txt')),
  'ta'= read.gsm(file.path(extPath(),'GSM3021564_Trichoplax_adhaerens_UMI_table.txt')),
  'sm'= read.gsm(file.path(extPath(),'schmidtea_dge.txt'))
)

save_data_set(sc,  overwrite=T)


# genomes (gene locations on chromosomes) parsing
chrom.files <- sapply(sp.names, function(sp) file.path(extPath(), paste(sp,'.chrom',sep='')))
genomes <- sapply(chrom.files, parse_chrom_file, remove_prefix=T)
save_data_set(genomes, overwrite=T)

# parse the orthologous groups so as to allow categorizing genes with orthology and without
genes.with.orthology <- parse_orthogroups_genes(file.path(extPath(), 'orthogroups.txt'), orthologs.only=T)
save_data_set(genes.with.orthology,  overwrite=T)

# read the annotations of the single cell data on the metacell level, where available
readMC <- function(filename) {  df <- read.table(filename, sep='\t', col.names = c('cell', 'metacell')); rownames(df) <- df$cell; df }
metacells <- list('aq'=readMC(file.path(extPath(), 'GSM3021561_Amphimedon_queenslandica_adult_metacell_definition.txt')),
                  'ta'=readMC(file.path(extPath(), 'GSM3021564_Trichoplax_adhaerens_metacell_definition.txt')),
                  'ml'=readMC(file.path(extPath(), 'GSM3021563_Mnemiopsis_leidyi_metacell_definition.txt')),
                  'nv'=readMC(file.path(extPath(), 'Nematostella_adult_metacell_assignments.nohead'))
                  )

# read the cell type level annotations
readC <- function(filename, ...) read.table(filename, sep='\t', col.names = c('metacell', 'cellstate'))
readCnem <- function(filename, ...) read.table(filename, sep='\t', col.names = c('metacell', 'cellstate','none'))[,c(1,2)]
cells <- list('aq'=readC(file.path(extPath(), 'aq_adult.metacells.txt')),
              'ta'=readC(file.path(extPath(), 'ta.metacells.txt')),
              'ml'=readC(file.path(extPath(), 'ml.metacells.txt')),
              'nv'=readCnem(file.path(extPath(), 'Nematostella_adult_metacell_annotation.nohead'))
             )

cells.metacells <- lapply(names(metacells), function(n) merge(metacells[[n]], cells[[n]]))
names(cells.metacells) <- names(metacells)

# add planaria manually since it does not have the same metacell format
planaria.suerat <- read.table(file.path(extPath(), 'Planaria_Seurat_annot.csv'), sep=',', comment.char = '', quote = '', header=T)
planaria.seurat.reformat <- planaria.suerat[,c('final_Id','final_Id','X..')]
colnames(planaria.seurat.reformat) <- c('metacell','cellstate','cell')
planaria.seurat.reformat$cellstate <- planaria.seurat.reformat$cellstate %>%
  { gsub(' [[:digit:]]+$', '',  .) } %>%
  { gsub(' cells', '', .)  } %>% 
  { gsub(' progenitors$', ' PCs', .) }  %>% 
  { gsub(' neoblast', ' NB', .)} %>% factor()

cells.metacells$sm <- planaria.seurat.reformat

orderMetacells <- function(mctab) mctab[with(mctab, order(cellstate,metacell)),]
ordered.cells.metacells <- lapply(cells.metacells, orderMetacells)
setRows <- function(mctab) { rownames(mctab) <- mctab$cell; mctab }
cellinfo <- lapply(ordered.cells.metacells, setRows)

save_data_set(cellinfo,  overwrite=T)

