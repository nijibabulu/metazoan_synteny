#! /usr/bin/env Rscript

library(devtools)
library(ggplot2)
library(cowplot)
load_all()

plot.fit <- function(data, fit, name, end=4) {
  fit_values <- sapply(0:end, function(x) dnbinom(x, fit$estimate[1], fit$estimate[2]))
  
  ggplot(data.frame(InterveningSize=unlist(data)[which(data <= 4)]), aes(InterveningSize, stat(density))) + 
    geom_histogram(bins=end+1) + lims(y=c(0,1)) +
    geom_point(data=data.frame(x=0:end,y=fit_values), aes(x=x,y=y),color='red') + ggtitle(name) +
    annotate('text', x=0, y=1, vjust=1, hjust=0, size=4.5,
             label=paste('r =', round(fit$estimate[1],2), '\np =', round(fit$estimate[2],2)))
}

intervening_sizes_data <- load_data_set("intervening_sizes_data")
fits <- load_data_set("fits")
plots <- mapply(plot.fit, intervening_sizes_data, fits, longnames, SIMPLIFY = F) 
plot_grid(plotlist=plots[c('aq','ml','ta','nv','sm')], nrow=1)
