#! /usr/bin/env Rscript
library(tidyverse)
library(ggthemes)
library(cowplot)
library(devtools)
library(broom)
load_all('.')


#' Modified version of the ggbiplot function from vqv/ggbiplot
#'
#' @param var.labels a list of labels to show from the biplot
#' @param minor.alpha the alpha value to use for the variables that are not
#'   included in var.labels
#' @param label.samples show the names of the samples
#' @param plot.points plot points even if labels are present
#' @param top.var plot the top N var loadings (mutually exclusive to var.labels)
#' @param equal.coords plot such that the scale of x and y are equal
#'
ggbiplot.mod <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
          obs.scale = 1 - scale, var.labels=NULL, var.scale = scale, groups = NULL,
          ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3,
          alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69,
          varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE,
          minor.alpha=.4, minor.color='gray', label.samples=F, plot.points=F,
          top.var=NULL, equal.coords=T,
          ...)
{
  library(ggplot2)
  library(scales)
  library(grid)
  library(ggrepel)
  stopifnot(length(choices) == 2)
  if(!xor(is.null(top.var), is.null(var.labels))) {
    stop("You may only provide one of top.var or var.labels")
  }
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord),
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  #if(circle) {
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
    v.scale <- rowSums(v^2)
    df.v <- r * df.v/sqrt(max(v.scale))
  #}
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% var.)",
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  } else if(label.samples) {
    df.u$labels <- rownames(df.u)
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) +
    ylab(u.axis.labs[2])
  if(equal.coords) {
    g <- g + coord_equal()
  }
  if(!is.null(top.var)) {
    df.v$size <- sqrt(df.v$xvar^2 + df.v$yvar^2)
    var.labels <- rownames(df.v[order(df.v$size, decreasing = T),])[1:top.var]
  }
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi,
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r *
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"),
                         size = 1/2, alpha = 1/3)
    }
    df.v$major = if(is.null(var.labels)) T else df.v$varname %in% var.labels
    #arrow.alpha <- if(is.null(var.labels)) 1 else ifelse(df.v$varname %in% var.labels, 1, minor.alpha)
    #arrow.color <- if(is.null(var.labels)) 1 else ifelse(df.v$varname %in% var.labels, muted('red'), minor.color)
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0,
                                           xend = xvar, yend = yvar, alpha=major, color=major),
                          arrow = arrow(length = unit(1/2,  "picas")),
                          #color = muted("red"),
                          show.legend = F) +
      scale_color_manual(values=c(minor.color, muted('red'))) +
      scale_alpha_manual(values=c(minor.alpha, 1))
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text_repel(aes(label = labels, color = groups),
                               size = labels.size)
    }
    else {
      g <- g + geom_text_repel(aes(label = labels), size = labels.size)
    }
  }
  if(is.null(df.u$labels) || plot.points) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2,
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    df.vl <- if(is.null(var.labels)) df.v else with(df.v, df.v[varname %in% var.labels,])
    g <- g + geom_text(data = df.vl,
                       aes(label = varname,  x = xvar, y = yvar,  hjust = hjust),
                       color = "darkred",  size = varname.size, fontface="bold")
  }
  return(g)
}

plot_tau_correlation <- function(ca, cb, data, name_map) {
  ggplot(blocks_tau_spread, aes_string(x=ca, y=cb)) +
      geom_point(color='red', size=.5) +
      geom_smooth(color='black', size=.5, method=lm, se=TRUE) + lims(x=c(0,1),y=c(0,1)) +
      annotate("text", x=.1, y=.9, hjust=0, parse=T,
               label=lm(as.formula(str_glue("{ca}~{cb}")), data=blocks_tau_spread) %>% glance() %>%
                 str_glue_data("R^{2}=={signif(r.squared,2)}~p=={signif(p.value,2)}")) +
      labs(x=name_map[[str_sub(ca, start=-2)]],
           y=name_map[[str_sub(cb, start=-2)]])
}

main <- function() {
  blocks_tau <- load_data_set('blocks_tau')

  blocks_tau_gathered <- blocks_tau %>%
    gather("tau_type", "tau", tau, tau_unnamed)

  longnames <- load_data_set('longnames')

  clusters.reducedA <- load_data_set('clusters.reducedA')
  common <- query.synClustSet(clusters.reducedA, c('aq','nv','ta'), strict=F)
  common.names <- common %>% getSpeciesBlocks('aq') %>% names()
  common.tbl <- common %>%
    getBlocks %>%
    map(~tibble(name=name(.),
                species=species(.),
                cluster_name=clusterName(.))) %>%
    bind_rows()

  blocks_tau_spread <- blocks_tau_gathered %>% inner_join(common.tbl) %>%
    filter(block_type=='observed') %>%
    select(-name, block_type) %>%
    unite(tau_species, tau_type, species) %>%
    spread(tau_species, tau)

  combn(c('tau_aq','tau_ta','tau_nv'), 2, simplify=F) %>%
    map(~plot_tau_correlation(.[1],.[2], data=blocks_tau_spread, name_map=longnames))  %>%
    plot_grid(plotlist=., nrow=1, labels=letters)
  ggsave(filename='SupplementaryFigure8.pdf', width=9, height=3)

  combn(c('tau_unnamed_aq','tau_unnamed_ta','tau_unnamed_nv'), 2, simplify=F) %>%
    map(~plot_tau_correlation(.[1],.[2], data=blocks_tau_spread, name_map=longnames))  %>%
    plot_grid(plotlist=., nrow=1, labels=letters)
  ggsave(filename='tauplot_w_unnamed.pdf', width=9, height=3)

  blocks_tau_gathered[blocks_tau_gathered$name=='147(aq)',]$name = '*147(aq)'
  blocks_tau_gathered[blocks_tau_gathered$name=='71(aq)',]$name = '*71(aq)'

  blocks_tau_gathered %>%
    filter(species=='aq', block_type=='observed', name %in% c('*147(aq)', '*71(aq)',common.names)) %>%
    mutate(name=fct_reorder(.$name, .$tau)) %>%
    ggplot(aes(x=name,y=tau,fill=tau_type)) + geom_bar(position='dodge',stat='identity') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave('tau_allmeasures.pdf')

  blocks_tau_gathered %>%
    filter(species=='aq', block_type=='observed',
           tau_type=='tau_genewise',
           name %in% c('*147(aq)', '*71(aq)',common.names)) %>%
    mutate(name=fct_reorder(.$name, .$tau)) %>%
    ggplot(aes(x=name,y=tau,fill=tau_type)) + geom_bar(position='dodge',stat='identity') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave('tau_genewise.pdf')

  blocks_tau_gathered %>%
    filter(species=='aq', block_type=='observed',
           tau_type=='tau_genewise_geomean',
           name %in% c('*147(aq)', '*71(aq)',common.names)) %>%
    mutate(name=fct_reorder(.$name, .$tau)) %>%
    ggplot(aes(x=name,y=tau,fill=tau_type)) + geom_bar(position='dodge',stat='identity') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave('tau_genewise_geomean.pdf')


  biases <- load_data_set('scaled_biases')

  biases_common <- biases %>% inner_join(common.tbl)

  biases_split <- biases_common %>% group_split(species) %>%
    set_names(biases_common %>% group_keys(species) %>% flatten_chr())

  pcas <- biases_split %>%
    map(~#filter(., cellstate != "Sperm") %>%
          select(., -species, -cluster_name) %>%
          spread(name, xi) %>%
          column_to_rownames("cellstate") %>%
          prcomp())


  blocks_of_interest <- c('71(aq)', '498(ta)', '70(nv)', '37(aq)', '350(ta)', '38(nv)')


  map(pcas[c('aq','ta','nv')], function(p) {
    margin <- if(identical(p,pcas[['nv']])) 3 else 2
    ggbiplot.mod(p, obs.scale=1, var.scale=0, var.labels = blocks_of_interest,
                 label.samples = T, plot.points=T, equal.coords = F) +
      theme(axis.title = element_text(size=10),
            axis.text = element_text(size=9),
            plot.margin = unit(c(0,margin,0,0), 'lines')) +
      scale_x_continuous(expand=expand_scale(mul=.3)) +
      scale_y_continuous(expand=expand_scale(mul=.3))
  }) %>% plot_grid(plotlist=., nrow=1)
  ggsave('biplot.cur.pdf', width=10, height=10/3)


  # for the final figure, use "B" and "C" rather than the block names
  c_blocks <- c('71(aq)', '498(ta)', '70(nv)')
  b_blocks <- c('37(aq)', '350(ta)', '38(nv)')
  map(pcas[c('aq','ta','nv')], function(p) {
    rownames(p$rotation)[which(rownames(p$rotation) %in% b_blocks)] <- "b"
    rownames(p$rotation)[which(rownames(p$rotation) %in% c_blocks)] <- "c"
    margin <- if(identical(p,pcas[['nv']])) 3 else 2
    ggbiplot.mod(p, obs.scale=1, varname.size=4, var.scale=0, var.labels = c("b","c"),
                 label.samples = T, plot.points=T, equal.coords = F) +
      theme(axis.title = element_text(size=10),
            axis.text = element_text(size=9),
            plot.margin = unit(c(0,margin,0,0), 'lines')) +
      scale_x_continuous(expand=expand_scale(mul=.3)) +
      scale_y_continuous(expand=expand_scale(mul=.3))
  }) %>% plot_grid(plotlist=., nrow=1)
  ggsave('biplot.scaled.longloadings.cur.pdf', width=10, height=10/3)
}

if(sys.nframe() == 0) {
  main()
}
