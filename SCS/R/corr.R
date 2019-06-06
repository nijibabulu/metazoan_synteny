expr_corr <- function(expr) {
  texpr <-Matrix::t(expr)
  form <- reformulate(termlabels=c(colnames(texpr)[2:ncol(texpr)], '-1'), response=colnames(texpr)[1])
  fit <- lm(form, data=as.data.frame(as.matrix(texpr)))
  sumfit <- summary(fit)
  sumfit$adj.r.squared
}

expr_glm_pois_corr <- function(expr) {
  texpr <-Matrix::t(expr)
  texpr.df <- as.data.frame(as.matrix(texpr))
  attach(texpr.df)
  form <- reformulate(termlabels=c(colnames(texpr.df)[2:ncol(texpr)], '-1'), response=colnames(texpr.df)[1])
  fit <- glm(form, family=poisson)
  rsquared <- rsq(fit)
  detach(texpr.df)
  rsquared
}


expr_glm_nbinom_corr <- function(expr) {
  texpr <-Matrix::t(expr)
  texpr.df <- as.data.frame(as.matrix(texpr))
  attach(texpr.df)
  form <- reformulate(termlabels=c(colnames(texpr.df)[2:ncol(texpr)], '-1'), response=colnames(texpr.df)[1])
  fit <- glm.nb(form)
  rsquared <- rsq(fit)
  detach(texpr.df)
  rsquared
}

expr_maxcorr <- function(expr) {
  texpr <- Matrix::t(expr)
  corr <- cor(as.matrix(texpr))
  max(corr[row(corr) != col(corr)])
}



fisherz <- function(m) {
  0.5 * log((1 + m)/(1 - m))
}

fisherz.inv <- function (z)
{
  (exp(2 * z) - 1)/(1 + exp(2 * z))
}

# take a correlation matrix and calculate the fisher-z corrected mean
fischer_transformed <- function(m, func) {
  z <- fisherz(m)
  mu <- func(m[row(m) != col(m)])
  fisherz.inv(mu)
}

mean_corr <- function(m) fischer_transformed(m, mean)
var_corr <- function(m) fischer_transformed(m, var)

expr_varcorr <- function(expr) {
  texpr <-Matrix::t(expr)
  corr <- suppressWarnings(cor(as.matrix(texpr)))
  var_corr(corr)
}

expr_meancorr <- function(expr) {
  texpr <-Matrix::t(expr)
  corr <- suppressWarnings(cor(as.matrix(texpr)))
  mean_corr(corr)
}

expr_meancorr_spearman <- function(expr) {
  texpr <-Matrix::t(expr)
  corr <- suppressWarnings(cor(as.matrix(texpr), method='spearman'))
  mean_corr(corr)
}

expr_max_bicorr <- function(expr) {
  texpr <-Matrix::t(expr)
  corr <- suppressWarnings(bicor(texpr))
  max(corr[row(corr) != col(corr)])
}

expr_mean_bicorr <- function(expr) {
  texpr <-Matrix::t(expr)
  corr <- suppressWarnings(bicor(texpr))
  mean_corr(corr)
}

expr_var_bicorr <- function(expr) {
  texpr <-Matrix::t(expr)
  corr <- suppressWarnings(bicor(texpr))
  var_corr(corr)
}

