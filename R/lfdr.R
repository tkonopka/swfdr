# Computation of local false discovery rate
#


#' Estimation of local false discovery rate
#'
#' Local FDR is a comparison of a p-value distribution with a flat/uniform distribution.
#' The implementation is similar to the code for lfdr in package qvalue.
#' An important difference is that the monotone transformation is carried out
#' before the adjustment by pi0. 
#' 
#' @param p numeric vector of p-values
#' @param pi0 numeric vector of pi0 estimates (same length as p)
#' @param eps numeric, regularization of p-values, same effect as in qvalue::lfdr
#' @param adj numeric, adjustment of smoothing bandwidth in density estimte, same effect
#' as in qvalue::lfdr
#'
#' @return numeric vector of local fdr values, truncated into interval [0, 1]
lfdr <- function(p, pi0, eps=1e-8, adj=1.5) {
  
  # estimation of density of empirical p-values in units of the normal distribution
  p <- pmin(pmax(p, eps), 1-eps)
  x <- qnorm(p)
  xd <- density(x, adjust = adj)
  xs <- smooth.spline(x = xd$x, y = xd$y)
  y <- predict(xs, x)$y
  result <- dnorm(x) / y

  # apply the monotone transformation
  o <- order(p, decreasing = FALSE)
  ro <- order(o)
  result <- cummax(result[o])[ro]
  
  pmin(pi0*result, 1)
}

