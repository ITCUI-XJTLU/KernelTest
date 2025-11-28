#' KernelTest: High-performance kernel-based tests for ChIP-seq data
#'
#' This package implements kernel-based nonparametric tests for detecting
#' differential histone enrichment in ChIP-seq data.
#'
#' @importFrom graphics abline axis hist legend lines mtext par points text
#' @importFrom stats pnorm qnorm quantile
#'
#' @useDynLib KernelTest, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @keywords internal
"_PACKAGE"
