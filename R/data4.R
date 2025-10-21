#' Example ChIP-seq signal matrix: data4
#'
#' Normalized ChIP-seq signal intensities for the same five genes across
#' 280 genomic bins under *condition 2* (e.g., control group). The gene order and
#' genomic bins are identical to those in \code{data1}, allowing direct row-wise comparison.
#'
#' @format A numeric matrix with 5 rows and 280 columns:
#' \describe{
#'   \item{Rows (1–5)}{Genes, in the same order as in \code{data1}.}
#'   \item{Columns (V3–V282)}{Genomic bins or positions along the gene body, representing normalized ChIP-seq signal intensity.}
#' }
#'
#' @usage data(data4)
#' @examples
#' data(data4)
#' dim(data4)
#' image(data4, main = "ChIP-seq signal (condition 2)")
#'
#' @keywords datasets
#' @docType data
#' @name data4
NULL
