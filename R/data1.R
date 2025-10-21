#' Example ChIP-seq signal matrix: data1
#'
#' Normalized ChIP-seq signal intensities for five genes across 280 genomic bins
#' under *condition 1* (e.g., case group). Each row corresponds to one gene, and
#' each column represents a genomic position (bin) along the gene body.
#'
#' @format A numeric matrix with 5 rows and 280 columns:
#' \describe{
#'   \item{Rows (1–5)}{Genes. Each row corresponds to one gene, with the same order as in \code{data4}.}
#'   \item{Columns (V3–V282)}{Genomic bins or positions along the gene body, representing normalized ChIP-seq signal intensity.}
#' }
#'
#' @usage data(data1)
#' @examples
#' data(data1)
#' dim(data1)
#' image(data1, main = "ChIP-seq signal (condition 1)")
#'
#' @keywords datasets
#' @docType data
#' @name data1
NULL
