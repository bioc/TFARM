#' Contains genomic regions in the first chromosome of the MCF-7 human breast 
#' adenocarcinoma cell line at the ranges side, and the presence indexes of 
#' transcription factors in such regions at the metadata side.
#'
#' MCF7_chr1 is a Large GRanges in which metadata columns identify
#' transcription factors and genomic coordinates of regions in the first 
#' chromosome of the MCF-7 human breast adenocarcinoma cell line are represented
#' in the left-hand side of the GRanges, therefore each row is a different 
#' genomic region.
#'
#' @docType data
#'
#' @usage data("MCF7_chr1")
#'
#' @format An object of class \code{"GRanges"}
#'
#' @keywords datasets
#'
#' @examples
#' data("MCF7_chr1")
#' head(MCF7_chr1)
"MCF7_chr1"
