#' Contains the Importance Index associated with each co-regulator
#' which is present in at least one association rule.
#'
#' IMP_Z is a list of 12 elements, containig the Importance Indexes for all
#' transcription factors (\code{p}) in the set of relevant association rules
#' extracted. It is included in the \code{data_man} data collection.
#'
#' @docType data
#'
#' @usage data("data_man")
#'
#' @format An object of class \code{"list"}
#'
#' @keywords datasets
#'
#' @examples
#' # IMP_Z is found in the data_man collection of datasets:
#' data("data_man")
#' head(IMP_Z[[1]])
"IMP_Z"
