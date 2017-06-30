#' Contains the mean Importance Index of pairs of transcription factors which
#' are present in at least one association rule.
#'
#' Within the \code{data_man} data collection, the dataset I_c_2 has 2 columns
#' and 78 rows:
#' the fist column (\code{TF}) contains the transcription factors
#' (single or in pairs) and the second column (\code{imp}) lists the Importance
#' Indexes associated with each trascription factor or pair of transcription
#' factors.
#'
#' @docType data
#'
#' @usage data("data_man")
#'
#' @format An object of class \code{"data.frame"}
#'
#' @keywords datasets
#'
#' @examples
#' # I_c_2 is found in the data_man collection of datasets:
#' data("data_man")
#' head(I_c_2$imp)
"I_c_2"
