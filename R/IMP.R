#' Contains the mean Importance Index of each co-regulator.
#'
#' Within the \code{data_man} data collection, the dataset IMP has 4 columns 
#' and 12 rows:
#' transcription factors present in at least one association rule are listed
#' in the fist column (\code{TF}), the second column (\code{imp}) contains the
#' means of the Importance Index associated with each trascription factor,
#' the standard deviations for each transcription factor are reported in the
#' third column (\code{sd}) and finally the number of rules found for each
#' transcription factor are transcribed in the fourth column (\code{nrules}).
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
#' # IMP is found in the data_man collection of datasets:
#' data("data_man")
#' head(IMP$imp)
"IMP"
