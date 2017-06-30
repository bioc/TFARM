#' Contains the candidate co-regulators and the number of rules
#' associated with them.
#'
#' Within the \code{data_man} data collection, the dataset I has 3 columns 
#' and 12 rows:
#' the fist column contains the transcription factors (\code{IMP.TF}),
#' the Importance Indexes associated with each trascription factor
#' are listed in the second column (\code{IMP.imp}) and the third column
#' contains the number of rules found for each transcription factor
#' (\code{IMP.nrules}).
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
#' # I is found in the data_man collection of datasets:
#' data("data_man")
#' head(I$IMP.TF)
"I"
