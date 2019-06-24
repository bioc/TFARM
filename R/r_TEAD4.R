#' Contains the association rules for the prediction of the presence
#' of the transcription factor TEAD4 in the considered genomic regions, 
#' i.e., with TEAD4 in the right-hand-side of the association rules. 
#'
#' Within the \code{data_man} data collection, the dataset r_TEAD4 has 5 columns 
#' and 28 rows:
#' the first column contains the left-hand-side of the rules (\code{lhs}),
#' the second column is the right-hand-side of the rules (TEAD4=1) (\code{rhs}),
#' the third column reports the support measures (\code{support}),
#' the fourth column contains the confidence measures (\code{confidence}) and
#' the lift measures are listed in the fifth column (\code{lift}).
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
#' # Load r_TEAD4 from the data_man collection of datasets:
#' data("data_man")
#' head(r_TEAD4$lhs)
"r_TEAD4"
