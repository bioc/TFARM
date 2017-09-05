#' Represents an example of rulesTF output, i.e. the subset of rules whose
#' left-hand-sides contain FOSL2, and the correspondent quality measures.
#'
#' Within the \code{data_man} data collection, the dataset r_FOSL2 has 5 columns
#' and 28 rows:
#' the first column contains the left-hand-side of the rules (\code{lhs}),
#' the second column is the right-hand-side of the rules (\code{rhs}),
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
#' # r_FOSL2 is found in the data_man collection of datasets:
#' data("data_man")
#' head(r_FOSL2$lhs)
"r_FOSL2"
