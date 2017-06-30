#' Contains the delta variations of support,
#' confidence and lift.
#'
#' DELTA is a list of 12 elements and each element has three columns
#' representing support (\code{diff_supp_Z}), confidence (\code{diff_conf_Z})
#' and lift (\code{diff_lift_Z}) respectively. It is included in the 
#' \code{data_man} collection.
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
#' # DELTA is found in the data_man collection of datasets:
#' data("data_man")
#' head(DELTA[[1]])
"DELTA"
