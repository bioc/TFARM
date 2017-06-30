#' Builds the itemset in the left-hand-side of an association rule.
#'
#' The function is used in the construction of the left-hand-side of association
#' rules to search for a rule of interest.
#'

#' @param items a string vector.
#' @return Itemset of the left-hand-side of a searched rule, with the items
#' in \code{items}.
#' @export
#'
#' @examples
#' itemset(c("TAF1=1","EP300=1"))
#' # the output is: "TAF1=1,EP300=1"


itemset <- function(items){
    itemset.0 <- paste("{",items[1],sep="")
    if (length(items) > 1){
        for (i in 2:length(items)){
            itemset.0 <- paste(itemset.0,items[i],sep=",")
            }
    }
	itemset.1 <- paste(itemset.0,"}",sep="")
	return(itemset.1)
}
