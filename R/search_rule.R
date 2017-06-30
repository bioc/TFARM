#' Finds an association rule in a dataset.
#'
#' The function looks for an association rule in a dataset \code{data}
#' with a fixed left-hand-side (\code{LHS}) and a fixed right-hand-side
#' (\code{RHS}). The function is used in the \code{rulesTF0}.
#' @param data a binary matrix or data.frame in which the rule is searched
#' @param LHS a string with the left-hand-side of the searched rule
#' @param RHS a string with the right-hand-side of the searched rule
#'
#' @return A vector with five elements: the left-hand-side of the searched rule,
#' the right-hand-side of the searched rule, the support, confidence and lift
#' of the found rule.
#' @export
#' @import arules

search_rule <- function(data, LHS, RHS){

    m <- dim(data)[2]
    data.f <- data
    # Transactions extraction
    for (i in 1:m) {
        data.f[,i] <- as.factor(data.f[,i])
    }
    # Rules extraction
    rule <- apriori(data.f, parameter=list(supp=0, conf=0, target="rules"),
                    appearance=list(lhs=items(LHS), rhs=RHS, default="none"))
    if (length(rule) ==  0) {
        print("Rule not found")
        r <- NA
    }
    else if (length(rule) > 0) {
        r_0 = data.frame(lhs = labels(lhs(rule)),
                         rhs = labels(rhs(rule)), rule@quality)
        r <- r_0[r_0$lhs==LHS,]
        if (dim(r)[1] == 0) {
            print("Rule not found")
            r <- NA
        }
    }
    return(r)
}

