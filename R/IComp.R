#' Computes the Importance Index of a transcription factor
#' in a set of association rules.
#'
#' Given an association rule and a transcription factor \code{TFi},
#' it is evaluated the contribution of \code{TFi} in the rule for the prediction
#' of the presence of the item in the right-hand-side of the rule.
#' Since this contribution is evaluated based on the variations of
#' support, confidence and lift of the rule, the user can visualize
#' such variations by setting the parameter \code{figures = TRUE}.
#'
#' @param TFi string or string vector: the transcription factor (or combination
#' of transcription factors) whose importance distribution is evaluated.
#' @param rules_TF a set of rules in which \code{TFi} is present.
#' @param rules_noTF a set of rules obtained from rules_TF removing
#' the transcription factor (or combination of transcription factors)
#' in TFi (obtained with the function \code{\link{rulesTF0}}).
#' @param figures logical; if \code{figures = TRUE}, graphics with support,
#' lift and confidence distributions of the rules in \code{rulesTF} and
#' \code{rulesTF0} are returned.
#'
#' @return A list of four elements: the \code{imp} element is a vector of
#' doubles with the importances of TFi in each rule in \code{rulesTF};
#' the \code{delta} element of the list is a list with variations of
#' standardized distributions of the three measures of support,
#' confidence and lift. This output is used in the function \code{\link{IPCA}}
#' for the Principal Component Analysis of such distributions. The \code{rwi}
#' element is a data.frame with the rules in \code{rulesTF}
#' in which the transcription factor \code{TFi} has Importance Index greater
#' than one and the \code{rwo} element is a data.frame with rules in \code{rwi}
#' obtained removing the transcription factor \code{TFi}. Furthermore, if the
#' input argument \code{figures} is set to TRUE, also the plots of the distributions
#' of support, confidence and lift of the rules before and after removing the
#' transcription factor TFi are provided.
#' @export
#' @importFrom stats var
#' @importFrom graphics par matplot legend title
#'
#' @examples
#' # Load r_FOSL2 and r_noFOSL2 from the data_man collection of datasets:
#' data('data_man')
#'
#' # The Importance Indexes of FOSL2=1 in the set of rules r_FOSL2 are given by:
#' IComp('FOSL2=1', r_FOSL2, r_noFOSL2, figures=TRUE)


IComp <- function(TFi, rules_TF, rules_noTF, figures) {
    if (length(TFi) > 1) {
        TFi <- TFi[[1]]
        for (i in 2:length(TFi)) TFi <- paste(TFi, TFi[[i]], sep = ",")
    }
    both <- cbind(rules_TF, rules_noTF)
    if (all(is.na(both)) == TRUE) {
        return(list(imp = NA, delta = NA))
    } else {
        Z <- data.frame(both)
        if (dim(Z)[1] > 1) {
            if (all(colSums(both[, c(8, 9, 10)]) == 0)) {
                Z[, c(8, 9, 10)] <- 0
            }
            index <- c(3, 4, 5, 8, 9, 10)
            vars <- vapply(both[, index], var, numeric(1))
            means <- vapply(both[, index], mean, numeric(1))
            for (i in index) {
                Z[, i] <- (both[, i] - means[which(index == i)])/sqrt(vars[which(index ==
                                                                                     i)])
            }
            Z[, index[which(vars == 0)]] <- both[, index[which(vars ==
                                                                   0)]]
        }
    }


    # matrix of variations of the standardized measures
    Z <- data.frame(Z)
    colnames(Z) <- colnames(both)
    diff_supp_Z <- Z[, 3] - Z[, 8]
    diff_conf_Z <- Z[, 4] - Z[, 9]
    diff_lift_Z <- Z[, 5] - Z[, 10]
    diffs_Z <- data.frame(Z[, 1], Z[, 6], diff_supp_Z, diff_conf_Z, diff_lift_Z)
    # evaluation of the Importance Index of TF in each rule
    imp_Z_rule_0 <- apply(diffs_Z[, 3:5], 1, sum)
    # consider only positives Importance Indexes
    imp_Z_rule <- imp_Z_rule_0[imp_Z_rule_0 > 0]
    m <- max(rules_TF$lift)
    n <- dim(Z)[1]
    if (figures == TRUE) {
        # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
        TF_i <- paste0(strsplit(TFi, "=")[[1]][1], "=0")
        par(mfrow = c(1, 3))
        matplot(both[, c(3, 8)], type = "l", lwd = 2, ylab = "support",
                xlab = "rules ID", ylim = c(0, 1))
        legend("topright", c(paste("with", TFi), paste("with", TF_i)),
               col = c("black", "red"), lty = c(1, 2), lwd = c(2, 2))
        title(main = "Rules support distribution")
        matplot(both[, c(4, 9)], type = "l", lwd = 2, ylab = "confidence",
                xlab = "rules ID", ylim = c(0, 1))
        title(main = "Rules confidence distribution")
        legend("topright", c(paste("with", TFi), paste("with", TF_i)),
               col = c("black", "red"), lty = c(1, 2), lwd = c(2, 2))
        matplot(both[, c(5, 10)], type = "l", lwd = 2, ylab = "lift", xlab = "rules ID",
                ylim = c(0, m + 3))
        title(main = "Rules lift distribution")
        legend("topright", c(paste("with", TFi), paste("with", TF_i)),
               col = c("black", "red"), lty = c(1, 2), lwd = c(2, 2))
    }
    # Consider only the rules with positive Importance Indexes
    d_Z <- diffs_Z[imp_Z_rule_0 > 0, 3:5]
    rwi <- rules_TF[imp_Z_rule_0 > 0, ]
    rwo <- rules_noTF[imp_Z_rule_0 > 0, ]
    return(list(imp = imp_Z_rule, delta = d_Z, rwi = rwi, rwo = rwo))
}

