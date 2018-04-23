#' Principal Components Analysis for distributions of variations of support and 
#' confidence obtained removing a transcription factor
#' from a set of rules.
#'
#' The function is used for validation of the defined Importance Index.
#' It is defined as the linear combination of variations of
#' support and confidence measures.
#' The Principal Component Analysis lets the user evaluate
#' if in the 1D reference system defined by such linear combination,
#' it is possible to describe the variability of the data.
#'
#' @param delta_list list of variations of distributions
#' of support and confidence measures, obtained using the \code{IComp}.
#' @param IMP the importance matrix with the mean Importance Index
#' of every candidate co-regulator transcription factor and the number of rules
#' in which each of them appears.
#' @return Variance explained by every principal component (\code{summary}),
#' scores (i.e., the coordinates) of data in delta_list in the reference system
#' defined by the principal components (\code{scores}) and loadings
#' (i.e., the coefficinets) of the linear combination that defines
#' each principal component (\code{loadings}).
#' The plots of the variability, the cumulate percentage of variance explained
#' by each principal component and loadings of every principal component are
#' also returned.
#'
#' @export
#' @importFrom stats princomp
#' @importFrom graphics layout barplot plot title abline box
#'
#' @examples
#' # Load IMP, DELTA and TF_Imp from the data_man collection of datasets:
#' data('data_man')
#'
#' colnames(IMP)
#' TF_Imp <- data.frame(IMP$TF, IMP$imp, IMP$nrules)
#' i.pc <- IPCA(DELTA, TF_Imp)
#' names(i.pc)


IPCA <- function(delta_list, IMP) {
    D_2 <- unlist(delta_list)[which(!is.na(unlist(delta_list)))]
    Delta_0 <- data.frame(TF = rep(IMP[, 1], IMP[, 3]), matrix(D_2,
                                                                 ncol = 2, byrow = TRUE))
    colnames(Delta_0)[2:3] <- c("delta s", "delta c")
    Delta <- Delta_0[, 2:3]
    pc <- princomp(Delta, cor = FALSE, scores = TRUE)
    layout(matrix(c(1, 2, 3, 3, 4, 4, 5, 5), 4, 2, byrow = TRUE), heights = c(1,
                                                                              0.6, 0.6, 0.6))
    barplot(pc$sdev^2, las = 2, main = "Variances of the principal components",
            ylab = "Variances", ylim = c(0, 1), cex.main = 1)
    plot(cumsum(pc$sdev^2)/sum(pc$sde^2), type = "b", axes = TRUE,
         xlab = "number of components", ylab = "contribute to total variance")
    box()
    title(main = "Cumulate variances of the principal components", cex.main = 1)
    scores <- pc$scores
    load <- pc$loadings
    for (i in 1:2) {
        barplot(load[, i], ylim = c(-1, 1), ylab = paste("Comp.", i,
                                                           sep = ""))
        abline(h = 0)
        if (i == 1)
            title(main = "Loadings of the principal components")
    }
    return(list(summary = summary(pc), scores = scores, loadings = load))
}
