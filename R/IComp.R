
#' Computes the Importance Index of a transcription factor in a set of association rules
#'
#' Given an association rule and a transcription factor \code{TFi}, it is evaluated the contribution of \code{TFi} in the rule for the prediction of the presence of the item in the right-hand-side of the rule.
#' Since this contribution is evaluated based on the variations of support, confidence and lift of the rule, the user can visualize such variations setting the parameter \code{figures = TRUE}.
#'
#' @param TFi string or string vector: the transcription factor (or combination of transcription factors) whose importance distribution is evaluated
#' @param rules_TF a set of rules in which \code{TFi} is present
#' @param rules_noTF a set of rules obtained from rules_TF removing transcription factor (or combination of transcription factors) in TFi (obrained with the function \code{\link{rulesNTF}})
#' @param figures logical if \code{figures = TRUE}, graphics with support, lift and confidence distributions of the rules in \code{rules_TF} and \code{rules_noTF} are returned
#'
#' @return A list of two elements: the \code{imp} element is a vector of integers with the importances of TFi in each rule in \code{rules_TF}; the \code{delta} element of the list is a list with variations of standardized distributions of the three measures of support, confidence and lift. This output is used in the function \code{\link{IPCA}} for the Principal Component Analysis of such distributions.
#' @examples
#' # once defined rulesTAF1, rulesTAF1not from
#' # rulesTF and rulesNTF, IComp computes the importance of TAF1 = 1
#' # IComp("TAF1=1", rulesTAF1, rulesTAF1not, figures=TRUE)
#' @export
#'

IComp <- function(TFi, rules_TF, rules_noTF, figures){
  #TFi_0 <- lapply(items(TFi), function(x){return(paste(x,"1",sep="="))})
  #if(length(TFi_1) == 1) {TFi <- TFi}
  if(length(TFi) > 1) {
    TFi <- TFi[[1]]
    for (i in 2:length(TFi))
      TFi<- paste(TFi, TFi[[i]], sep=',')}
  both <- cbind(rules_TF,rules_noTF)
  if (all(is.na(both)) == TRUE) {
    return(list("imp"=NA,"delta"=NA)) }
  else {
      Z <- data.frame(both)
      if (dim(Z)[1] > 1) {
          for (i in c(3,4,5,8,9,10)){
              if (all(colSums(both[,c(8,9,10)])== 0)){
                  Z[,c(8,9,10)] <- 0
              }
              else if (var(both[,i]) == 0) Z[,i] <- both[,i]
              else  Z[,i] <- (both[,i]-mean(both[,i]))/(sqrt(var(both[,i])))
          }
      }
    Z <- data.frame(Z)
    colnames(Z) <- colnames(both)
    diff_supp_Z <- abs((Z[,3])-(Z[,8]))
    # check:
    diff_conf_Z <- abs((Z[,4])-(Z[,9]))
    diff_lift_Z <- abs((Z[,5])-(Z[,10]))
    diffs_Z <- data.frame(Z[,1], Z[,6], diff_supp_Z, diff_conf_Z, diff_lift_Z)
    imp_Z_rule <- apply(diffs_Z[,3:5], 1, sum)
    m <- max(Z[,5])
    n <- dim(Z)[1]
    if (figures == TRUE){
      #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
      par(mfrow=c(1,3))
      matplot(both[,c(3,8)], type='l', lwd=2, ylab='support', xlab='rules ID',ylim=c(0,1))
      legend('topright', c(paste('with',TFi),paste('without',TFi)), col=c('black','red'), lty=c(1,2),lwd=c(2,2))
      title(main='Rules support distribution')
      #plot.new()
      matplot(both[,c(4,9)], type='l', lwd=2, ylab='confidence', xlab='rules ID', ylim=c(0,1))
      title(main='Rules confidence distribution')
      legend('topright', c(paste('with',TFi),paste('without',TFi)), col=c('black','red'), lty=c(1,2),lwd=c(2,2))
      #plot.new()
      matplot(both[,c(5,10)], type='l', lwd=2, ylab='lift', xlab='rules ID', ylim=c(0,m+5))
      title(main='Rules lift distribution')
      legend('topright', c(paste('with',TFi),paste('without',TFi)), col=c('black','red'), lty=c(1,2),lwd=c(2,2))
    }
    return(list("imp"=imp_Z_rule,"delta"=diffs_Z[,3:5]))
  }
}

