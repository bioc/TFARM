#' Visualizes with boxplots the importance index distribution of a set of transcription factors.
#'
#'  For a set of candidate co-regulator transcription factors, the Importance Index distribution of each transcription factor is plotted in boxplots.
#'  The shape of every boxplot depends on: the dimension of the distribution, which is equal to the number of rules in which every transcription factor appear (the higher is such number, the larger is a boxplot), the variability of
#'  the distribution (the higher is the variability of the importance index distribution, the longer is a boxplot). Moreover, The higher is the median of the Importance Index distribution for a candidate co-regulator transcription factor, the higher the boxplot is aligned with respect to the y axis.
#'
#' @param I a list of importance distributions
#' @param TFs string vector with the names of the transcription factors whose importance distributions are in I
#'
#' @examples
#' # To visualize the Importance Index distribution of TCF12=1, TAF1=1 and EP300=1:
#' # TFs <- c("TCF12=1", "TAF1=1", "EP300=1")
#' # IMP_TCF12 <- IComp("TCF12=1", rules_TCF12, rules_noTCF12, figures=FALSE)
#' # IMP_TAF1 <- IComp("TAF=1", rules_TAF1, rules_noTAF1, figures=FALSE)
#' # IMP_EP300 <- IComp("EP300=1", rules_EP300, rules_noEP300, figures=FALSE)
#' # IMP <- list(IMP_TCF12, IMP_TAF1, IMP_EP300)
#' # distribViz(IMP, TFs)
#' @export
#'
#'

distribViz = function(I,TFs){
  lungh <- sapply(I, function(x){length(x)})
  ll = unlist(I)
  df <- data.frame(importance=ll, TF=rep(TFs,lungh))
  boxplot(df$importance ~ df$TF, varwidth =TRUE, cex.lab=0.6, las=2)
  title(main='Transcription Factor Importance index distributions')
}
