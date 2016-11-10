#' Visualizes with boxplots the importance index distribution of a set of transcription factors.
#'
#'  For a set of candidate co-regulator transcription factors, the Importance Index distribution of each transcription factor is plotted in boxplots.
#'  The shape of every boxplot depends on: the dimension of the distribution, which is equal to the number of rules in which every transcription factor appear (the higher is such number, the larger is a boxplot), the variability of
#'  the distribution (the higher is the variability of the importance index distribution, the longer is a boxplot). Moreover, The higher is the median of the Importance Index distribution for a candidate co-regulator transcription factor, the higher the boxplot is aligned with respect to the y axis.
#'
#' @param I a list of importance distributions
#' @param TFs string vector with the names of the transcription factors whose importance distributions are in I
#'
#' @export
#'
#' @examples
#' # Load the data:
#' data("data_man")
#'
#' # Plot the Importance Index distributions of transcription factors in p:
#' distribViz(IMP_Z,p)

distribViz = function(I,TFs){
  lungh <- sapply(I, function(x){length(x)})
  ll = unlist(I)
  df <- data.frame(importance=ll, TF=rep(TFs,lungh))
  boxplot(df$importance ~ df$TF, varwidth =TRUE, cex.lab=0.6, las=2)
  title(main='Transcription Factor Importance index distributions')
}
