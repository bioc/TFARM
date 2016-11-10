#' Generates an heatmap visualization of the mean importances of pairs of transcription factors.
#'
#' A square matrix (TFs x TFs) is built. The element (i,j) of such matrix contains the mean Importance Index of
#' the couple of transcription factors (TFi, TFj). This is rapresented in a heatmap visualization in which the scale color from red (dark) to white (light) indicates
#' low mean importance of the couple (red), or high mean importance of the couple (white).
#'
#' @param TFs a string vector with names of transcription factors
#' @param I a vector of mean importances of pairs of transcription factors in TFs
#'
#' @export
#'
#' @import fields
#' @examples
#' # Load the data:
#' data("data_man")
#'
#' # Heatmap visualization of the mean importances of transcription factors in p
#' # and their combinations in two elements:
#' heatI(p, I_c_2)

heatI <- function(TFs,I){
  l <- length(TFs)
  estract<- function(x)
  {
    vec <- rep(0, length(TFs))
    for (i in 1:length(vec))
    {
      vec[i] <- gregexpr(TFs[i], x)[[1]][1] != -1
    }
    return (vec)
  }
  aa = t(sapply(I[,1],estract))
  matrix_imp <- matrix(0, l, l)
  for ( i in 1:dim(aa)[1])
  {
    pos = rep(NA,2)
    pos[1] <- which(aa[i,] == 1)[1]
    pos[2] <- which(aa[i,] == 1)[2]
    if (is.na(pos[2]))
    {
      matrix_imp[pos[1], pos[1]] = I[i,]$imp
    }
    if (!is.na(pos[2]))
    {
      matrix_imp[pos[1], pos[2]]  = I[i,]$imp
      matrix_imp[pos[2], pos[1]] = I[i,]$imp
    }
  }
  rownames(matrix_imp) = TFs
  colnames(matrix_imp) = TFs
  par(mar = c(4,4,4,10))
  image(matrix_imp, axes=FALSE,xlab='', ylab='',col=heat.colors(20),ylim=c(1.04,-0.04))
  axis(1, at=seq(0,1, length=l), labels=rownames(matrix_imp), las=2, cex.axis=0.8)
  axis(2, at=seq(0,1, length=l), labels=rownames(matrix_imp), las=2, cex.axis=0.8)
  image.plot(matrix_imp,col=heat.colors(20),legend.only=TRUE,ylim=c(1.04,-0.04))
  title(main='Importance index of pairs of Transcription Factors')
}
