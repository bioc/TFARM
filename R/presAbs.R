#' Splits a set of transcription factors in 'present' or 'absent' transcription factors in a set of rules.
#'
#' The function is used to find candidate co-regulator transcription factors with a given transcription factor, by finding the transcription factors that contribute to the prediction of its presence in a set of association rules.
#'
#'
#' @param TFs a string vector: set of transcription factors to search in the set of rules.
#' @param rules data.frame: rules and their quality measures (support, confidence and lift).
#' @param type logical parameter: if \code{type = TRUE} only rules with present transcription factors in the left-hand-side are considered.
#'
#' @return A list of two string vectors: the list \code{pres} contains all the transcription factors in \code{TFs} that are present in \code{rules}, and the list \code{abs} contains all the transcription factors in \code{TFs} that are absent in \code{rules}.
#'
#' @export
#' @examples
#' # Load the data:
#' data("data_man")
#' data("MCF7_chr1")
#'
#' # Transcription factors present in at least one of the regions in the considered dataset:
#' m <- dim(MCF7_chr1)[2]
#' c <- colnames(MCF7_chr1)[2:m]
#'
#' names(presAbs(c, r_TEAD4, TRUE))
#'
#' # Transcription factors present in at least one of the association rules:
#' p <- presAbs(c, r_TEAD4, TRUE)$pres
#' p


presAbs = function(TFs, rules, type){
  TFs_1 <- unlist(sapply(TFs, function(x){paste(x,"1",sep="=")}))
  TFs_0 <- unlist(sapply(TFs, function(x){paste(x,"0",sep="=")}))
  TFs_v <- cbind(t(TFs_1),t(TFs_0))
  # To find TFs present in at least a rule
  ALL_pres <- sapply(TFs_v, function(y){
    pres <- sapply(rules[,1], function(x){
      if(items(y) %in% items(x)) return(1)
    })
    tot_pres <- sum(unlist(pres))
    return(tot_pres)
  })
  if (type == 'TRUE'){
      pp <- as.vector(TFs_1[which(ALL_pres>0)])
      aa <- as.vector(TFs_1[which(!TFs_1%in%pp)])
  }
  else {
  # Present TFs:
  pp <- TFs_v[which(ALL_pres>0)]
  # Absent TFs:
  aa <- TFs_v[which(!TFs_v%in%pp)]}
  return(list("pres"=pp,"abs"=aa))
}



