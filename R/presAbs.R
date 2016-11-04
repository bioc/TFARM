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
#' @export
#'

presAbs = function(TFs, rules, type){
 # source('items.R')
  TFs_1 <- unlist(sapply(TFs, function(x){paste(x,"1",sep="=")}))
  TFs_0 <- unlist(sapply(TFs, function(x){paste(x,"0",sep="=")}))
  TFs_v <- cbind(t(TFs_1),t(TFs_0))
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
  #pp_2 <- unlist(strsplit(pp_1,"=1"))
  #pp <- names(pp_2[which(pp_2 != 'TEAD4=1')])
  # Absent TFs:
  aa <- TFs_v[which(!TFs_v%in%pp)]}
  #aa <- aa_2[which(aa_2 != 'TEAD4')]
  return(list("pres"=pp,"abs"=aa))
}



