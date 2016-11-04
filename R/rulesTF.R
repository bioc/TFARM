#' Extracts a subset of rules that contain a certain transcription factor (or a combination of transcription factors) in their left-hand-side.
#'
#' From a set of relevant association rules, only the ones containing \code{TFi} in their left hand side are subsetted, together with their quality measures of support, confidence and lift. The function is then used for the evaluation of the importance distribution of \code{TFi}.
#'
#' @param TFi a string with the name of the transcription factor (or combination of transcription factors) wanted in the left-hand-side of the rules to find
#' @param rules a data.frame with association rules and their quality measures of support, confidence and lift
#' @param verbose a logical parameter. If \code{verbose = TRUE}, a console message warns the user when the set of rules containing \code{TFi} is empty.
#'
#' @return A data.frame with association rules containing \code{TFi} in their left-hand-side, and with their quality measures of support, confidence and lift.
#' @export
#'
rulesTF <- function(TFi, rules, verbose){
  #TFi <- lapply(items(TF_i), function(x){return(paste(x,"1",sep="="))})
  if(length(TFi) == 1) {TF_i <- TFi}
  if(length(TFi) > 1) {
    TF_i <- TFi[[1]]
    for (i in 2:length(TFi))
      TF_i<- paste(TF_i, TFi[[i]], sep=',')}
  n_rules <- dim(rules)[1]
  rules_subs <- NULL
  for (i in 1:n_rules){
    if ( all(items(TF_i) %in% items(rules$lhs[i])) )
      rules_subs <- cbind(rules_subs, paste(rules$lhs[i]))
  }
  ## apply() + function

  if (length(rules_subs) == 0) {
      if (verbose == 'TRUE') {
    print(paste('None of the rules contains', TF_i, sep=' '))
    return(NA)}
      else if (verbose =='FALSE') (return(NA))
  }

  else {  #loof for quality measures of the rule
    n_subs <- length(rules_subs)
    n_all <- dim(rules)[1]
    all_TF <- matrix(0, n_subs, 5)
    all_TF <- data.frame(all_TF)
    colnames(all_TF) <- c('lhs','rhs', 'support', 'confidence', 'lift')
    for (i in 1:n_subs){
      d <- items(rules_subs[i])
      for (j in 1:n_all){
        c <- items(rules$lhs[j])
        if (all(c %in% d)){
          all_TF[i,3:5] <- rules[j,c(3,4,5)]
          all_TF[i,2] <- paste(rules[j,2])
          all_TF[i,1] <- paste(rules[j,1])
        }
      }
    }

  }
  return(all_TF)
}

