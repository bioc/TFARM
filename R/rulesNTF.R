#' Removes a transcription factor (or a combination of transcription factors) from the left-hand-side of a set of rules.
#'
#' The function removes a given transcription factor \code{TFi} (or a combination of transcription factors) chosen by the user, from the subset of relevant association rules extracted with the function \code{\link{rulesTF}}. Then it searches for the obtained rules and their quality measures of support, confidence and lift
#' in the set of most relevant associations extracted with the function \code{\link{rulesGen}}. If a rule is not found in such set, its quality measures are set to 0.
#'
#' @param TFi a string, or a string vector: transcription factor (or combination of transcription factors) to remove from the set of rules.
#' @param sub_rules a data.frame with a subset of rules containing \code{TFi}, and their quality measures of support, confidence and lift (i.e., rules from which the user wants to remove \code{TFi}).
#' @param all_rules a data.frame with a set of rules and their quality measures of support, confidence and lift.
#'
#' @return a data.frame with all the rules in the set \code{rules} without the transcription factor (or combination of transcription factors) in \code{TFi}, and their quality measures of support, confidence and lift.
#' @export
#' @examples
#' # Load the data:
#' data("data_man")
#'
#' r_noFOSL2 <- rulesNTF("FOSL2=1", r_FOSL2, r_TEAD4)
#'

rulesNTF <- function(TFi, sub_rules, all_rules){
  if(length(TFi) == 1) {TF_i <- TFi}
  if(length(TFi) > 1) {
    TF_i <- TFi[[1]]
    for (i in 2:length(TFi))
      TF_i<- paste(TF_i, TFi[[i]], sep=',')}
  if (all(is.na(sub_rules)) == TRUE) {
    return(NA)}
  else {
    K <- length(items(TF_i))
    # remove TF_i
    TF_vp <- unlist(lapply(items(TF_i),function(x){return(paste("", x, sep=","))}))
    TF_vd <- unlist(lapply(items(TF_i),function(x){return(paste(x, "", sep=","))}))
    sub <- list()
    # remove TF_i from lhs
    subs_noTF <- sapply(sub_rules$lhs, function(x){
      sub[[1]] <- gsub(items(TF_i)[1],"",gsub(TF_vp[1],"",gsub(TF_vd[1],"",x)))
      if (K == 1) return(sub[[1]])
      else if (K > 1){
        for (i in 2:K){
          sub[[i]] <- gsub(items(TF_i)[i],"",gsub(TF_vp[i],"",gsub(TF_vd[i],"",sub[[i-1]])))
        }
        return(sub[[K]])
      }
    })

    n_subs_2 <- length(subs_noTF)
    n_all <- dim(all_rules)[1]
    all_noTF <- matrix(0, n_subs_2, 5)
    all_noTF <- data.frame(all_noTF)
    colnames(all_noTF) <- c('lhs','rhs', 'support', 'confidence', 'lift')
    all_noTF$lhs <- paste(as.vector(subs_noTF))
    all_noTF$rhs <- sub_rules$rhs
    for (i in 1:n_subs_2){
      b <- items(subs_noTF[i])
      for (j in 1:n_all){
        a <- items(all_rules$lhs[j])
        if (all(a %in% b)){
          all_noTF[i,3:5] <- all_rules[j,c(3,4,5)]
          all_noTF[i,2] <- paste(all_rules[j,2])
        }
      }
    }
    return(all_noTF)
  }
}
