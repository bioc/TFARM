#' Substitutes the presence of a transcription factor (or a combination of
#' transcription factors) in the left-hand-side of a set of rules,
#' with its absence.
#'
#' The function substitutes the presence of a given transcription factor
#' \code{TFi} (or a combination of transcription factors) chosen by the
#' user with its absence, in the subset of relevant association rules
#' extracted with the function \code{\link{rulesTF}}. Then it searches
#' for the obtained rules and their quality measures of
#' support, confidence and lift in the set of most relevant associations
#' extracted with the function \code{\link{rulesGen}}.
#' If a rule is not found in such set, it is searched in all the association
#' rules generable from the considered dataset using the function
#' \code{\link{search_rule}}. If a rule is not found in the set of all
#' the association rules, its quality measures are set to 0.
#'
#' @param TFi a string, or a string vector: transcription factor
#' (or combination of transcription factors) to remove from the set of rules.
#' @param sub_rules a data.frame with a subset of rules containing \code{TFi},
#' and their quality measures of support, confidence and lift (i.e., rules from
#' which the user wants to remove \code{TFi}).
#' @param all_rules a data.frame with a set of all the rules and their quality
#' measures of support, confidence and lift, to be considered for the search of
#' the obtained rules and their quality measures.
#' @param data a GRanges object which contains the Indicator of presence matrix
#' i.e., a matrix with 1 and 0 values representing presence or absence,
#' respectively (in case other values different from 0 are present, all of them
#' are considered as representing presence).
#' @param RHS the right-hand-side of the considered association rules.
#'
#' @return A data.frame with all the rules in the set \code{sub_rules}
#' in which the transcription factor (or combination of transcription factors)
#' \code{TFi} is absent, and their quality measures of support, confidence
#' and lift.
#' @export
#' @importFrom GenomicRanges elementMetadata as.data.frame
#'
#' @examples
#' # Load r_TEAD4 and r_FOSL2 from the data_man collection of datasets:
#' data("data_man")
#' # Load MCF7_chr1:
#' data("MCF7_chr1")
#'
#'
#' r_noFOSL2 <- rulesTF0("FOSL2=1", r_FOSL2, r_TEAD4, MCF7_chr1, "TEAD4=1")
#'


rulesTF0 <- function(TFi, sub_rules, all_rules, data, RHS){

    # Selection of the Indicator of presence matrix, where other values
    # different from 0 are considered as representing presence and are set to 1
    data<-as.data.frame(elementMetadata(data))
    1->data[data!=0]
    # Analysis on the Indicator of presence matrix
    if(length(TFi) == 1) {TF_i <- TFi}
    if(length(TFi) > 1) {
        TF_i <- TFi[[1]]
        for (i in 2:length(TFi))
            TF_i<- paste(TF_i, TFi[[i]], sep=',')
        }
    K <- length(items(TF_i))

    TF_vp <- paste("", items(TF_i), sep=",")
    TF_vd <- paste(items(TF_i), "", sep=",")


    # substitute "TF=1" with "TF=0"
    # case with 1 TF
    if (length(items(TFi)) == 1) {
        rule_noTF <- lapply(sub_rules$lhs, function(x){
            r <- items(x)
            TF <- unlist(strsplit(TF_i,"="))[1]
            r[r%in%TF_i] <-  paste(TF,0,sep="=")
            return(r)
        })
    }


	 # case with 2 or more TF
	 else {
	     rule_noTF <- lapply(sub_rules$lhs, function(x){
	         r <- items(x)
	         TF <- lapply(items(TFi), function(x){
	             return(unlist(strsplit(x, "="))[1])
	         })
	         TFs_new <- paste(TF,0,sep="=")
	         r[r%in%items(TFi)] <-  paste(TFs_new)
	         return(r)
	     })
	 }

	  # inverse function of items
	  subs_noTF <- lapply(rule_noTF, itemset)

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
            if (all(all_noTF[i,c(3,4,5)] == 0)) {
                out <- search_rule(data, subs_noTF[i], RHS)
                if ((length(out) == 1 && is.na(out)) || (length(out) > 1 &&
                                                         all(out == 'NA'))){
                    all_noTF[i,1] <- paste(subs_noTF[i])
                    all_noTF[i,2] <- paste("{", "}", sep=RHS)
                    all_noTF[i,3:5] <- c(0,0,0)
                }
                else {
                    all_noTF[i,1] <- paste(out$lhs)
                    all_noTF[i,2] <- paste(out$rhs)
                    all_noTF[i,3:5] <- search_rule(data, subs_noTF[i], RHS)[3:5]
                }
            }
        }
        return(all_noTF)
}

