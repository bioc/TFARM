#' Exctracts relevant association rules.
#'
#' From the dataset in \code{data}, the function extracts a set of association rules, with a certain item in the right-hand-side. Each rule extracted has support greater than \code{minsupp} and confidence greater than \code{minconf}. The extraction is made using the \link[arules]{apriori} function implemented in the \code{arules} package.
#' \code{minsupp} and \code{minconf} thresholds are set by the user in order to extract a limited number of most relevant association rules.
#' @param data a binary matrix or data frame
#' @param TF a string with the name of the trancription factor wanted in the right-hand-side of the extracted rules
#' @param minsupp an integer, the minimal support of the extracted rules
#' @param minconf an integer, the minimal confidence of the extracted rules
#' @param type a logical parameter; if \code{type = TRUE} only rules with present transcription factor in the left-hand-side are extracted (i.e., the left-hand-side of the extracted rules is of the type {TF1=1, TF2=1, TF3=1}). If \code{type = FALSE} also rules with absent transcription factors in the left-hand-side are extracted (i.e., the left-hand-side of the extracted rules can be of the type {TF1=1, TF2=0, TF3=1} or {TF1=0, TF2=0, TF3=0}).
#'
#' @return A data frame with the association rules extracted and their quality measures of support, confidence and lift.
#' @export
#' @import arules
#'
#' @examples
#' # To extract association rules from data, with TEAD4=1 in the left-hand-side and
#' # support greater than 0.2 and confidence greater than 0.8:
#' # r.TEAD4 <- rulesGen(data, "TEAD4=1", minsupp=0.2, minconf=0.8, type = TRUE)
#'
#' @seealso \link[arules]{apriori}
rulesGen <- function(data, TF, minsupp, minconf, type){

  m <- dim(data)[2]
  data <- data[,2:m]
  m <- dim(data)[2]
  data.f <- data
  # Transactions extraction
  for (i in 1:m) {
    data.f[,i] <- as.factor(data.f[,i])
  }
  trans <- as(data.f, "transactions")
  names <- colnames(data.f)
  #TF_0 <- paste(TF,"0",sep="=")
  #TF_1 <- paste(TF,"1",sep="=")
  TF_a <- unlist(strsplit(TF,'='))[1]
  TF_1 <- paste(TF_a, "1", sep="=")
  names_1 <- unlist(sapply(names, function(x){paste(x,"1", sep="=")}))
  names_1 <- names_1[which(names_1 != TF_1)]
  # Rules extraction
  if (type == 'TRUE'){
     rules <- apriori(data.f, parameter=list(supp=minsupp, conf=minconf, target="rules"),
                       appearance=list(lhs=names_1, rhs=TF, default="none"))
  }
  else {
     rules <- apriori(data.f, parameter=list(supp=minsupp, conf=minconf, target="rules"),
                       appearance=list(rhs=TF, default="lhs"))

  }
  if (length(rules) ==  0) print("No rules found. Change parameter thresholds!")
  else if (length(rules) > 0) {
  r.TF <- subset(rules, subset=rhs %ain% TF)
  rules.TF = data.frame(lhs = labels(lhs(r.TF)),
                        rhs = labels(rhs(r.TF)), r.TF@quality)
  return(rules.TF) }
}

