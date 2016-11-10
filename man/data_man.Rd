\name{data_man}
\alias{data_man}
\alias{DELTA}
\alias{I}
\alias{IMP_Z}
\alias{I_c_2}
\alias{imp_FOSL2}
\alias{p}
\alias{r_TEAD4}
\alias{r_noFOSL2}
\alias{r_FOSL2}
\docType{data}
\title{
    data_man
}
\description{
    It contains the datasets used in all the examples of the \code{TFARM} package.
}
\usage{data("data_man")}
\format{
    It contains four datasets:
        \describe{
            \item{\code{r_TEAD4}}{A data.frame with 30 rows (i.e., association rules) and 5 colums: the first column contains the left-hand-side of each rule, the second column contains the right-hand-side of each rule, the third, fourth and fifth column contain respectively support, confidence and lift of each rule. It is obtained with the \code{rulesGen} function, applied to the dataset \code{MCF7_chr1} contained in the package.}
            \item{\code{p}}{A string vector containg the 12 transcription factors present in the rules \code{r_TEAD4}.}
            \item{\code{r_FOSL2}}{A data.frame with 28 rows and 5 columns obtained with the function \code{rulesTF}, subsetting the \code{r_TEAD4} data.frame and extracting only the rules that contain \code{FOSL2=1} in the left-hand-side}
            \item{\code{r_noFOSL2}}{A data.frame with 28 rows and 5 columns obtained with the function \code{rulesNTF}, removing \code{FOSL2=1} from the left-hand-side of the rules in \code{r_FOSL2} and searching the obtained rules in the data.frame \code{r_TEAD4}.}
            \item{\code{imp_FOSL2}}{A vector containing the 28 Importance Indexes of \code{FOSL2=1} in the 28 rules contained in \code{r_FOSL2}. It is obtained with the function \code{IComp}.}
            \item{\code{DELTA}}{A list of 12 data.frame. Each element of the list is a data.frame with 3 columns correspondent to the variations of standardized support, confidence and lift obtained removing a transcription factor in \code{p} from the rules in \code{r_FOSL2}, respectively.}
            \item{\code{I}}{A data.frame with 12 rows and 3 columns containing the mean importance of the 12 transcription factors in \code{p} in the rules in \code{r_FOSL2}, and the number of rules in \code{r_FOSL2} in which each transcription factor in \code{p} appears.}
            \item{\code{IMP_Z}}{A list of 12 numerical vectors. Each element of the list contains the Importance Indexes of the transcription factors in \code{p} in the rules in \code{r_FOSL2}.}
            \item{\code{I_c_2}}{A data.frame with 67 rows and 2 columns that contains the mean importances of the transcription factors in \code{p} and their combination in couples, in the rules in \code{r_FOSL2}.}
        }
}

\examples{
    data(data_man)
}
\keyword{datasets}
