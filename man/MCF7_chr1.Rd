\name{MCF7_chr1}
\alias{MCF7_chr1}
\docType{data}
\title{
MCF7_chr1
}
\description{
Example of dataset for the use of the method proposed in TFARM package.
}
\usage{data("MCF7_chr1")}
\format{
  A data frame with 2944 rows and 26 columns. Each row represents a genomic region and each column is the variable presence/absence of a transcription factor in a genomic region. The element (i,j) of the matrix (with j>1) is equal to 1 if the transcription factor j is present in the region i, otherwise is equal to 0.
}
\details{
The dataset is obtained from computations of ENCODE ChIP-seq data of the localization of transcription factors binding sites and histone modifications in DNA, as well as RefSeq data of the localization of promotorial regions in DNA (https://www.ncbi.nlm.nih.gov/refseq/). Such data have been processed and extracted with GMQL (GenoMetric Query Language http://www.bioinformatics.deib.polimi.it/GMQL/) queries.}
\examples{
data(MCF7_chr1)

# To visualize the indicator of presence matrix as a heatmap:
m <- dim(MCF7_chr1)[2]
image(as.matrix(MCF7_chr1[,2:m]))
}
\keyword{datasets}
