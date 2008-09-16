\name{s3utr}
\alias{s3utr}
\docType{data}
\title{ Test sequence data}
\description{
A vector of 3' UTR sequence data, the names correspond to Entrez Gene
IDs and the data were extracted using biomaRt.
}
\usage{data(s3utr)}
\format{
 A character vector, the values are the 3' UTR for a set of genes, the
names are Entrez Gene Identifiers.
}
\details{
The data were downloaded using the \code{getSequence} function in the
biomaRt package and duplicate strings removed. There remain some duplicated 
Entrez IDs but the reported 3' UTRs are different.
}
\examples{
data(s3utr)
}
\keyword{datasets}
