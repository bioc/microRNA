\name{seedRegions}
\alias{seedRegions}
\title{ A function to retrieve the seed regions from microRNA sequences}
\description{
The seed region of a microRNA consists of a set of nucleotides at the
5' end of the microRNA, typically bases 2 through 7, although some times 8
is used.
}
\usage{
seedRegions(x, start = 2, stop = 7)
}
\arguments{
  \item{x}{ A vector of microRNA sequences. }
  \item{start}{ The start locations, can be a vector. }
  \item{stop}{ The stop locations, can be a vector. }
}
\details{
We use \code{\link{substr}} to extract these sequences.
}
\value{
 A vector of the same length as \code{x} with the substrings.
}
\author{ R. Gentleman }
\seealso{ \code{\link{substr}}}
\examples{
 data(hsSeqs)
 seedRegions(hsSeqs[1:5])
 seedRegions(hsSeqs[1:3], start=c(2,1,2), stop=c(8,7,9))
}
\keyword{ manip }
