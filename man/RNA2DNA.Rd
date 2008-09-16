\name{RNA2DNA}
\alias{RNA2DNA}
\title{ A Function to translate RNA sequences into DNA sequences. }
\description{
RNA and DNA differ in that RNA uses uracil (U) and DNA uses thiamine (T),
this function translates an RNA sequence into a DNA sequence by translating
the characters.
}
\usage{
RNA2DNA(x)
}
\arguments{
  \item{x}{ A valid RNA sequence. }
}
\details{
No checking for validity of sequence is made, and the input sequence is
translated to upper case.
}
\value{
 A character vector, of the same length as \code{x} where all characters
are in upper case, and any instance of \code{U} in \code{x} is replaced
by a \code{T}.
}
\author{ R. Gentleman }
\seealso{ \code{\link{chartr}}}
\examples{
 input = c("AUCG", "uuac")
 RNA2DNA(input)
}
\keyword{ manip }
