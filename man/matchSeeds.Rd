\name{matchSeeds}
\alias{matchSeeds}
\title{ A function to match seed regions to sequences. }
\description{
 Given an input set of seed regions and a set of sequences all locations
of the seed regions (exact matches) within the sequences are found.
}
\usage{
matchSeeds(seeds, seqs)
}
\arguments{
  \item{seeds}{ The seeds, or short sequences, to match. }
  \item{seqs}{ The sequences to find matches in. }
}
\details{
We presume that the problem is an exact matching problem and that
all sequences are in the correct orientation for that.  If, for
example, you start with seed regions from a microRNA (for \code{seeds})
 and 3'UTR sequences (for \code{seqs}), then you would want to reverse
complement one of the two sequences. And make sure all sequences
are either DNA or RNA.

 Names from either \code{seeds} or \code{seqs} are propogated, as much
as is possible.
}
\value{
A list containing one entry for each element of \code{seeds} that had
at least one match in one entry of \code{seqs}. Each element of this list
is a named vector containing the elements of \code{seqs} that the corresponding
seed has an exact match in.
}
\author{ R. Gentleman }
\seealso{\code{\link{seedRegions}}}
\examples{
data(hsSeqs)
data(s3utr)
hSeedReg = seedRegions(hsSeqs)
comphSeed = lapply(hSeedReg,RNAString)
comphSeed = lapply(comphSeed,reverseComplement)
comphSeed = unlist(lapply(comphSeed, function(x){as.character(x)}))
comph = RNA2DNA(comphSeed)
mx = matchSeeds(comph, s3utr)
}
\keyword{ manip }
