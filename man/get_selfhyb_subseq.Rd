\name{get_selfhyb_subseq}
\alias{get_selfhyb_subseq}
\alias{show_selfhyb_counts}
\alias{show_selfhyb_lengths}

\title{Get Self-Hybridizing Subsequences}
\description{
  This function finds the longest self-hybridizing subsequences present
  in RNA or DNA sequences.
}
\usage{
get_selfhyb_subseq(seq, minlen, type = c("RNA", "DNA"))
show_selfhyb_counts(L)
show_selfhyb_lengths(L)
}

\arguments{
  \item{seq}{character vector of RNA or DNA sequences}
  \item{minlen}{an integer specifying the minimum length in bases of the
    self-hybridizing subsequences.  Subsequences with length less than
    \code{minlen} will be ignored.}
  \item{type}{one of \code{"RNA"} or \code{"DNA"} depending on the type of
    sequences provided in \code{seq}.  Note that you cannot mix RNA and
    DNA sequences.}
  \item{L}{The output of \code{get_selfhyp_subseq}.}
}

\details{
  \code{get_selfhyb_subseq} finds the longest self-hybridizing
  subsequences of the specified minimum length.  It does this using
  suffix trees and the \code{getLongestSubstring} function provided by
  the Rlibstree package.

  These are defined to be the longest string that is found in both
  the input sequence, \code{seq}, and in its reverse complement.
}

\value{
  A list with an element for each sequence in \code{seq}.  The list will
  be named using \code{names(seq)}.

  Each element is itself a list with an element for each longest
  self-hybridizing subsequence (there can be more than one).  Each such
  element is yet another list with components:

  \item{starts}{integer vector giving the character start positions for
    the self-hybridizing subsequence in the sequence.}
  
  \item{rcstarts}{integer vector giving the character start positions
    for the reverse complement of the self-hybridizing subsequence in
    the sequence.}
}

\author{Seth Falcon}

\examples{
seqs = c(a="UGAGGUAGUAGGUUGUAUAGUU", b="UGAGGUAGUAGGUUGUGUGGUU",
         c="UGAGGUAGUAGGUUGUAUGGUU")

ans = get_selfhyb_subseq(seqs, minlen=3, type="RNA")
length(ans)

ans[["a"]]

show_selfhyb_counts(ans)
show_selfhyb_lengths(ans)

}

\keyword{ manip }

