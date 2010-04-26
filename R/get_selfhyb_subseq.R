.buildEl = function(s, r, minlen) {
    both = c(s, r)
    slen = nchar(s, type="chars")
    pkg = "Rlibstree"
    library(pkg, character.only=TRUE)
    autoseqs = Rlibstree::getLongestSubstring(Rlibstree::SuffixTree(both),
                                   repeated=FALSE,
                                   range=c(minlen, 0))
    if (is.null(autoseqs))
      return(list())
    names(autoseqs) = autoseqs
    selfhybs = lapply(autoseqs,
      function(x) {
          locL = gregexpr(x, both, fixed=TRUE)
          starts = as.integer(locL[[1]]) # remove match.length attr
          mlen = attr(locL[[2]], "match.length")
          rcstarts = slen - (as.integer(locL[[2]]) + mlen - 1) + 1
          list(starts=starts, rcstarts=rcstarts)
          ## FIXME: add gap calculation here
      })
    ## FIXME: add base count/freq
    selfhybs
}

get_selfhyb_subseq = function(seq, minlen, type=c("RNA", "DNA")) {
    pkg = "Rlibstree"
    suppressWarnings(require(pkg, character.only=TRUE, quietly=TRUE)) ||
        stop("you need to install the 'Rlibstree' package in order to\n",
             "  use get_selfhyb_subseq()")
    type = match.arg(type)
    revcompfun = switch(type,
      RNA=RNAStringSet,
      DNA=DNAStringSet,
      stop("type must be RNA or DNA"))
    rcSeq = reverseComplement(revcompfun(seq))
    rcSeq = as.character(rcSeq)
    ans = mapply(.buildEl, seq, rcSeq, minlen, USE.NAMES=FALSE)
    names(ans) = names(seq)
    ans
}

show_selfhyb_lengths = function(L) {
    lens = sapply(L, function(x) {
        if (length(x))
          nchar(names(x)[1], type="char")
        else
          0
    })
    table(lens, dnn="SHS length")
}

show_selfhyb_counts = function(L) {
    cnts = sapply(L, length)
    table(cnts, dnn="SHS counts")
}
