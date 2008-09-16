 RNA2DNA = function(x) chartr("U", "T", toupper(x))

 seedRegions = function(x, start=2, stop=7)
	substr(x, start, stop)


 matchSeeds = function(seeds, seqs) {
        m1 = lapply(seeds, gregexpr, text=seqs, ignore.case=TRUE)
        names(m1) = names(seeds)
        hasMatch = lapply(m1, function(x) sapply(x, function(y) y[y>0]))
        noMatches = sapply(hasMatch, function(x) sum(sapply(x, length)) > 0 )
        hasM = hasMatch[noMatches]
        whichM = lapply(hasM, function(x) {
                            hasMM = sapply(x, length) > 0
                            y = x[hasMM]
                            names(y) = names(seqs)[hasMM]
                            y})
        names(whichM) = names(seeds)[noMatches]
        return(whichM)
  }


