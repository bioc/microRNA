
## sample script that computes the maximum ``sequence similarity''
## between each non-human miRNA in miRBase and human miRNAs.


library(Biostrings)

read.fa <- function(file)
{
    mirna <- readFASTA(file, strip.desc = TRUE)
    mirna.seqs <- sapply(mirna, "[[", "seq")
    names(mirna.seqs) <-
        sapply(strsplit(sapply(mirna, "[[", "desc"), " "), "[", 1)
    mirna.seqs
}

fa.file <- "miRBase-12.0-mature.fa.gz"

if (!file.exists(fa.file))
    download.file("ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/12.0/mature.fa.gz",
                  destfile = fa.file)


##mirna.seqs <- read.fa("miRBase-12.0-mature.fa")

mirna.seqs <- read.fa(gzfile(fa.file))

human.mirna <- mirna.seqs[substring(names(mirna.seqs), 1, 3) == "hsa"]
nonhuman.mirna <- mirna.seqs[substring(names(mirna.seqs), 1, 3) != "hsa"]


sim.with.human <-
    function(pattern, human.seqs = human.mirna,
             type = c("match", "align"), verbose = FALSE)
{
    type <- match.arg(type)
    sub.mat <-
        switch(type,
               match = {
                   ans <- matrix(-10000L, 4, 4,
                                 dimnames = list(c("A", "U", "C", "G"),
                                                 c("A", "U", "C", "G")))
                   diag(ans) <- 1L
                   ans
               },
               align = {
                   ans <- matrix(-1L, 4, 4,
                                 dimnames = list(c("A", "U", "C", "G"),
                                                 c("A", "U", "C", "G")))
                   diag(ans) <- 1L
                   ans
               })
    gap <- 
        switch(type,
               match = -10000L,
               align = -1L)
    sim <-
        pairwiseAlignment(human.seqs, pattern,
                          type = "local",
                          substitutionMatrix = sub.mat,
                          gapOpening = gap,
                          gapExtension = -10000L,
                          scoreOnly = FALSE)
    scores <- score(sim)
    w <- which.max(scores)
    if (verbose) {
        cat(sprintf("%s : %s -> %3g",
                    human.seqs[w], pattern, scores[w]),
            fill = TRUE)
        cat(paste(as.character(sim[w]), collapse = "\n"), fill = TRUE)
    }
    scores[w]
}


computeSimilarity <-
    function(query.seqs = nonhuman.mirna,
             target.seqs = human.mirna,
             ...)
{
    ans <- numeric(length(query.seqs))
    names(ans) <- names(query.seqs)
    ntot <- length(ans)
    for (i in seq_along(ans))
    {
        cat(sprintf("\r[%5g/%5g]", i, ntot))
        ans[i] <- 
            sim.with.human(query.seqs[i], target.seqs, ...)
    }
    ans
}


similarity.match <- 
    computeSimilarity(nonhuman.mirna,
                      human.mirna,
                      type = "match",
                      verbose = TRUE)

similarity.align <- 
    computeSimilarity(nonhuman.mirna,
                      human.mirna,
                      type = "align",
                      verbose = TRUE)

save(similarity.align, similarity.match,
     file = "similarity.rda")

