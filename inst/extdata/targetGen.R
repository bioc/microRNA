 ### USAGE: go to the data directory and source this script...
 ### ie.  source("../inst/extdata/targetGen.R")
 ### If you extend it to other species, don't forget that you will have to also add man pages for the other data structures you make.


### Generic preprocessing of the file
hmat = read.delim("../../targetData/v5.txt.homo_sapiens", header = TRUE, skip=3, sep = "\t")
mmat = read.delim("../../targetData/v5.txt.mus_musculus", header = TRUE, skip=3, sep = "\t")


### Generic filtering function
testSpec = function(x, str){
     if(length(grep(str,as.character(x)))>0){
         1 }
     else { 0 }
}

#Process human
hsInd = apply(hmat, 1, function(x) testSpec(x, str ="hsa-miR"))

hsMat = hmat[hsInd>0,]
hsMat = as.matrix(hsMat)

## hsMap = list()
## hsMap = hsMat[,12]
## names(hsMap) = hsMat[,2]
##library(Biobase)
## hsTargets = l2e(as.list(hsMap))
## hsTargets = hsMap

hsTargets = data.frame(name = hsMat[,"SEQ"], target = hsMat[,"TRANSCRIPT_ID"], chrom = hsMat[,"CHR"], start = hsMat[,"START"], end = hsMat[,"END"], strand = hsMat[,"STRAND"])

rownames(hsTargets) = NULL 

save(hsTargets, file = "hsTargets.rda")



#Process mouse
mmInd = apply(mmat, 1, function(x) testSpec(x, str ="mmu-miR"))

mmMat = mmat[mmInd>0,]
mmMat = as.matrix(mmMat)

## mmMap = list()
## mmMap = mmMat[,12]
## names(mmMap) = mmMat[,2]
## ## mmTargets = l2e(as.list(mmMap))
## mmTargets = mmMap

mmTargets = data.frame(name = mmMat[,"SEQ"], target = mmMat[,"TRANSCRIPT_ID"], chrom = mmMat[,"CHR"], start = mmMat[,"START"], end = mmMat[,"END"], strand = mmMat[,"STRAND"])

rownames(mmTargets) = NULL

save(mmTargets, file = "mmTargets.rda")



