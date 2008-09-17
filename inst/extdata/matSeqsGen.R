 ### USAGE: go to the data directory and source this script...
 ### ie.  source("../inst/extdata/matSeqsGen.R")
 ### If you extend it to other species, don't forget that you will have to also add man pages for the other data structures you make.


 ### Generic preprocessing of the file
mat = read.delim("../inst/extdata/mature.fa", header = FALSE, sep = " ")

oInd = seq(1,dim(mat)[1],by=2)
eInd = seq(2,dim(mat)[1],by=2)

matO = mat[oInd,]
matE = mat[eInd,]

cMat = cbind(matO,matE)
cMat = cMat[,1:6]

fMat = apply(cMat, 1, function(x) gsub("^>","", x))
fMat = t(fMat)


 ### Generic filtering function
testSpec = function(x, str){
     if(length(grep(str,as.character(x)))>0){
         1 }
     else { 0 }
}



 ### Then Filter out the human rows
hsaDim = apply(fMat, 1, function(x) testSpec(x, str ="sapiens"))
hsaDim = as.vector(unlist(hsaDim))
hsaMat = fMat[hsaDim>0,]

 ### Make a list and save it:
hsSeqs = list()
hsSeqs = hsaMat[,6]
names(hsSeqs) = hsaMat[,1]

 ### Then just use the save() to save this thing as an .Rda file.
save(hsSeqs, file = "hsSeqs.rda")



 ### Then Filter out the mouse rows
mmuDim = apply(fMat, 1, function(x) testSpec(x, str ="musculus"))
mmuDim = as.vector(unlist(mmuDim))
mmuMat = fMat[mmuDim>0,]

 ### Make a list and save it:
mmSeqs = list()
mmSeqs = mmuMat[,6]
names(mmSeqs) = mmuMat[,1]

 ### Then just use the save() to save this thing as an .Rda file.
save(mmSeqs, file = "mmSeqs.rda")

