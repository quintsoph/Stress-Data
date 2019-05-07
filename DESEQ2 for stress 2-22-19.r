#nohup R CMD BATCH ./myprog.R &
library("RColorBrewer")
 library("DESeq2")
 library("geneplotter")
 library("EDASeq")
 library("genefilter")
 library ("locfit")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("EDASeq")

#Load gene file
	setwd("/DataDrives/dd12/stress/L1HSoutput")
	batchnormal<- read.table("stressgenematrix.txt", sep="\t", fill=TRUE, header=TRUE, row.names=1)
	batchnormal2<- as.matrix(batchnormal)
	
#setup covairate tissue type							  
	setwd("/DataDrives/dd12/stress/L1HSoutput")
	attrib<-read.table("attributematrix.txt", fill=TRUE, sep="\t", row.names=1, header=TRUE)
	attribs<- as.matrix(attrib)
	
#create the R class

dds <- DESeqDataSetFromMatrix(countData = batchnormal2,
							colData = as.data.frame(attribs),
                              design = ~ treatment + rep + time + treatment:time)

# step 1
  FDR <- 0.1
  floorPDEG <- 0.05
  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
  dds <- dds[filter,]
  #size factors
ddsnormal <- estimateSizeFactors(dds)
    sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
    sf <- sizeFactors(ddsnormal) 
	setwd("/DataDrives/dd12/stress/L1HSoutput")
    write.table(sf, file="sizefactors.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
	sf <- read.table("sizefactors.1.txt", row.names=1, sep="\t", header=TRUE)
	  sf <- c(t(sf))
	  sizeFactors(ddsnormal) <- sf
	  colnormal = colData(ddsnormal)
	
#skip step 2 and 3 since it is the same tissue type
###Input TEtranscripts
	setwd("/DataDrives/dd12/stress/L1HSoutput")
totaltissue<- read.table("stressL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
ttissue<- t(totaltissue)
library("stringr")
t_TEcntmatrix = ttissue
 t_cntmatrix = t(batchnormal2)
m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
new_patients = m[,1]
new_genenames = colnames(m[,-1])
new_cntmatrix = as.matrix(t(m[,-1]))
colnames(new_cntmatrix) = new_patients
  write.table(new_cntmatrix, file="cntmatrix.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t") 
  
  
  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
  colnames(new_cntmatrix) = NULL
  
  
  ddsnew <- DESeqDataSetFromMatrix(
    countData = new_cntmatrix,
    colData = colnormal,
    design = ~ treatment + rep + time + treatment:time)
  ddsnew$type <- droplevels(ddsnew$type)
  # set sizefactors 
  sizeFactors(ddsnew) <- sf
  save(ddsnew, file="ddsNew.DEGES.RData")
  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
  write.table(ddsnewcnts, file="normalizedcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
 
 
  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
  write.table(L1HSnormalized, file="L1HS.normalized.txt", quote=FALSE, row.names=TRUE, sep="\t")
  log2colsums <- log2(colSums(ddsnewcnts))
  write.table(log2colsums, file="log2colsums.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
#log transformation (vst)
 if (file.exists("VSTcnts.DEGES.RData")) {
  load(file = "VSTcnts.DEGES.RData")
 } else {
  vsd <- vst(ddsnew, blind = FALSE)
  VSTcnts <- assay(vsd)
  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.txt", quote=FALSE, row.names = TRUE)
  save(VSTcnts, file = "VSTcnts.DEGES.RData")
  write.table(VSTcnts, file="VSTcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
  #genenames = rownames(ddsnew)
  #genenames <- gsub("ALR/Alpha", "ALR_Alpha", genenames)
  #genenames <- gsub("BSR/Beta", "BSR_Beta", genenames)
  #rownames(ddsnew) <- genenames
  #coldataNew = colData(ddsnew)
  #for (i in 1:nrow(VSTcnts)) {
   # fname = paste0("VSTcnts.normal/",genenames[i], ".txt")
    #cnts <- matrix(VSTcnts[i,], ncol=1)
    #rownames(cnts) <- coldataNew$patient
    #colnames(cnts) <- c("patient\tgene")
    #write.table(cnts, file=fname, quote=FALSE, sep="\t", row.names=TRUE)
  }
  L1HSnormalized = read.table("L1HS.normalized.txt", row.names=1, sep="\t", check.names = FALSE)
  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
  colnames(L1HStable)[5:6] = c("normalizedcnt", "VSTcnt")
  write.table(L1HStable, file="L1HS.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
 #subset ddsnew
 #use deseq function then results function 
 