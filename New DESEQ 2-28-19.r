#nohup R CMD BATCH ./myprog.R &
library("RColorBrewer")
 library("DESeq2")
 library("geneplotter")
 library("EDASeq")
 library("genefilter")
 library ("locfit")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("EDASeq")

#Load gene file
	setwd("/DataDrives/dd12/stress/L1HSoutput")
	batchnormal<- read.table("stressgenematrix.txt", sep="\t", fill=TRUE, header=TRUE, row.names=1)
	batchnormal2<- as.matrix(batchnormal)
	
#setup covairate tissue type							  
	setwd("/DataDrives/dd12/stress/L1HSoutput")
	attrib<-read.table("attributematrix.txt", fill=TRUE, sep="\t", row.names=1, header=TRUE)
	attribs<- as.matrix(attrib)
#treatment separate
#h202
h202sub<- subset(attribs, attribs[,1] == "H2O2", select=c(1,2,3))
h2o2cut<- batchnormal2[,match(rownames(h202sub), colnames(batchnormal2))]
	#separate by time
	#0
	time0h2o2<- subset(h202sub, h202sub[,3] == "0hr", select=c(1,2,3))
	time0h2o2cut<- h2o2cut[,match(rownames(time0h2o2), colnames(h2o2cut))]
	#1
	time1h2o2<- subset(h202sub, h202sub[,3] == "1hr", select=c(1,2,3))
	time1h2o2cut<- h2o2cut[,match(rownames(time1h2o2), colnames(h2o2cut))]
	#4
	time4h2o2<- subset(h202sub, h202sub[,3] == "4hr", select=c(1,2,3))
	time4h2o2cut<- h2o2cut[,match(rownames(time4h2o2), colnames(h2o2cut))]
	#8
	time8h2o2<- subset(h202sub, h202sub[,3] == "8hr", select=c(1,2,3))
	time8h2o2cut<- h2o2cut[,match(rownames(time8h2o2), colnames(h2o2cut))]
#TM
tmsub<- subset(attribs, attribs[,1] == "Tm", select=c(1,2,3))
tmcut<- batchnormal2[,match(rownames(tmsub), colnames(batchnormal2))]
	#separate by time
	#0
	time0tm<- subset(tmsub, tmsub[,3] == "0hr", select=c(1,2,3))
	time0tmcut<- tmcut[,match(rownames(time0tm), colnames(tmcut))]
	#1
	time1tm<- subset(tmsub, tmsub[,3] == "1hr", select=c(1,2,3))
	time1tmcut<- tmcut[,match(rownames(time1tm), colnames(tmcut))]
	#4
	time4tm<- subset(tmsub, tmsub[,3] == "4hr", select=c(1,2,3))
	time4tmcut<- tmcut[,match(rownames(time4tm), colnames(tmcut))]
	#8
	time8tm<- subset(tmsub, tmsub[,3] == "8hr", select=c(1,2,3))
	time8tmcut<- tmcut[,match(rownames(time8tm), colnames(tmcut))]
setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
##H2O2 0 and 1 test
	newattribfile<- rbind(time0h2o2, time1h2o2)
	newcountdata<- cbind(time0h2o2cut, time1h2o2cut)

	dds <- DESeqDataSetFromMatrix(countData = newcountdata,
								colData = as.data.frame(newattribfile),
								  design = ~ time)
								  
	# step 1
	  FDR <- 0.1
	  floorPDEG <- 0.05
	  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
	  dds <- dds[filter,]
	  #size factors
	ddsnormal <- estimateSizeFactors(dds)
		sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
		sf <- sizeFactors(ddsnormal) 
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
		write.table(sf, file="sizefactorsh2o201.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
		sf <- read.table("sizefactorsh2o201.1.txt", row.names=1, sep="\t", header=TRUE)
		  sf <- c(t(sf))
		  sizeFactors(ddsnormal) <- sf
		  colnormal = colData(ddsnormal)
		
	#skip step 2 and 3 since it is the same tissue type
	###Input TEtranscripts
		setwd("/DataDrives/dd12/stress/L1HSoutput")
	totaltissue<- read.table("stressL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	ttissue<- t(totaltissue)
	library("stringr")
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
	t_TEcntmatrix = ttissue
	 t_cntmatrix = t(newcountdata)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	colnames(new_cntmatrix) = new_patients
	write.table(new_cntmatrix, file="cntmatrix.h2o201.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t") 
	  
	  
	  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
	  colnames(new_cntmatrix) = NULL
	  ddsnew <- DESeqDataSetFromMatrix(
		countData = new_cntmatrix,
		colData = colnormal,
		design = ~ time)
	  ddsnew$type <- droplevels(ddsnew$type)
	  # set sizefactors 
	  sizeFactors(ddsnew) <- sf
	  save(ddsnew, file="ddsNew.DEGES.RData")
	  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
	  write.table(ddsnewcnts, file="normalizedcnt.h2o201.txt", quote=FALSE, row.names=TRUE, sep="\t")
	 
	 
	  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
	  write.table(L1HSnormalized, file="L1HS.normalized.h2o201.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  log2colsums <- log2(colSums(ddsnewcnts))
	  write.table(log2colsums, file="log2colsums.h2o201.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
	#log transformation (vst)
	 vsd <- vst(ddsnew, blind = FALSE)
	  VSTcnts <- assay(vsd)
	  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.h2o201.txt", quote=FALSE, row.names = TRUE)
	  save(VSTcnts, file = "VSTcnts.h2o201.DEGES.RData")
	  write.table(VSTcnts, file="VSTcnt.h2o201.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
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

	  L1HSnormalized = read.table("L1HS.normalized.h2o201.txt", row.names=1, sep="\t", check.names = FALSE)
	  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
	  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
	  colnames(L1HStable)[5:6] = c("normalizedcnt", "VSTcnt")
	  write.table(L1HStable, file="L1HS.h2o201.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
######
#####
###
##H2O2 0 and 4 test
	newattribfile<- rbind(time0h2o2, time4h2o2)
	newcountdata<- cbind(time0h2o2cut, time4h2o2cut)

	dds <- DESeqDataSetFromMatrix(countData = newcountdata,
								colData = as.data.frame(newattribfile),
								  design = ~ time)
								  
	# step 1
	  FDR <- 0.1
	  floorPDEG <- 0.05
	  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
	  dds <- dds[filter,]
	  #size factors
	ddsnormal <- estimateSizeFactors(dds)
		sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
		sf <- sizeFactors(ddsnormal) 
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
		write.table(sf, file="sizefactorsh2o204.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
		sf <- read.table("sizefactorsh2o204.1.txt", row.names=1, sep="\t", header=TRUE)
		  sf <- c(t(sf))
		  sizeFactors(ddsnormal) <- sf
		  colnormal = colData(ddsnormal)
		
	#skip step 2 and 3 since it is the same tissue type
	###Input TEtranscripts
		setwd("/DataDrives/dd12/stress/L1HSoutput")
	totaltissue<- read.table("stressL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	ttissue<- t(totaltissue)
	library("stringr")
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
	t_TEcntmatrix = ttissue
	 t_cntmatrix = t(newcountdata)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	colnames(new_cntmatrix) = new_patients
	write.table(new_cntmatrix, file="cntmatrix.h2o204.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t") 
	  
	  
	  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
	  colnames(new_cntmatrix) = NULL
	  ddsnew <- DESeqDataSetFromMatrix(
		countData = new_cntmatrix,
		colData = colnormal,
		design = ~ time)
	  ddsnew$type <- droplevels(ddsnew$type)
	  # set sizefactors 
	  sizeFactors(ddsnew) <- sf
	  save(ddsnew, file="ddsNew.DEGES.RData")
	  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
	  write.table(ddsnewcnts, file="normalizedcnt.h2o204.txt", quote=FALSE, row.names=TRUE, sep="\t")
	 
	 
	  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
	  write.table(L1HSnormalized, file="L1HS.normalized.h2o204.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  log2colsums <- log2(colSums(ddsnewcnts))
	  write.table(log2colsums, file="log2colsums.h2o204.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
	#log transformation (vst)
	 vsd <- vst(ddsnew, blind = FALSE)
	  VSTcnts <- assay(vsd)
	  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.h2o204.txt", quote=FALSE, row.names = TRUE)
	  save(VSTcnts, file = "VSTcnts.h2o204.DEGES.RData")
	  write.table(VSTcnts, file="VSTcnt.h2o204.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
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

	  L1HSnormalized = read.table("L1HS.normalized.h2o204.txt", row.names=1, sep="\t", check.names = FALSE)
	  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
	  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
	  colnames(L1HStable)[5:6] = c("normalizedcnt", "VSTcnt")
	  write.table(L1HStable, file="L1HS.h2o204.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
######
#####
###
##H2O2 0 and 8 test
	newattribfile<- rbind(timeh2o2, time8h2o2)
	newcountdata<- cbind(time0h2o2cut, time8h2o2cut)

	dds <- DESeqDataSetFromMatrix(countData = newcountdata,
								colData = as.data.frame(newattribfile),
								  design = ~ time)
								  
	# step 1
	  FDR <- 0.1
	  floorPDEG <- 0.05
	  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
	  dds <- dds[filter,]
	  #size factors
	ddsnormal <- estimateSizeFactors(dds)
		sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
		sf <- sizeFactors(ddsnormal) 
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
		write.table(sf, file="sizefactorsh2o208.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
		sf <- read.table("sizefactorsh2o208.1.txt", row.names=1, sep="\t", header=TRUE)
		  sf <- c(t(sf))
		  sizeFactors(ddsnormal) <- sf
		  colnormal = colData(ddsnormal)
		
	#skip step 2 and 3 since it is the same tissue type
	###Input TEtranscripts
		setwd("/DataDrives/dd12/stress/L1HSoutput")
	totaltissue<- read.table("stressL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	ttissue<- t(totaltissue)
	library("stringr")
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
	t_TEcntmatrix = ttissue
	 t_cntmatrix = t(newcountdata)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	colnames(new_cntmatrix) = new_patients
	write.table(new_cntmatrix, file="cntmatrix.h2o208.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t") 
	  
	  
	  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
	  colnames(new_cntmatrix) = NULL
	  ddsnew <- DESeqDataSetFromMatrix(
		countData = new_cntmatrix,
		colData = colnormal,
		design = ~ time)
	  ddsnew$type <- droplevels(ddsnew$type)
	  # set sizefactors 
	  sizeFactors(ddsnew) <- sf
	  save(ddsnew, file="ddsNew.h2o208.DEGES.RData")
	  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
	  write.table(ddsnewcnts, file="normalizedcnt.h2o208.txt", quote=FALSE, row.names=TRUE, sep="\t")
	 
	 
	  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
	  write.table(L1HSnormalized, file="L1HS.normalized.h2o208.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  log2colsums <- log2(colSums(ddsnewcnts))
	  write.table(log2colsums, file="log2colsums.h2o208.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
	#log transformation (vst)
	 vsd <- vst(ddsnew, blind = FALSE)
	  VSTcnts <- assay(vsd)
	  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.h2o208.txt", quote=FALSE, row.names = TRUE)
	  save(VSTcnts, file = "VSTcnts.h2o208.DEGES.RData")
	  write.table(VSTcnts, file="VSTcnt.h2o208.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
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

	  L1HSnormalized = read.table("L1HS.normalized.h2o208.txt", row.names=1, sep="\t", check.names = FALSE)
	  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
	  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
	  colnames(L1HStable)[5:6] = c("normalizedcnt", "VSTcnt")
	  write.table(L1HStable, file="L1HS.h2o208.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
######
#####
###
##TM 0 and 1 test
	newattribfile<- rbind(time0tm, time1tm)
	newcountdata<- cbind(time0tmcut, time1tmcut)

	dds <- DESeqDataSetFromMatrix(countData = newcountdata,
								colData = as.data.frame(newattribfile),
								  design = ~ time)
								  
	# step 1
	  FDR <- 0.1
	  floorPDEG <- 0.05
	  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
	  dds <- dds[filter,]
	  #size factors
	ddsnormal <- estimateSizeFactors(dds)
		sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
		sf <- sizeFactors(ddsnormal) 
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
		write.table(sf, file="sizefactorstm01.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
		sf <- read.table("sizefactorstm01.1.txt", row.names=1, sep="\t", header=TRUE)
		  sf <- c(t(sf))
		  sizeFactors(ddsnormal) <- sf
		  colnormal = colData(ddsnormal)
		
	#skip step 2 and 3 since it is the same tissue type
	###Input TEtranscripts
		setwd("/DataDrives/dd12/stress/L1HSoutput")
	totaltissue<- read.table("stressL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	ttissue<- t(totaltissue)
	library("stringr")
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
	t_TEcntmatrix = ttissue
	 t_cntmatrix = t(newcountdata)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	colnames(new_cntmatrix) = new_patients
	write.table(new_cntmatrix, file="cntmatrix.tm01.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t") 
	  
	  
	  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
	  colnames(new_cntmatrix) = NULL
	  ddsnew <- DESeqDataSetFromMatrix(
		countData = new_cntmatrix,
		colData = colnormal,
		design = ~ time)
	  ddsnew$type <- droplevels(ddsnew$type)
	  # set sizefactors 
	  sizeFactors(ddsnew) <- sf
	  save(ddsnew, file="ddsNew.tm01.DEGES.RData")
	  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
	  write.table(ddsnewcnts, file="normalizedcnt.tm01.txt", quote=FALSE, row.names=TRUE, sep="\t")
	 
	 
	  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
	  write.table(L1HSnormalized, file="L1HS.normalized.tm01.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  log2colsums <- log2(colSums(ddsnewcnts))
	  write.table(log2colsums, file="log2colsums.tm01.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
	#log transformation (vst)
	 vsd <- vst(ddsnew, blind = FALSE)
	  VSTcnts <- assay(vsd)
	  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.tm01.txt", quote=FALSE, row.names = TRUE)
	  save(VSTcnts, file = "VSTcnts.tm01.DEGES.RData")
	  write.table(VSTcnts, file="VSTcnt.tm01.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
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

	  L1HSnormalized = read.table("L1HS.normalized.tm01.txt", row.names=1, sep="\t", check.names = FALSE)
	  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
	  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
	  colnames(L1HStable)[5:6] = c("normalizedcnt", "VSTcnt")
	  write.table(L1HStable, file="L1HS.tm01.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
######
#####
###
##TM 0 and 4 test
	newattribfile<- rbind(time0tm, time4tm)
	newcountdata<- cbind(time0tmcut, time4tmcut)

	dds <- DESeqDataSetFromMatrix(countData = newcountdata,
								colData = as.data.frame(newattribfile),
								  design = ~ time)
								  
	# step 1
	  FDR <- 0.1
	  floorPDEG <- 0.05
	  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
	  dds <- dds[filter,]
	  #size factors
	ddsnormal <- estimateSizeFactors(dds)
		sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
		sf <- sizeFactors(ddsnormal) 
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
		write.table(sf, file="sizefactorstm04.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
		sf <- read.table("sizefactorstm04.1.txt", row.names=1, sep="\t", header=TRUE)
		  sf <- c(t(sf))
		  sizeFactors(ddsnormal) <- sf
		  colnormal = colData(ddsnormal)
		
	#skip step 2 and 3 since it is the same tissue type
	###Input TEtranscripts
		setwd("/DataDrives/dd12/stress/L1HSoutput")
	totaltissue<- read.table("stressL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	ttissue<- t(totaltissue)
	library("stringr")
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
	t_TEcntmatrix = ttissue
	 t_cntmatrix = t(newcountdata)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	colnames(new_cntmatrix) = new_patients
	write.table(new_cntmatrix, file="cntmatrix.tm04.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t") 
	  
	  
	  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
	  colnames(new_cntmatrix) = NULL
	  ddsnew <- DESeqDataSetFromMatrix(
		countData = new_cntmatrix,
		colData = colnormal,
		design = ~ time)
	  ddsnew$type <- droplevels(ddsnew$type)
	  # set sizefactors 
	  sizeFactors(ddsnew) <- sf
	  save(ddsnew, file="ddsNew.tm04.DEGES.RData")
	  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
	  write.table(ddsnewcnts, file="normalizedcnt.tm04.txt", quote=FALSE, row.names=TRUE, sep="\t")
	 
	 
	  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
	  write.table(L1HSnormalized, file="L1HS.normalized.tm04.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  log2colsums <- log2(colSums(ddsnewcnts))
	  write.table(log2colsums, file="log2colsums.tm04.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
	#log transformation (vst)
	 vsd <- vst(ddsnew, blind = FALSE)
	  VSTcnts <- assay(vsd)
	  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.tm04.txt", quote=FALSE, row.names = TRUE)
	  save(VSTcnts, file = "VSTcnts.tm04.DEGES.RData")
	  write.table(VSTcnts, file="VSTcnt.tm04.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
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

	  L1HSnormalized = read.table("L1HS.normalized.tm04.txt", row.names=1, sep="\t", check.names = FALSE)
	  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
	  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
	  colnames(L1HStable)[5:6] = c("normalizedcnt", "VSTcnt")
	  write.table(L1HStable, file="L1HS.tm04.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
######
#####
###
##TM 0 and 8 test
	newattribfile<- rbind(time0tm, time8tm)
	newcountdata<- cbind(time0tmcut, time8tmcut)

	dds <- DESeqDataSetFromMatrix(countData = newcountdata,
								colData = as.data.frame(newattribfile),
								  design = ~ time)
								  
	# step 1
	  FDR <- 0.1
	  floorPDEG <- 0.05
	  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
	  dds <- dds[filter,]
	  #size factors
	ddsnormal <- estimateSizeFactors(dds)
		sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
		sf <- sizeFactors(ddsnormal) 
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
		write.table(sf, file="sizefactorstm08.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
		sf <- read.table("sizefactorstm08.1.txt", row.names=1, sep="\t", header=TRUE)
		  sf <- c(t(sf))
		  sizeFactors(ddsnormal) <- sf
		  colnormal = colData(ddsnormal)
		
	#skip step 2 and 3 since it is the same tissue type
	###Input TEtranscripts
		setwd("/DataDrives/dd12/stress/L1HSoutput")
	totaltissue<- read.table("stressL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	ttissue<- t(totaltissue)
	library("stringr")
	setwd("/DataDrives/dd12/stress/L1HSoutput/timeresults")
	t_TEcntmatrix = ttissue
	 t_cntmatrix = t(newcountdata)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	colnames(new_cntmatrix) = new_patients
	write.table(new_cntmatrix, file="cntmatrix.tm08.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t") 
	  
	  
	  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
	  colnames(new_cntmatrix) = NULL
	  ddsnew <- DESeqDataSetFromMatrix(
		countData = new_cntmatrix,
		colData = colnormal,
		design = ~ time)
	  ddsnew$type <- droplevels(ddsnew$type)
	  # set sizefactors 
	  sizeFactors(ddsnew) <- sf
	  save(ddsnew, file="ddsNew.tm04.DEGES.RData")
	  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
	  write.table(ddsnewcnts, file="normalizedcnt.tm08.txt", quote=FALSE, row.names=TRUE, sep="\t")
	 
	 
	  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
	  write.table(L1HSnormalized, file="L1HS.normalized.tm08.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  log2colsums <- log2(colSums(ddsnewcnts))
	  write.table(log2colsums, file="log2colsums.tm08.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
	#log transformation (vst)
	 vsd <- vst(ddsnew, blind = FALSE)
	  VSTcnts <- assay(vsd)
	  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.tm08.txt", quote=FALSE, row.names = TRUE)
	  save(VSTcnts, file = "VSTcnts.tm08.DEGES.RData")
	  write.table(VSTcnts, file="VSTcnt.tm08.txt", quote=FALSE, row.names=TRUE, sep="\t")
	  
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

	  L1HSnormalized = read.table("L1HS.normalized.tm08.txt", row.names=1, sep="\t", check.names = FALSE)
	  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
	  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
	  colnames(L1HStable)[5:6] = c("normalizedcnt", "VSTcnt")
	  write.table(L1HStable, file="L1HS.tm08.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
		 
	 	 		 
	 	 
	 