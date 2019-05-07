#gather genes from the file and format for DESEQ
setwd("/DataDrives/dd12/stress/TETranscripts")
trials<- read.table("stresspats.txt", sep="\t")
genefilematrix=NULL
trialcut2=NULL
for (i in trials[,1]){
	trialfile=NULL
	trialfile= paste(i,"/",i,".ele.cntTable", sep="")
	trialoutput<- read.table(trialfile)
	trialcut<- data.frame(trialoutput[2:27968,])
	trialcut2<- data.frame(trialcut[,2])
	trialcutgenes<- data.frame(trialcut[,1])
	rownames(trialcut2)=trialcutgenes[,1]
	colnames(trialcut2)=i
	trialcut2<- as.matrix(trialcut2)
	genefilematrix<- cbind(genefilematrix, trialcut2)
	}

	setwd("/DataDrives/dd12/stress/L1HSoutput")
	write.table(genefilematrix, "stressgenematrix.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#create attribute file
setwd("/DataDrives/dd12/stress/TETranscripts")
trials<- read.table("stresspats.txt", sep="\t")
test<- as.matrix(trials)
stringers<- strsplit(test, "_")
finaltreatment=NULL
for (i in 1:36){

	treatment<- stringers[[i]][3]
	finaltreatment<- rbind(finaltreatment, treatment)
}
finalrep=NULL
for (i in 1:36){

	repeation<- stringers[[i]][4]
	finalrep<- rbind(finalrep, repeation)
}
finaltime=NULL
for (i in 1:36){

	times<- stringers[[i]][5]
	finaltime<- rbind(finaltime, times)
}

atributes<- cbind(trials, finaltreatment, finalrep, finaltime)
colnames(atributes)=c("trial", "treatment", "rep", "time")
	setwd("/DataDrives/dd12/stress/L1HSoutput")
	write.table(atributes, "attributematrix.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

