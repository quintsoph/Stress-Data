###Gather data for L1HS of each pat and put into one file
#Test with one file
setwd("/DataDrives/dd12/stress/TETranscripts")
pat<- read.table("Total_RNA_H2O2_Rep1_0hr_7008199/Total_RNA_H2O2_Rep1_0hr_7008199.discount.ele.cntTable")
rownames(pat)=pat[,1]
L1HSpat<-pat[c("L1HS:L1:LINE"),]
final<-cbind(L1HSpat, "Total_RNA_H2O2_Rep1_0hr_7008199")
colnames(final)=c("transposon", "expression", "patient")

#create for all tissue pats
#stress data
data.frame(patientlist<- read.table("/DataDrives/dd12/stress/TETranscripts/stresspats.txt"))
setwd("/DataDrives/dd12/stress/TETranscripts")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/DataDrives/dd12/stress/L1HSoutput")
write.table(allpats, "stressL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
