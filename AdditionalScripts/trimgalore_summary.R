#!/usr/local/bin/Rscript

dir = commandArgs(TRUE)[1]
setwd(dir)
system("egrep -w \"Total Sequences|Sequence length\" */*/*/fastqc_data* > preprocessing_numbers.txt")

x=read.table("preprocessing_numbers.txt",sep="\t")
#x=x[-grep("trimm",x[,1]),]
names=strsplit(as.character(x[,1]),"/")
names=unique(sapply(names,"[",c(3)))
total_seq=x[grep("Total Sequences",x[,1]),]
m=matrix(NA,ncol=2,nrow=length(names))
rownames(m)=names
for (i in 1:length(names)){
  m[names[i],1] = as.character(total_seq[grep(names[i],total_seq[,1]),2])
}
seq_length=x[grep("Sequence length",x[,1]),]
for (i in 1:length(names)){
  m[names[i],2] = as.character(seq_length[grep(names[i],seq_length[,1]),2])
}
names=gsub("_all_lane","",rownames(m))
names=gsub("_fastqc","",names)
rownames(m)=names
colnames(m)=c("Total No of Reads", "Sequence Lenght")
system("rm preprocessing_numbers.txt")

#After trimming
system("egrep -w \"Total Sequences|Sequence length\" trimmedReads/*/fastqc_data* > trimgalore_summary.txt")

x=read.table("trimgalore_summary.txt",sep="\t")
names=strsplit(as.character(x[,1]),"/")
names=unique(sapply(names,"[",c(2)))
total_seq=x[grep("Total Sequences",x[,1]),]
m2=matrix(NA,ncol=2,nrow=length(names))
rownames(m2)=names
for (i in 1:length(names)){
  m2[names[i],1] = as.character(total_seq[grep(names[i],total_seq[,1]),2])
}
seq_length=x[grep("Sequence length",x[,1]),]
for (i in 1:length(names)){
  m2[names[i],2] = as.character(seq_length[grep(names[i],seq_length[,1]),2])
}
names=gsub("all_lane_","",rownames(m2))
names=gsub("_fastqc","",names)
names=gsub("_val_(\\d)","",names)
rownames(m2)=names
colnames(m2)=c("Total No of Reads", "Sequence Lenght")
system("rm trimgalore_summary.txt")

out=cbind(m2[,1],round(as.numeric(m2[,1])/as.numeric(m[,1]),digits=4),m2[,2])
colnames(out)=c("Total No of Reads","Percentage kept","Sequence Length")

write.csv(out,"trimmgalore_summary.csv",quote=F)

