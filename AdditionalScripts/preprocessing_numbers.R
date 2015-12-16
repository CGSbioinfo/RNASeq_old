#!/usr/local/bin/Rscript

dir = commandArgs(TRUE)[1]
setwd(dir)
system("egrep -w \"Total Sequences|Sequence length\" */*/*/fastqc_data* > preprocessing_numbers.txt")

x=read.table("preprocessing_numbers.txt",sep="\t")
if(!length(grep("trimm",x[,1]))==0){
 x=x[-grep("trimm",x[,1]),] 
}
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

names=gsub("all_lane_","",rownames(m))
names=gsub("_fastqc","",rownames(m))
rownames(m)=names
colnames(m)=c("Total No of Reads", "Sequence Lenght")
write.csv(m,"preprocessing_numbers.csv", quote=F)
system("rm preprocessing_numbers.txt")

