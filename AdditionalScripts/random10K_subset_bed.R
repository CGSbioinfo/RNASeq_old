#!/usr/local/bin/Rscript

dir = commandArgs(TRUE)[1]
file = commandArgs(TRUE)[2]

setwd(dir)

x=read.table(file, stringsAsFactors=FALSE)
l=dim(x)[1]
lines=sample(1:l,10000)
name=strsplit(file,".bed")[1]

out=x[lines,]
write.table(out,paste0(dir,"/",name,"_10k.bed"),quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
