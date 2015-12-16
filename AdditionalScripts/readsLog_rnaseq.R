#!/usr/local/bin/Rscript

dir = commandArgs(TRUE)[1]
setwd(dir)

tg=read.table("tg_Logout.txt",sep="\t", stringsAsFactors = FALSE)
tg[,2]=gsub("_","",tg[,2])
tg[,1]=gsub("_all.*","",tg[,1])

procreads=tg[tg$V2=="Processedreads:",]
procreads[,3]=gsub("_","",procreads[,3])
procreads=procreads[seq(from=1,to=(dim(procreads)[1]),by=2),]
rownames(procreads)=procreads[,1]

trimreads=tg[tg$V2=="Trimmedreads:",]
trimreads[,3]=gsub("_","",trimreads[,3])
trimreads[,3]=gsub("\\("," ",trimreads[,3])
trimreads[,3]=gsub("\\%)"," ",trimreads[,3])
trimreads=cbind(trimreads[,1:2],matrix(unlist(strsplit(trimreads[,3], " ")), ncol=2, byrow = TRUE))
trimreads[,1]=paste(trimreads[,1],rep(c("R1","R2"),(dim(trimreads)[1]/2)),sep="_")

star=read.table("star_logout.txt")
star[,1]=gsub("_all_lane_Log.final.out:","",star[,1])
star[,2]=gsub("__","",star[,2])
star[,2]=gsub("_ |","",star[,2])
star[,3]=gsub("%","",star[,3])
star=split(star,star$V2)

out=as.data.frame(cbind(Input_Reads=procreads[,3],After_Trimming=star[[1]][match(rownames(procreads), star[[1]][,1]),3],
                        Percentage_after_trimming=round(as.numeric(star[[1]][match(rownames(procreads), star[[1]][,1]),3])/as.numeric(procreads[,3])*100, digits=2),
                        Percentage_Uniquely_Mapped_Reads=as.numeric(star[[9]][match(rownames(procreads), star[[9]][,1]),3]),
                        Percentage_Multiple_Loci_Reads=as.numeric(star[[4]][match(rownames(procreads), star[[4]][,1]),3])+
                          as.numeric(star[[5]][match(rownames(procreads), star[[5]][,1]),3]),
                        Percentage_Unmapped_Reads=as.numeric(star[[6]][match(rownames(procreads), star[[6]][,1]),3])+
                          as.numeric(star[[7]][match(rownames(procreads), star[[7]][,1]),3])+
                          as.numeric(star[[8]][match(rownames(procreads), star[[8]][,1]),3])))
rownames(out)=rownames(procreads)
write.csv(out,file="reads_log.csv",row.names=TRUE,quote=FALSE)
