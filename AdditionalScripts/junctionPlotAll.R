#!/usr/local/bin/Rscript

library(reshape)
library(grid)
library(ggplot2)

dir = commandArgs(TRUE)[1]
setwd(dir)


y=read.csv('alignedReads/1pass_ann/QC/junctionSat_all.csv', row.names=1, header=FALSE)

y=cbind(rownames(y),y/1000)
colnames(y)=x=c('sample',seq(from=5,to=100,by=5))
y=melt(y)

pdf("alignedReads/1pass_ann/QC/junctionSaturationAll.pdf", width=12)
ggplot(y, aes(variable, value, group=sample, colour = sample)) + geom_line() + ylab("Number of splicing junctions (x1000)") + xlab("Percent of total reads") + theme( axis.title.x =element_text(size=14), axis.title.y =element_text(size=14), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key.height=unit(.8,"line")) 
dev.off()
