#!/usr/local/bin/Rscript

suppressMessages(library(RColorBrewer))
suppressMessages(library(matrixStats))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))

dir = commandArgs(TRUE)[1]

##########################
# HTSEQ COUNT OUTPUT
##########################

setwd('./countedReads/1pass_ann/')
files=list.files()
files=files[grep(pattern = "count",files)]

data=read.table(files[1], row.names=1)
#colnames(data)=strsplit(gsub("-","_",files[1]),"_\\.")[[1]][1]
colnames(data)=files[1]

for (i in 2:length(files)){
  name=files[i]
  name2=gsub("-","_",name)
  name2=strsplit(name2,"_\\.")[[1]][1]
  assign(name, read.table(name))
  sample=eval(as.name(name))
  index=match(rownames(data), sample[,1])
  data=cbind(data, sample[index,2])
  colnames(data)[i]=name
}

data.notcounted=data[which(rownames(data)=='__no_feature'):dim(data)[1],]
colnames(data.notcounted)=gsub('.count','',colnames(data.notcounted))

data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),]
colnames(data)=gsub('.count','',colnames(data))

#ord=c()
#for (i in 1:ncol(data)){
#  ord=c(ord,grep(paste0("^",i,"_"),colnames(data)))
#}

#data=data[,ord]

# Distribution of counts
n_reads_in_genes=colSums(data)

n_reads_total=read.csv('../../mapping_summary.csv', row.names=1)
sample_names=rownames(n_reads_total)
n_reads_total=n_reads_total$Mapped_num
names(n_reads_total)=sample_names

n_reads_no_feature=t(data.notcounted["__no_feature",])[,1]
n_reads_ambiguous=t(data.notcounted["__ambiguous",])[,1]

counting_summ=cbind(n_reads_total, n_reads_in_genes[sample_names], n_reads_no_feature[sample_names], n_reads_ambiguous[sample_names])
colnames(counting_summ)=c('Total','In_genes', 'Not_in_genes', 'Ambiguous')

# Numbers of genes detected
n_genes=sapply(data,function(x){table(x==0)['FALSE']})
names(n_genes)=gsub('.FALSE','',names(n_genes))
counting_summ=cbind(counting_summ,Num_genes=n_genes)
write.csv(counting_summ,'../../counts_summary.csv')

# Plot
data.melt=suppressMessages(melt(data))
setwd("../../")
pdf("countsDistributionHist.pdf", width=12, height=7)
#par(mar=c(6.1,4.1,4.1,2.1))
#suppressMessages(ggplot(data.melt, aes(log(value), fill = variable)) +  geom_density(alpha=0.3))
ggplot(data.melt, aes(x=log(value), fill=variable )) + stat_bin(biwidth=1) + theme(legend.text=element_text(size=12))
suppressMessages(dev.off())
pdf("countsDistributionBox1.pdf", width=12, height=7)
par(mar=c(9.7,4.1,4.1,2.1))
boxplot(data,las=2)
suppressMessages(dev.off())
pdf("countsDistributionBox2.pdf", width=12, height=7)
par(mar=c(9.7,4.1,4.1,2.1))
boxplot(data,las=2, outline=FALSE)
suppressMessages(dev.off())

### heatmap correlation samples based on count profile
data_norm=t(t(data)/rowSums(t(data)))
pdf('sample_counts_correlation_heatmap.pdf')
heatmap.2(cor(data_norm), scale=c('none'), margins=c(10,10),density.info='density', trace='none')
dev.off()




