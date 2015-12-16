#!/usr/local/bin/Rscript

suppressMessages(library(RColorBrewer))
suppressMessages(library(matrixStats))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
require(grid)

dir = commandArgs(TRUE)[1]

##########################
# HTSEQ COUNT OUTPUT
##########################

setwd('./countedReads/1pass_ann/')
files=list.files()
files=files[grep(pattern = "count",files)]

data=read.table(files[1], row.names=1)
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
colnames(data.melt)[1]='Sample'
ggplot(data.melt, aes(x=log(value), fill=Sample )) + stat_bin(biwidth=1) + theme(legend.text=element_text(size=10), legend.key.size = unit(.45, "cm"))
suppressMessages(dev.off())

# Distribution lines
distribution=as.data.frame(matrix(0,ncol=length(files),nrow=14))
colnames(distribution)=colnames(data)
rownames(distribution)=c('0', '1-100', '101-200','200-500','501-1,000', '1,001-5,000', '5,001-10,000', '10,001-20,000', '20,001-30,000', '30,000-40,000', '40,001-50,000', '50,001-100,000', '100,001-50,000', '>500,000')
for (i in 1:length(colnames(data))) {
  sample_name=colnames(data)[i]
  x=data[,sample_name]
  zero=table(x==0)['TRUE']
  bin1=table(x>0 & x<=100)['TRUE']
  bin2=table(x>100 & x<=200)['TRUE']
  bin3=table(x>200 & x<=500)['TRUE']
  bin4=table(x>500 & x<=1000)['TRUE']
  bin5=table(x>1000 & x<=5000)['TRUE']
  bin6=table(x>5000 & x<=10000)['TRUE']
  bin7=table(x>10000 & x<=20000 )['TRUE']
  bin8=table(x>20000 & x<=30000 )['TRUE']
  bin9=table(x>30000 & x<=40000 )['TRUE']
  bin10=table(x>40000 & x<=50000 )['TRUE']
  bin11=table(x>50000 & x<=100000 )['TRUE']
  bin12=table(x>100000 & x<=500000 )['TRUE']
  bin13=table(x>500000)['TRUE']
  sample=c(zero,bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13)
  sample[is.na(sample)]=0
  distribution[,sample_name]=sample
}
distribution=distribution/colSums(distribution)*100
distribution.melt=melt(cbind(bins=factor(row.names(distribution), levels=row.names(distribution)), distribution))
pdf('countDistribution_lines.pdf', width=12)
ggplot(distribution.melt,aes(x=bins,y=value,color=variable,group=variable)) + geom_path() + 
  xlab('Number of counts') + ylab('Percentage of genes') + 
  scale_y_continuous(breaks=seq(0,70,10), labels=c('0%','10%','20%','30%','40%','50%','60%','70%')) + 
  scale_colour_discrete(name='Sample') + 
  theme(axis.text.x=element_text(angle=90), axis.text.y=element_text(),legend.text=element_text(size=10), legend.key.size = unit(.38, "cm"))
dev.off()

# Distribution lines removing all zeros
data_filtered=data[!rowSums(data)==0,]
distribution=as.data.frame(matrix(0,ncol=length(files),nrow=14))
colnames(distribution)=colnames(data_filtered)
rownames(distribution)=c('0', '1-100', '101-200','200-500','501-1,000', '1,001-5,000', '5,001-10,000', '10,001-20,000', '20,001-30,000', '30,000-40,000', '40,001-50,000', '50,001-100,000', '100,001-50,000', '>500,000')
for (i in 1:length(colnames(data_filtered))) {
  sample_name=colnames(data_filtered)[i]
  x=data_filtered[,sample_name]
  zero=table(x==0)['TRUE']
  bin1=table(x>0 & x<=100)['TRUE']
  bin2=table(x>100 & x<=200)['TRUE']
  bin3=table(x>200 & x<=500)['TRUE']
  bin4=table(x>500 & x<=1000)['TRUE']
  bin5=table(x>1000 & x<=5000)['TRUE']
  bin6=table(x>5000 & x<=10000)['TRUE']
  bin7=table(x>10000 & x<=20000 )['TRUE']
  bin8=table(x>20000 & x<=30000 )['TRUE']
  bin9=table(x>30000 & x<=40000 )['TRUE']
  bin10=table(x>40000 & x<=50000 )['TRUE']
  bin11=table(x>50000 & x<=100000 )['TRUE']
  bin12=table(x>100000 & x<=500000 )['TRUE']
  bin13=table(x>500000)['TRUE']
  sample=c(zero,bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13)
  sample[is.na(sample)]=0
  distribution[,sample_name]=sample
}
distribution=distribution/colSums(distribution)*100
distribution.melt=melt(cbind(bins=factor(row.names(distribution), levels=row.names(distribution)), distribution))
pdf('countDistribution_lines_nozeros.pdf', width=12)
ggplot(distribution.melt,aes(x=bins,y=value,color=variable,group=variable)) + geom_path() + 
  xlab('Number of counts') + ylab('Percentage of genes') + 
  scale_y_continuous(breaks=seq(0,70,10), labels=c('0%','10%','20%','30%','40%','50%','60%','70%')) + 
  scale_colour_discrete(name='Sample') + 
  theme(axis.text.x=element_text(angle=90), axis.text.y=element_text(),legend.text=element_text(size=10), legend.key.size = unit(.38, "cm"))
dev.off()


### heatmap correlation samples based on count profile
data_norm=t(t(data)/rowSums(t(data)))
pdf('sample_counts_correlation_heatmap.pdf')
heatmap.2(cor(data_norm), scale=c('none'), margins=c(10,10),density.info='density', trace='none')
dev.off()




