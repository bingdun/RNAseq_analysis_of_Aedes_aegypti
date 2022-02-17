
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("DESeq2")
#biocLite("pasilla")
#install.packages("pheatmap")
library("DESeq2")
#library("edgeR")
library(splitstackshape)
#library(factoextra)
library(pheatmap)
library(ggfortify)
library(ggplot2)
library(ggrepel)
library(grid)
library(reshape2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

setwd("~/salmon_out_aaeg")  

remove(list = ls())

a1<-read.csv("mosquito_salmon_reads.csv", row.names=1)

a1.1 = round(a1, 0)
countData <- as.matrix(a1.1)
head(countData)
colData <- read.csv("mosquito_salmon_colData.csv", row.names=1)
head(colData)
ncol(countData)
colnames(colData)
#colnames(countData) <- NULL

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

dds$condition
all(rownames(colData) %in% colnames(countData))
all(rownames(colData) == colnames(countData))

#dds <- estimateSizeFactors(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds


#Normalization: estimate the effective library size.
dds_sub = dds
#dds_sub = dds
#dds_sub = dds[, grepl("Tissue*", dds$condition ) ]

dds_sub$condition <- droplevels(dds_sub$condition)
dds_sub
dds_sub = estimateSizeFactors(dds_sub) 
#This function obtains dispersion estimates for Negative Binomial distributed data.
cdsBlind = estimateDispersions(dds_sub)  #It doesn't change PCA results

#These two tranformations are similar, the rlog might perform a bit better when the size factors vary widely, and the varianceStabilizingTransformation is much faster when there are many samples.
#vsd= rlog(cdsBlind)
vsd= varianceStabilizingTransformation(cdsBlind)

#rowVars: Variance Estimates For Each Row (Column) In A Matrix
rv = rowVars(assay(vsd))
#select the first 500 most variant genes 
select = order(rv, decreasing=TRUE)[seq_len(500)]
#select = order(rv, decreasing=TRUE)

pca = prcomp(t(assay(vsd)[select,]))



vsd <- varianceStabilizingTransformation(dds_sub) 
graph <- plotPCA(vsd, intgroup=c("Treatment","Tissue"), returnData=TRUE)
percentVar <- round(100 * attr(graph, "percentVar"))
cols = brewer.pal(12, "Paired")
list(color = brewer.pal(12, "Paired"))
levels(graph$Treatment)
#plot based on bacteria
ggplot(graph, aes(PC1, PC2, color=Tissue,shape= as.character(Tissue))) +
  geom_point(size=3) +
  #scale_color_manual(values = cols) +
  #guides(fill=FALSE, shape= FALSE, colour=FALSE)+
  #stat_ellipse() + 
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  #scale_colour_brewer(palette = "Palette_Name")+
  theme_set(theme_bw())+
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave("mosquito_all_tissue_PCA.png",width = 8, height = 7)
ggsave("mosquito_all_tissue_PCA.pdf",width = 8, height = 7)


