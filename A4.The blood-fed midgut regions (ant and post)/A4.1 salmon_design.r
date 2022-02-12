#source("http://bioconductor.org/biocLite.R")
#biocLite()
#a- libraries required for survival analysis and plotting
library(survival)
library(ggplot2)
library(grid)
library(xlsxjars)
library(xlsx)
library(reshape2)
library(multcomp)
library(ARTool)
library(PMCMR)
library(MASS)
library(dplyr)
library(tidyr)
# to use if you prefer csv tables

setwd("e:/my_jianguoyun/friends/cornell_data/bretta_tissue/salmon_out_aaeg")  
#setwd("~/Dropbox/01_Xiaoli Bing/1- data/01 Genome_analysis/mosquito_Breta/RNA_seq/vectorbase/")  
#6- d1 analyses####
d1 <- read.csv("mosquito.Aaeg.salmon_vector.counts.csv", header=T)     
head(d1)[0:3]
d1.1 = setNames(data.frame(t(d1[,-1])), d1[,1])
head(d1.1)[0:4]
d1.1$ID_sample = rownames(d1.1)

d2 = read.csv("design.csv", header=T) 
d2 = dplyr::rename(d2, ID_sample = Barcode)
head(d2)
dim(d2)
#d2$ID_sample = paste0("X",d2$ID,"_ID")
d2$condition <- paste0(d2$Tissue,"_",d2$Treatment)

d2.1 = merge(d2, d1.1, by = "ID_sample")
head(d2.1) [1:10]


a1 <- d2.1

colnames(a1)[0:16]
c1 <- a1[,1:8]
b1 = dplyr::select(a1, ID_sample, 9:length(colnames(a1)))

b1.1 = setNames(data.frame(t(b1[,-1])), b1[,1])

head(b1.1)[1:4]
colnames(c1)
dim(b1)

#c3 needsto delete twhe first row mannually
write.csv(b1.1, "mosquito_salmon_reads.csv", row.names = T)
write.csv(c1, "mosquito_salmon_colData.csv", row.names = F)


head(b1.1)
a1.1 = colSums(b1.1)
#calculate TPMs for each sample
a1.2 = 1000000 * sweep(b1.1, 2, colSums(b1.1),`/`)
head(a1.2)
colSums(a1.2)
write.csv(a1.2, "mosquito.Aaeg.salmon_TPMs.csv", row.names = T)




