#source("http://bioconductor.org/biocLite.R")
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
#biocLite("devtools")    # only if devtools not yet installed
#biocLite("COMBINE-lab/wasabi")
#biocLite("pachterlab/sleuth")

library(wasabi)
library(sleuth)
library(tximport)
library(readr)


#base_dir <- "e:/Cornell/mosquito_Breta/RNA_seq/20180915_rnaseq/"
base_dir <- "e:/my_jianguoyun/friends/cornell_data/bretta_tissue/salmon_out_aaeg"


#base_dir <- "/home/xiaoli/Documents/bretta_rnaseq/kallisto"
#set the number of CPU’s to use for the analysis; set it according to your instance size
#options(mc.cores = 4L)
#he first step in a sleuth analysis is to specify where the kallisto results are stored
sample_id <- dir(file.path(base_dir,"salmon_output"))
sample_id
#A list of paths to the kallisto results indexed by the sample IDs
sal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "salmon_output",id))
sal_dirs

#this step is from "wasabi"
prepare_fish_for_sleuth(sal_dirs)

#read file path and define path
b1 = read.csv(file.path(base_dir, "list.txt"), sep = ".", header = TRUE, stringsAsFactors=FALSE)
b1$X00
b1.1 <- dplyr::mutate(b1, path = sal_dirs)
b1.2 = dplyr::select(b1.1, X00, path)


files <- file.path(sal_dirs, "quant.sf")
#names(files) <- paste0("sample", 1:6)
#all(file.exists(files))
a0 = read.csv(file.path(base_dir, "t2g.csv"), header = TRUE, stringsAsFactors=FALSE)
head(a0)
tx2gene = dplyr::select(a0, TXNAME = transcript, GENEID = gene)
head(tx2gene)

#assign count from trascripts to genes
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$abundance)[,1:3]
head(txi$counts)[,1:3]

#export counts per gene
b2 = txi$counts
colnames(b2) = b1$X00
head(b2)
getwd()
setwd("e:/my_jianguoyun/friends/cornell_data/bretta_tissue/salmon_out_aaeg") 

write.csv(b2, "mosquito.Aaeg.salmon_vector.counts.csv")


