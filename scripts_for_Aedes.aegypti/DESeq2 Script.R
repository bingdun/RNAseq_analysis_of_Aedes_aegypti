#a - load DESeq2 pacakge

library(DESeq2)

#b- set a directory to work in
setwd("") #active folder

#c- open and read your data
cts <- read.table("Read Data.csv", header=TRUE, sep=",", dec=".",na.strings=".",row.names = 1) # call the data
coldata <- read.table("Column Data.csv", header=TRUE, sep=",", dec=".",na.strings=".")

print(rownames(coldata))

row.names(coldata)
row.names(coldata) <- coldata$batch




all(rownames(coldata) %in% colnames(cts))



dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds<-DESeq(dds)
res<-results(dds)

res
write.csv(res, file ="X Results.csv")

