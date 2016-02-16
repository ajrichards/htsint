#!/usr/bin/Rscript
## install notes
# > source("http://bioconductor.org/biocLite.R")
# > biocLite("DESeq2")

## samples notes
#    A = endurant
#    B = non-endurant
#    C = non-endurant
#    D = endurant
#    E = endurant
#    F = non-endurant
#    G = endurant
#    H = non-endurant

library("DESeq2")

args <-commandArgs(TRUE)

if (length(args)==0){
    print("ERROR: Did not specify counts file e.g 'Rscript run-deseq.R raw-counts.csv'") 
    q()
}

countsFile = args[1]
outFile = args[2]

if (!file.exists(countsFile)){
    print("ERROR: invalid counts file")    
    q()
}

if (is.na(outFile)){
    outFile = "deseq.csv"
}

#### from counts file (count matrix)
countData <- read.csv(countsFile,header=TRUE,row.names=1,com='')
countData <- round(countData)
colData <- data.frame(
    row.names=colnames(countData),
    condition=c("endurant", "non-endurant", "non-endurant", "endurant", "endurant", "non-endurant","endurant","non-endurant"),
    libType=c("paired", "paired", "paired","paired","paired","paired","paired","paired")
)

## filter the data
dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
dds$condition <- factor(dds$condition,levels=c("non-endurant","endurant"))

dds <- DESeq(dds)
rld <- rlog(dds, blind=FALSE)
res <- results(dds)
use <- res$baseMean>attr(res,"filterThreshold")
print(table(use))

## run DESeq2 using filtered data
countData <- countData[use& ! is.na(res$pvalue),]

## create a DESeqDataSet
colData <- data.frame(
    row.names=colnames(countData),
    condition=c("endurant", "non-endurant", "non-endurant", "endurant", "endurant", "non-endurant","endurant","non-endurant"),
    libType=c("paired", "paired", "paired","paired","paired","paired","paired","paired")
)

dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
dds$condition <- factor(dds$condition,levels=c("non-endurant","endurant"))

## set 'non-endurant' as the control
dds <- DESeq(dds)
rld <- rlog(dds, blind=FALSE)
res <- results(dds)
resOrdered <-res[order(res$padj),]
print(head(resOrdered))

## export results to csv
write.csv(as.data.frame(resOrdered),file=outFile,quote=FALSE)
write.csv(as.data.frame(assay(rld)),file=gsub("\\.csv","-samples.csv",outFile),quote=FALSE)
