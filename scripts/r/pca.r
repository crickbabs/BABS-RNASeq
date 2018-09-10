#!/usr/bin/Rscript

library(DESeq2)

# script arguments
args <- commandArgs(trailingOnly=T)
rda_to_load <- args[1]
rda_to_save <- args[2]

# loading
load(rda_to_load)
rlog <- data

# let DESeq do the PCA for us
pca <- DESeq2::plotPCA(rlog, intgroup="sample", ntop=1000, returnData=T)

# saving
data <- pca
save(data, file=rda_to_save)

