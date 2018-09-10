#!/usr/bin/Rscript

library(DESeq2)

# script arguments
args <- commandArgs(trailingOnly=T)
rda_to_load <- args[1]
rda_to_save <- args[2]

# loading
load(rda_to_load)
dds <- data

# regularised logarithm
rlog <- DESeq2::rlog(dds)

# saving
data <- rlog
save(data, file=rda_to_save)

