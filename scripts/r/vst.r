#!/usr/bin/Rscript

library(DESeq2)

# script arguments
args <- commandArgs(trailingOnly=T)
rda_to_load <- args[1]
rda_to_save <- args[2]

# loading
load(rda_to_load)
dds <- data

# variance stabilising transformation
if ( nrow(dds) < 1000 ) {
	#vst <- DESeq2::vst( dds , nsub=nrow(dds) )
	vst <- DESeq2::varianceStabilizingTransformation( dds )
} else {
	#vst <- DESeq2::vst(dds)
	vst <- DESeq2::varianceStabilizingTransformation( dds )
}

# saving
data <- vst
save(data, file=rda_to_save)

