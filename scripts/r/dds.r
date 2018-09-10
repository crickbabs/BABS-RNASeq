#!/usr/bin/Rscript

library(optparse)
library(DESeq2)
library(BiocParallel)

###############################################################################
### SCRIPT ARGUMENTS

# arguments
rda.opt <- make_option( c("-r", "--rda-se"), type="character", default=NULL )
fml.opt <- make_option( c("-f", "--formula"), type="character", default="~ 1" )
cor.opt <- make_option( c("-c", "--ncores"), type="integer", default=30 )
out.opt <- make_option( c("-o", "--rda-out"), type="character", default=NULL )

# parsing
options <- list(rda.opt, fml.opt, cor.opt, out.opt)
parser <- OptionParser(option_list=options)
args <- parse_args(parser)

# affectation
rda_se <- args[["rda-se"]]
formula <- args[["formula"]]
ncores <- args[["ncores"]]
rda_out <- args[["rda-out"]]

###############################################################################
### PARALLELIZATION

register(MulticoreParam(ncores))

###############################################################################
### LOADING DATA

load(rda_se)
se <- data

###############################################################################
### DESEQ2 MODELING

SummarizedExperiment::assay(se) <- round( SummarizedExperiment::assay(se) )
dds <- DESeq2::DESeqDataSet(se, design=as.formula(formula))
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, parallel=T)

###############################################################################
### SAVING

data <- dds
save(data, file=rda_out)

