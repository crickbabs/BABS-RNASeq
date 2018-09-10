#!/usr/bin/Rscript

# nourdine.bah@crick.ac.uk
# philip.east@crick.ac

library(optparse)
library(rtracklayer)
library(SummarizedExperiment)

############
## ARGUMENTS

# the directory that contains the counts files
results.opt <-
	make_option(
					c("-r", "--results-dir"),
					type="character",
					default=NULL,
					help="directory that contains the count files"
					)

# the design file that maps sample names to conditions
design.opt <-
	make_option(
					c("-d", "--design-file"),
					type="character",
					default=NULL,
					help="design file absolute path"
					)

# the gene annotations
gtf.opt <-
	make_option(
					c("-g", "--gtf-file"),
					type="character",
					default=NULL,
					help="gtf file absolute path"
					)

# the design file that maps sample names to conditions
output.opt <-
	make_option(
					c("-o", "--output-name"),
					type="character",
					default=NULL,
					help="output basename for csv and rda"
					)

# create the options and parse them
options <- list(results.opt, design.opt, gtf.opt, output.opt)
parser <- OptionParser(option_list=options)
args <- parse_args(parser)
directory <- args[["results-dir"]]
design_filepath <- args[["design-file"]]
gtf_filepath <- args[["gtf-file"]]
name <- args[["output-name"]]
############
############

# the design file
design <- read.csv(design_filepath)

# annotation
gtf <- rtracklayer::import(gtf_filepath)

# create the matrix
files <- sort(list.files( path=directory , full.names=T, pattern="*.tsv$" ))
countsL <- lapply(files, read.table, header=T, row.names=1)
counts <- as.matrix( do.call(cbind, countsL) )
#counts <- counts[ rowMeans(counts) > 0 , ]
coldata <- design[ match( colnames(counts) , design[,"sample"] ) , ]
rowdata <- gtf[match(rownames(counts), S4Vectors::mcols(gtf)[,"gene_id"] ),]
se <- SummarizedExperiment::SummarizedExperiment(
																 assays=list(counts=counts),
																 colData=coldata,
																 rowData=rowdata
																 )

# export rda
data <- se
rda <- paste0( name , ".rda" )
save(data, file=rda)

# export csv
csv <- paste0( name , ".csv" )
write.csv( as.data.frame(SummarizedExperiment::assay( se )) , csv )

