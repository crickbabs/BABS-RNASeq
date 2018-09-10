#!/usr/bin/Rscript

library(optparse)
library(DESeq2)
library(jsonlite)

####################
### SCRIPT ARGUMENTS

rda.opt <- make_option( c("-r", "--rda"), type="character", default=NULL )
des.opt <- make_option( c("-d", "--design"), type="character", default=NULL )
options <- list(rda.opt, des.opt)
parser <- OptionParser(option_list=options)
args <- parse_args(parser)
rda <- args[["rda"]]
design_path <- args[["design"]]

##########
## LOADING

# pca
load(rda)
pca <- data

# design file
design <- read.csv(design_path)

##############
## PLOTLY INFO

# experimental conditions
design <- design[ , - grep( "^file$|^file1$|^file2$" , names(design) ) ]
conditions <- names(design)[ - grep( "^sample$", names(design) ) ]
conditions <- sort(conditions)

# the output data frame
df <- merge(pca, design, by.x="name", by.y="sample")
df <- df[ , c("name", "PC1", "PC2", conditions) ]

# at the moment, the plot is a scatter plot with points whose features are
# color, shape, line color. So only three experimental variables can be
# represented

# return a color for each value of a factor
get_colors <-
	function(x) {

		# generate colors for each value of the experimental condition
		values <- as.numeric(as.factor( x ))
		m <- length(unique( values ))
		colors <- paste0( "#" , toupper(as.hexmode( sample( 2^24 , m ) )) )

		# match the colors to their value
		map <- colors[ values ]

		return(map)
	}

# map each factor to its marker feature
map <- data.frame("color"="", "symbol"="", "line"="")

# if the user haven't provided experimental conditions
if ( length(conditions) > 0 ) {

	# add marker features for the first three experimental variables
	n <- min(3, length(conditions))
	for (i in 1:n) {
	
		if (i==1) {
			df[,"color"] <- get_colors( df[ , conditions[i] ] )
			map[,"color"] <- conditions[i]
		}
	
		else if (i==2) {
			df[,"symbol"] <- as.numeric(as.factor( df[ , conditions[i] ] ))
			map[,"symbol"] <- conditions[i]
		}
	
		else if (i==3) {
			df[,"line"] <- get_colors( df[ , conditions[i] ] )
			map[,"line"] <- conditions[i]
		}
	}
	
	# add hover text
	if ( length(conditions) == 1 ) {
		df[,"hover"] <- df[,conditions]
	} else {
		df[,"hover"] <- apply(df[,conditions], 1 ,paste, collapse = "-")
	}

	df[,"hover"] <- paste(df[,"name"], df[,"hover"], sep=":")

} else {

	# just set the hover text
	df[,"hover"] <- df[,"name"]

}

# percentage of variance
pvar <- round( attr(pca, "percentVar")*100 , 1 )

#########
## OUTPUT

# export to json
l <- list( percent_var=pvar , map=map , pca=df )
json <- toJSON(l, pretty=T)
print(json)

