#!/bin/sh

java -Xmx10g -jar \$EBROOTRNAMINSEQC/RNA-SeQC_v[0-9].[0-9].[0-9].jar \
	-d 1000000 \
	-rRNA $rrna_list \
	-r $fasta \
	-t $rnaseqc_gtf \
	-o $metrics_filename \
	-singleEnd \
	-s "$name|$bam|$name"

