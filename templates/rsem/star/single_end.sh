#!/bin/sh

rsem-calculate-expression \
	--temporary-folder "tmp" \
	--star \
	--num-threads $cpus \
	--strandedness $strandedness \
	--estimate-rspd \
	--seed 1 \
	--star-output-genome-bam \
	--star-gzipped-read-file \
	$fastq \
	$index \
	$name
