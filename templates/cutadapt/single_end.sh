#!/bin/sh

cutadapt \
	-a AGATCGGAAGAGC \
	-o $output \
	-e 0.1 \
	-q 10 \
	-m 25 \
	-O 1 \
	$fastq > $logfile

md5sum $fastq > md5.txt
