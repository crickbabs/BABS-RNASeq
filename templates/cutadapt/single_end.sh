#!/bin/sh

cutadapt \
	-a AGATCGGAAGAGC \
	-o $output \
	-e 0.1 \
	-q 10 \
	-m 25 \
	-O 1 \
	$fastq > $logfile
