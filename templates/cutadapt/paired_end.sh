#!/bin/sh

cutadapt \
	-a AGATCGGAAGAGC \
	-A AGATCGGAAGAGC \
	-o $output1 \
	-p $output2 \
	-e 0.1 \
	-q 10 \
	-m 25 \
	-O 1 \
	$fastq1 $fastq2 > $logfile

