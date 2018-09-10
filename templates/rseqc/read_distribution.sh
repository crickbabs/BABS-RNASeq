#!/bin/sh

read_distribution.py \
	-i $bam \
	-r $bed \
	> $metrics_filename

