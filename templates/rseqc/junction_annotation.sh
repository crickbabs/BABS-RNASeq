#!/bin/sh

junction_annotation.py \
	-i $bam \
	-r $bed \
	-o $metrics_filename

