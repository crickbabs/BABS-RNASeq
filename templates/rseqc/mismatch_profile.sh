#!/bin/sh

mismatch_profile.py \
	-i $bam \
	-l $rough_read_length \
	-o $metrics_filename

