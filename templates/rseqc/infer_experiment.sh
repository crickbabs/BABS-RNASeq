#!/bin/sh

infer_experiment.py \
	-i $bam \
	-r $bed \
	> $metrics_filename

