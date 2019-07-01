#!/bin/sh

samtools sort \
	--threads $task.cpus \
	-o $filename \
	$bam

samtools index $filename
