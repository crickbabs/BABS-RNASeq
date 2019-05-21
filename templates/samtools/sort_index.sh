#!/bin/sh

samtools sort \
	--threads $task.cpus \
	-o $filename \
	$bam

samtools index $filename

samtools stats $bam  | grep "insert size " | sed --expression "s/SN/$name/" > insert_size.txt
