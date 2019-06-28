#!/bin/sh

samtools sort \
	--threads $task.cpus \
	-o $filename \
	$bam

samtools index $filename

samtools stats $bam  | grep "insert size " | sed --expression "s/SN/$name/" | awk -v FS='\t' -v OFS='\t' 'NR % 2 == 1 { o=\$0 ; next } { print o, \$3 }' > insert_size.txt
