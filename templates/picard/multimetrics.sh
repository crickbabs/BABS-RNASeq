#!/bin/sh

java -Xmx10g -Djava.io.tmpdir=$tmp_dirname \
	-jar \$EBROOTPICARD/picard.jar CollectMultipleMetrics \
	VALIDATION_STRINGENCY=SILENT \
	INPUT=$bam \
	OUTPUT=$metrics_filename \
	PROGRAM=CollectAlignmentSummaryMetrics \
	R=$fasta \
	TMP_DIR=$tmp_dirname

