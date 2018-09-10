#!/bin/sh

java -Xmx10g -Djava.io.tmpdir=$tmp_dirname \
	-jar \$EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
	VALIDATION_STRINGENCY=SILENT \
	INPUT=$bam \
	OUTPUT=$metrics_filename \
	REF_FLAT=$refflat \
	STRAND_SPECIFICITY=NONE \
	RIBOSOMAL_INTERVALS=$rrna_interval_list \
	TMP_DIR=$tmp_dirname

