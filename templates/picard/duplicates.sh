#!/bin/sh

java -Xmx10g -Djava.io.tmpdir=$tmp_dirname \
	-jar \$EBROOTPICARD/picard.jar MarkDuplicates \
	VALIDATION_STRINGENCY=SILENT \
	INPUT=$bam \
	OUTPUT=$filename \
	METRICS_FILE=$metrics_filename \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=false \
	TMP_DIR=$tmp_dirname

