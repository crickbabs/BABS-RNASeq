#!/bin/sh

java -Xmx10g -Djava.io.tmpdir=$tmp_dirname \
	-jar $EBROOTPICARD/picard.jar EstimateLibraryComplexity \
	VALIDATION_STRINGENCY=SILENT \
	INPUT=$bam \
	OUTPUT=$metrics_filename \
	TMP_DIR=$tmp_dirname
