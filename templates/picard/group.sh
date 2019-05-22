#!/bin/sh

java -Xmx10g -Djava.io.tmpdir=$tmp_dirname \
	-jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
	VALIDATION_STRINGENCY=SILENT \
	INPUT=$bam \
	OUTPUT=$filename \
	RGID=$name \
	RGLB=$name \
	RGPU=$name \
	RGSM=$name \
	RGCN=TheFrancisCrickInsitute \
	RGPL=Illumina

