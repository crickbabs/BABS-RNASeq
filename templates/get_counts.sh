#!/bin/sh

cat $filename \
	| cut -d\$'\\t' -f 1,5 \
	| sed 's/expected_count/$name/' \
	> $counts

