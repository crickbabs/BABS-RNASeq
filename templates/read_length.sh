#!/bin/sh

zcat $fastq \
	| head -n 16 \
	| sed -n "2~4p" \
	| awk '{ print length }' \
	| sort \
	| uniq \
	| sed -n "1p" \
	| tr -d "\n"

