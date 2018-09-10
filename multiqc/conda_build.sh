#!/bin/sh

module load Anaconda2/5.1.0

# output
BUILD_DIR=./build
ARCHIVE_DIR=./archive

# build
conda build \
	--python 3.6 \
	--croot $BUILD_DIR \
	--output-folder $ARCHIVE_DIR \
	multiqc_plugins

# clean
mv -v $ARCHIVE_DIR/linux-64/multiqc_plugins-1.0-py36_2.tar.bz2 ./
rm -rfv $ARCHIVE_DIR
rm -rfv $BUILD_DIR

