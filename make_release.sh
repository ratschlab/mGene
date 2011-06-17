#!/bin/bash

CURRENT=`pwd`

#rm -rf release/src
mkdir -p release/src
mkdir -p release/src/auxiliary_data/
mkdir -p release/src/shared_tools/

sh copy_files.sh

# remove binaries
cd release
	find . -name \*.mexa64 -exec rm -rf '{}' ';'
	find . -name \*.mex -exec rm -rf '{}' ';'
	find . -name \*.oct -exec rm -rf '{}' ';'
	find . -name \*.o -exec rm -rf '{}' ';'
	rm src/samtools/samtools
cd $CURRENT

#RDIR=mgene-0.2.0 ;
#mkdir -p $RDIR
#rsync -a release/ $RDIR
#cd $RDIR && find . -name .svn -exec rm -rf '{}' ';'

