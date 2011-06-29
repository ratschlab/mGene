#!/bin/bash

CURRENT=`pwd`

#rm -rf release/src
mkdir -p release/src
mkdir -p release/src/auxiliary_data/
mkdir -p release/src/shared_tools/
mkdir -p release/src/parsegff/

sh copy_files.sh

# remove binaries
find release -name \*.mexa64 -exec rm -rf '{}' ';'
find release -name \*.mex -exec rm -rf '{}' ';'
find release -name \*.oct -exec rm -rf '{}' ';'
find release -name \*.o -exec rm -rf '{}' ';'
rm src/samtools/samtools

RDIR=mgene-0.2.0 ;
mkdir -p $RDIR
rsync -a release/ $RDIR
find $RDIR -name .*.swp -exec rm -rf '{}' ';'
find $RDIR -name .svn -exec rm -rf '{}' ';'
find $RDIR -name cachegrind.* -exec rm -rf '{}' ';'
find $RDIR -name result.* -exec rm -rf '{}' ';'

