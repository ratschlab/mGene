#!/bin/bash

#rm -rf release/src
mkdir -p release/src
mkdir -p release/src/auxiliary_data/
mkdir -p release/src/shared_tools/
mkdir -p release/src/parsegff/

#sh copy_files.sh

# remove binaries
#find release -name \*.mexa64 -exec rm -rf '{}' ';'
#find release -name \*.mex -exec rm -rf '{}' ';'
#find release -name \*.oct -exec rm -rf '{}' ';'
#find release -name \*.o -exec rm -rf '{}' ';'
#rm src/samtools/samtools

RDIR=mgene-0.2.0-beta ;
mkdir -p $RDIR
rsync -a release/ $RDIR
find $RDIR -name .*.swp -exec rm -rf '{}' ';'
find $RDIR -name .svn -exec rm -rf '{}' ';'
find $RDIR -name cachegrind.* -exec rm -rf '{}' ';'
find $RDIR -name result.* -exec rm -rf '{}' ';'
find $RDIR -name \*.mexa64 -exec rm -rf '{}' ';'
find $RDIR -name \*.mex -exec rm -rf '{}' ';'
find $RDIR -name \*.oct -exec rm -rf '{}' ';'
find $RDIR -name \*.o -exec rm -rf '{}' ';'
rm $RDIR/src/samtools/samtools
rm -rf $RDIR/examples/results*
rm $RDIR/src/shared_tools/ngs/*.txt
rm $RDIR/src/shared_tools/ngs/*.pickle

tar -czvf $RDIR.tar.gz $RDIR

