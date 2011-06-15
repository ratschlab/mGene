#!/bin/bash

rm -rf release/src
mkdir -p release/src
mkdir -p release/src/auxiliary_data/
mkdir -p release/src/shared_tools/

sh copy_files.sh

export RDIR=mgene-0.2.0 ;
mkdir -p $RDIR
rsync -a release/ $RDIR

