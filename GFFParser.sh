#!/bin/bash
## USAGE: GFFParser.sh <GFF file> <Result file name> <gio>

source `dirname $0`/../mgene_config.sh

### TODO: either make sure that the right python/scipy version is available in the configure file
### 		or adjust the savemat call (the problem is that the "oned_as" keyword does not 
### 		exist in all scipy versions) to work with different python versions
export PYTHONPATH=/cbio/grlab/home/vipin/lib/python2.5/site-packages/

python `dirname $0`/GFFParser.py $1 $2

if test "$INTERPRETER" = "matlab"; then
	$MATLAB_BIN_PATH $MATLAB_OPTS -r "addpath $MGENE_SRC_PATH; paths; arrange_genes('$2');quit;"
	$MATLAB_BIN_PATH $MATLAB_OPTS -r "addpath $MGENE_SRC_PATH; paths; CurateGenes('$3', '$2'); quit;"
else
	$OCTAVE_BIN_PATH $OCTAVE_OPTS --eval "addpath $MGENE_SRC_PATH; paths; arrange_genes('$2');quit;"
	$OCTAVE_BIN_PATH $OCTAVE_OPTS --eval "addpath $MGENE_SRC_PATH; paths; CurateGenes('$3', '$2'); quit;"
fi
