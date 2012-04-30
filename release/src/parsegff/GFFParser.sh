#!/bin/bash
## USAGE: GFFParser.sh <GFF file> <Result file name> <gio>

source `dirname $0`/../mgene_config.sh

### TODO: either make sure that the right python/scipy version is available in the configure file
### 		or adjust the savemat call (the problem is that the "oned_as" keyword does not 
### 		exist in all scipy versions) to work with different python versions
#export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/vipin/lib/python2.5/site-packages/
export PYTHONPATH=/home/galaxy/lib/lib/python2.6/dist-packages/:/home/galaxy/lib/python2.6/site-packages/:/home/galaxy/lib/python2.5/site-packages/:/home/galaxy/svn/projects/genefinding/parsegff/:/home/galaxy/python/lib/python2.5/site-packages/:/fml/ag-raetsch/home/jonas/shogun/python/lib/python2.5/site-packages:/fml/ag-raetsch/home/jonas/shogun/trunk/src:/fml/ag-raetsch/home/jonas/shogun/trunk/src/python_modular

if ! test -f $2
then
	python `dirname $0`/GFFParser.py $1 $2
fi

mkdir -p `dirname $2`

echo python_call: python `dirname $0`/GFFParser.py $1 $2
python `dirname $0`/GFFParser.py $1 $2 2>&1

if test "$INTERPRETER" = "matlab"; then
	echo MGENE_SRC_PATH:$MGENE_SRC_PATH
	$MATLAB_BIN_PATH $MATLAB_OPTS -r "addpath $MGENE_SRC_PATH; dbstop error; paths; arrange_genes('$2');quit;"
	$MATLAB_BIN_PATH $MATLAB_OPTS -r "addpath $MGENE_SRC_PATH; paths; CurateGenes('$3', '$2'); quit;"
else
	echo MGENE_SRC_PATH:$MGENE_SRC_PATH
	$OCTAVE_BIN_PATH $OCTAVE_OPTS --eval "addpath $MGENE_SRC_PATH; paths; arrange_genes('$2');quit;"
	$OCTAVE_BIN_PATH $OCTAVE_OPTS --eval "addpath $MGENE_SRC_PATH; paths; CurateGenes('$3', '$2'); quit;"
fi
