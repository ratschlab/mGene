#!/bin/bash

source `pwd`/mgene_config.sh

if test "$INTERPRETER" = matlab; then
	LD_LIBRARY_PATH=$SHOGUN_PATH $MATLAB_BIN_PATH $MATLAB_OPTS -r "addpath $SHOGUN_PATH; try, x = sg('get_version'); fprintf('\nsg_test: success\n'); catch, fprintf('\nsg_test: failure\n'); end, exit" | grep "sg_test"
else
	LD_LIBRARY_PATH=$SHOGUN_PATH echo "addpath $SHOGUN_PATH; try, x = sg('get_version'); fprintf('\nsg_test: success\n'); catch, fprintf('\nsg_test: failure\n'); end, exit" |$OCTAVE_BIN_PATH $OCTAVE_OPTS | grep "sg_test"
fi
