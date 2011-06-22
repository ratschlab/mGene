#!/bin/bash

source `pwd`/mgene_config.sh

if test "$INTERPRETER" = matlab; then
	$MATLAB_BIN_PATH $MATLAB_OPTS -r "paths; try, lpenv = cplex_license(); fprintf('\nqp_test: success\n'); catch, fprintf('\nqp_test: failure\n'); end, exit"  2>/dev/null | grep "qp_test"
else
	$OCTAVE_BIN_PATH --eval "paths; try, lpenv = cplex_license(); fprintf('\nqp_test: success\n'); catch, fprintf('\nqp_test: failure\n'); end, exit"  2>/dev/null | grep "qp_test"
fi
