#!/bin/bash

source `pwd`/mgene_config.sh

if test "$INTERPRETER" = matlab; then
	$MATLAB_BIN_PATH $MATLAB_OPTS -r "paths; try, lpenv = cplex_license(); fprintf('\nqp_test_succ\n'); catch, fprintf('\nqp_test_failed\n'); end, exit"  2>/dev/null | grep "qp_test"
else
	echo qp_test_failed
fi
