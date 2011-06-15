#!/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2005-2009 Gunnar Raetsch, Gabriele Schweikert, Jonas Behr, Soeren Sonnenburg, Alexander Zien, Georg Zeller, Andre Noll, Cheng Soon Ong, Petra Philips
# Copyright (C) 2005-2009 Max Planck Society and Fraunhofer Institute FIRST 
#

set -e

. `dirname $0`/mgene_config.sh 

export LD_LIBRARY_PATH=${SHOGUN_PATH}:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${SHOGUN_PATH}:${LD_LIBRARY_PATH}

export PATH=${SGE_BIN_PATH}:${PATH}

if [ "$OPTIMIZER" == 'mosek' ];
  then
	export LD_LIBRARY_PATH=${OPTIMIZER_PATH}:${LD_LIBRARY_PATH}
    export DYLD_LIBRARY_PATH=${OPTIMIZER_PATH}:${DYLD_LIBRARY_PATH}

	if [ -f ${OPTIMIZER_PATH}/mosek.lic ];
    then 
		export MOSEKLM_LICENSE_FILE=${OPTIMIZER_PATH}/mosek.lic
    fi 
fi

export MGENE_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}

export MATLAB_RETURN_FILE=`tempfile`

if [ "$INTERPRETER" == 'octave' ];
  then
  #${OCTAVE_BIN_PATH} --eval "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; addpath $MGENE_SRC_PATH; mgene_config; try, $1($2); catch; fprintf(2, '$1 failed\n') ; global MATLAB_RETURN_FILE; MATLAB_RETURN_FILE=-1 ; end ; exit ;" || (echo starting Octave failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
  ${OCTAVE_BIN_PATH} --eval "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; addpath $MGENE_SRC_PATH; mgene_config; $1($2); exit ;"
fi

if [ "$INTERPRETER" == 'matlab' ];
  then
  ${MATLAB_BIN_PATH} -r "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; addpath $MGENE_SRC_PATH; mgene_config; try, $1($2); catch; fprintf(2, '$1 failed\n') ; global MATLAB_RETURN_FILE; MATLAB_RETURN_FILE=-1 ; end ; exit ;" || (echo starting Matlab failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
  #${MATLAB_BIN_PATH} -r "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; addpath $MGENE_SRC_PATH; mgene_config; dbstop error; $1($2) ; exit ;" || (echo starting Matlab failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
fi

test -f $MATLAB_RETURN_FILE || exit 0
ret=`cat $MATLAB_RETURN_FILE` ;
rm -f $MATLAB_RETURN_FILE
exit $ret


