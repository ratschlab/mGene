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

. `dirname $0`/../src/mgene_config.sh 

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
  ${OCTAVE_BIN_PATH} --eval "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; addpath `dirname $0`; run_matlab_command('$1', '$MGENE_SRC_PATH', $2); global MATLAB_RETURN_FILE; MATLAB_RETURN_FILE=-1 ; exit ;" || (echo starting Octave failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
fi

if [ "$INTERPRETER" == 'matlab' ];
  then
  	${MATLAB_BIN_PATH} $MATLAB_OPTS -r "dbstop error; addpath `dirname $0`; run_matlab_command('$1', '$MGENE_SRC_PATH', $2)" || (echo starting Matlab failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
fi

test -f $MATLAB_RETURN_FILE || exit -1
#ret=`cat $MATLAB_RETURN_FILE` ;
rm -f $MATLAB_RETURN_FILE
exit 0


