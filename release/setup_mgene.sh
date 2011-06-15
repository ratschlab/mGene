#/bin/bash

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

. ./bin/mgene_config.sh

echo mGene setup script \(version 0.0.1\) 
echo ==================================
echo

echo mGene base directory \(currently set to \"$MGENE_PATH\", suggest to set to \"`pwd`\"\)
read MGENE_PATH       
echo '=>' Setting mGene base directory to \"$MGENE_PATH\"
echo

echo mGene source code directory \(currently set to \"$MGENE_SRC_PATH\", suggest to set to \"`pwd`/src\"\)
read MGENE_SRC_PATH       
echo '=>' Setting mGene base directory to \"$MGENE_SRC_PATH\"
echo

echo Shogun binary directory \(currently set to \"$SHOGUN_PATH\", suggest to set to \"`pwd`/shogun\"\)
read SHOGUN_PATH
echo '=>' Setting Shogun binary directory to \"$SHOGUN_PATH\"
echo

echo Enable support for the Sun Gridengine '(y/n)'? [This option is still experimental and is not fully supported yet.]
read USE_SGE
echo

if [ "$USE_SGE" != 'y'  -a "$USE_SGE" != 'n' ];
then
	echo Unrecognized choice: \"$USE_SGE\"
	echo Aborting
	false
fi

SGE_BIN_PATH=

if [ "$USE_SGE" != 'n' ];
then
	echo Sun Gridengine directory \(currently set to \"$SGE_BIN_PATH\", suggest to set to \"/usr/local/sge\"\)
	read SGE_BIN_PATH
	echo '=>' Setting Sun gridengine directory to \"$SGE_BIN_PATH\"
	echo
fi


echo Which interpreter should be used \(\"octave\" or \"matlab\"\)
read INTERPRETER  

if [ "$INTERPRETER" != 'octave' -a  "$INTERPRETER" != 'matlab' ];
then
	echo Unrecognized choice: \"$INTERPRETER\"
	echo Aborting
	false
fi
echo '=>' Setting interpreter to \"$INTERPRETER\"
echo

if [ "$INTERPRETER" == 'octave' ];
then
	echo Please enter the full path to octave \(currently set to \"$OCTAVE_BIN_PATH\"\)
	read OCTAVE_BIN_PATH
	echo '=>' Setting octave\'s path to \"$OCTAVE_BIN_PATH\"
	echo

	MATLAB_BIN_PATH=

	echo Which optimizer should be used \(octave\'s builtin \"glpk\" or \"mosek\", currently set to \"$OPTIMIZER\"\)?
	read OPTIMIZER
	if [ "$OPTIMIZER" != 'glpk' -a "$OPTIMIZER" != 'mosek' ];
		then
		echo Unrecognized choice: \"$OPTIMIZER\"
		echo Aborting
		false;
	fi
	echo '=>' Setting the optimizer to \"$OPTIMIZER\"
	echo 
	
	OPTIMIZER_PATH=
	if [ "$OPTIMIZER" == 'mosek' ];
		then
		echo Please enter the path to the Mosek support files \(suggest to set to \"$MGENE_PATH/mosek\"\)
		read OPTIMIZER_PATH
		export OPTIMIZER_TOOLBOX_PATH=${OPTIMIZER_PATH}/octave
		echo '=>' Setting the optimizer path to \"$OPTIMIZER_PATH\";
		echo '=>' Setting the optimizer toolbox path to \"$OPTIMIZER_TOOLBOX_PATH\";
	fi
fi

if [ "$INTERPRETER" == 'matlab' ];
then
	echo Please enter the full path to matlab \(currently set to \"$MATLAB_BIN_PATH\"\)
	read MATLAB_BIN_PATH
	echo '=>' Setting octave\'s path to \"$MATLAB_BIN_PATH\"
	echo 

	OCTAVE_BIN_PATH=

	echo Which optimizer should be used \(matlabs builtin \"quadprog\" command or \"mosek\", currently set to \"$OPTIMIZER\"\)?
	read OPTIMIZER
	if [ "$OPTIMIZER" != 'quadprog' -a "$OPTIMIZER" != 'mosek' ];
		then
		echo Unrecognized choice: \"$OPTIMIZER\"
		echo Aborting
		false
	fi
	echo '=>' Setting the optimizer to \"$OPTIMIZER\"
	echo
	
	OPTIMIZER_PATH=
	if [ "$OPTIMIZER" == 'mosek' ];
		then
		echo Please enter the path to the Mosek support files \(suggest to set to \"$MGENE_PATH/mosek\"\)
		read OPTIMIZER_PATH
		echo '=>' Setting the optimizer path to \"$OPTIMIZER_PATH\"
		export OPTIMIZER_TOOLBOX_PATH=${OPTIMIZER_PATH}/matlab
		echo '=>' Setting the optimizer toolbox path to \"$OPTIMIZER_TOOLBOX_PATH\";
	fi
fi

cp bin/mgene_config.sh bin/mgene_config.sh.bak
grep -v -e OCTAVE_BIN_PATH -e OCTAVE_BIN_PATH -e MATLAB_BIN_PATH -e MGENE_PATH -e MGENE_SRC_PATH \
	-e SHOGUN_PATH -e SGE_BIN_PATH -e INTERPRETER bin/mgene_config.sh.bak \
	-e OPTIMIZER_PATH -e OPTIMIZER_TOOLBOX_PATH -e OPTIMIZER > bin/mgene_config.sh

# appending the relevant lines to mgene_config.sh
echo export MGENE_PATH=$MGENE_PATH >> bin/mgene_config.sh
echo export MGENE_SRC_PATH=$MGENE_SRC_PATH >> bin/mgene_config.sh
echo export SHOGUN_PATH=$SHOGUN_PATH >> bin/mgene_config.sh
echo export SGE_BIN_PATH=$SGE_BIN_PATH >> bin/mgene_config.sh
echo export INTERPRETER=$INTERPRETER >> bin/mgene_config.sh
echo export MATLAB_BIN_PATH=$MATLAB_BIN_PATH >> bin/mgene_config.sh
echo export OCTAVE_BIN_PATH=$OCTAVE_BIN_PATH >> bin/mgene_config.sh
echo export OPTIMIZER=$OPTIMIZER >> bin/mgene_config.sh
echo export OPTIMIZER_PATH=$OPTIMIZER_PATH >> bin/mgene_config.sh
echo export OPTIMIZER_TOOLBOX_PATH=$OPTIMIZER_TOOLBOX_PATH >> bin/mgene_config.sh

echo
echo Done.
echo 
