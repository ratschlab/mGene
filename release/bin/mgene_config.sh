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

export MGENE_PATH=/fml/ag-raetsch/home/raetsch/svn/releases/mGeneToolbox-0.1.0/release
export MGENE_SRC_PATH=/fml/ag-raetsch/home/raetsch/svn/releases/mGeneToolbox-0.1.0/release/src
export SHOGUN_PATH=/fml/ag-raetsch/home/raetsch/svn/releases/mGeneToolbox-0.1.0/release/shogun
export SGE_BIN_PATH=/usr/local/sge
export INTERPRETER=matlab
export MATLAB_BIN_PATH=/fml/ag-raetsch/share/software/matlab-7.6/bin/matlab
export OCTAVE_BIN_PATH=/fml/ag-raetsch/home/raetsch/bin/octave
export OPTIMIZER=mosek
export OPTIMIZER_PATH=/fml/ag-raetsch/home/raetsch/svn/releases/mGeneToolbox-0.1.0/release/mosek
export OPTIMIZER_TOOLBOX_PATH=/fml/ag-raetsch/home/raetsch/svn/releases/mGeneToolbox-0.1.0/release/mosek/matlab
