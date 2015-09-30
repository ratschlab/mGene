#!/bin/bash
#export PYTHONPATH=/cbio/grlab/home/jonas/python/usr/local/lib/python2.5/site-packages/
export PYTHONPATH=/cbio/grlab/home/jonas/python/usr/local/lib/python2.6/dist-packages/
#export PYTHONPATH=/cbio/grlab/home/cwidmer/lib/python2.5/site-packages/
#export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/jonas/shogun/trunk/src:/cbio/grlab/home/jonas/shogun/trunk/src/python_modular
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/cwidmer/svn/projects/multitask/python/base/
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/cwidmer/svn/projects/arts2/
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/jonas/svn/tools/python/
#export LD_LIBRARY_PATH=/cbio/grlab/home/jonas/shogun/python/lib/
#export LD_LIBRARY_PATH=/cbio/grlab/home/cwidmer/lib/ ## chris shogun fuer pygene
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cbio/grlab/home/cwidmer/lib/ ## chris shogun fuer pygene
export LD_LIBRARY_PATH=/cbio/grlab/home/jonas/shogun/trunk/src/libshogun:/cbio/grlab/home/jonas/shogun/trunk/src/libshogunui

python pygene.py

