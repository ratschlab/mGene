#!/bin/bash
export PYTHONPATH=/cbio/grlab/home/jonas/python/usr/local/lib/python2.6/dist-packages/
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/cwidmer/svn/projects/multitask/python/base/
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/cwidmer/svn/projects/arts2/
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/jonas/svn/tools/python/
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/jonas/shogun/trunk/applications/asp/
export PYTHONPATH=$PYTHONPATH:/cbio/grlab/home/jonas/svn/tools/ngs/
export LD_LIBRARY_PATH=/cbio/grlab/home/jonas/shogun/trunk/src/libshogun:/cbio/grlab/home/jonas/shogun/trunk/src/libshogunui


fafname=/cbio/grlab/nobackup/projects/sequencing_runs/worm-Dec09/genomes/japonica2/japonica2_condensed.fasta

~/shogun/trunk/applications/asp/asp --svm_type primal $fafname 
