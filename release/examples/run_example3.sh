#/bin/bash

set -e 

PROG=`basename $0`

echo ${PROG}: This program is part of mGene version 0.1 beta.
echo Please make sure you read the DISCLAIMER.
echo
echo This script trains all six signal, five content predictors, and the gene predictor.
echo Moreover, it performs gene prediction on the DNA and performs an evaluation
echo \(This version uses the monolithic workflows to achieve this goal\)
echo

if [ -z "$1" -o "$1" == '--help' ];
then
  echo Usage: $0 small\|nGASP
  echo "   or:" $0 --help
  false
fi 
if [ "$1" != 'small' -a "$1" != 'nGASP' ];
then
  echo invalid parameter
  false
fi

if [ "$1" == 'small' ];
then
  echo Note: Running this script takes about 60 minutes \(on a single CPU\).
  FASTA_INPUT=data/elegans.fasta 
  GFF3_INPUT=data/elegans.gff3 
fi
if [ "$1" == 'nGASP' ];
then
  echo Note: Running this script takes 10-20 hours \(on a single CPU\).
  FASTA_INPUT=data/nGASP-Train.fasta 
  GFF3_INPUT=data/nGASP-Train.gff3 
fi

BINDIR=../bin
RESULTDIR=./results-3-$1
mkdir -p $RESULTDIR

echo All results can be found in $RESULTDIR
echo

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 1. mGene Train Workflow %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo
 
echo trains classifiers for six signals, 5 contents and the predictor for gene prediction \[Log file in ${RESULTDIR}/elegans-train_workflow.log\]
rm -rf ${RESULTDIR}/elegans.tmgp
${BINDIR}/mgene_train_workflow ${FASTA_INPUT} ${GFF3_INPUT} wormbase ${RESULTDIR}/elegans.tmgp > ${RESULTDIR}/elegans-train_workflow.log

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. mGene Predict Workflow %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo uses the trained mGene predictor to predict genes on piece of DNA \[Log file in ${RESULTDIR}/elegans-predict_workflow.log\]
${BINDIR}/mgene_predict_workflow ${FASTA_INPUT} ${RESULTDIR}/elegans.tmgp ${RESULTDIR}/elegans-prediction.gff3 > ${RESULTDIR}/elegans-predict_workflow.log

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. mGene Evaluation Workflow %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo reads the reference genome annotation and compares it against the gene prediction
${BINDIR}/mgene_eval_workflow ${FASTA_INPUT} ${GFF3_INPUT} wormbase ${RESULTDIR}/elegans-prediction.gff3 0

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
