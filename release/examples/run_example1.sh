#/bin/bash

set -e 

PROG=`basename $0`

echo ${PROG}: This program is part of mGene version 0.1 beta.
echo Please make sure you read the DISCLAIMER.
echo
echo This script trains a single signal and content predictor, predicts and performs an evaluation
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
  echo Note: Running this script takes about 5 minutes \(on a single CPU\).
  FASTA_INPUT=data/elegans.fasta 
  FASTA_DIR=data/ 
  GFF3_INPUT=data/elegans.gff3 
fi
if [ "$1" == 'nGASP' ];
then
  echo Note: Running this script takes about 1h minutes \(on a single CPU\).
  FASTA_INPUT=data/nGASP-Train.fasta 
  GFF3_INPUT=data/nGASP-Train.gff3 
fi

BINDIR=../bin
RESULTDIR=./results-1-$1
mkdir -p $RESULTDIR

echo All results can be found in $RESULTDIR
echo

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Genome and Annotation %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo
 
echo 1a. read DNA fasta file and create a genome information object \(GIO\) \[Log file in ${RESULTDIR}/elegans-genometool.log\]
echo "${BINDIR}/genometool ${FASTA_INPUT} ${RESULTDIR}/elegans.gio "
${BINDIR}/genometool ${FASTA_INPUT} ${FASTA_DIR} ${RESULTDIR}/elegans.gio/log ${RESULTDIR}/elegans.gio #> ${RESULTDIR}/elegans-genometool.log

echo 1b. load the genome annotation in GFF3 format \(format version = wormbase\), create an annotation object \[Log file in ${RESULTDIR}/elegans-gff2anno.log\]
echo "${BINDIR}/gff2anno ${RESULTDIR}/elegans.gio ${GFF3_INPUT} ${RESULTDIR}/elegans.anno wormbase"
${BINDIR}/gff2anno ${RESULTDIR}/elegans.gio ${GFF3_INPUT} ${RESULTDIR}/elegans.anno/logfile ${RESULTDIR}/elegans.anno wormbase "-" "-" "-" || exit -1
#${BINDIR}/gff2anno ${RESULTDIR}/elegans.gio ${GFF3_INPUT} ${RESULTDIR}/elegans.anno/logfile ${RESULTDIR}/elegans.anno wormbase "-" "-" "-" #> ${RESULTDIR}/elegans-gff2anno.log 

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. Signal Prediction %
echo %%%%%%%%%%%%%%%%%%%%%%%%
echo

echo 2a. generate labels usable for training the signal detector for acceptor splice sites \[Log file in ${RESULTDIR}/elegans-anno2signallabel-acc.log\]
${BINDIR}/anno2signallabel ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-acc-label.bspf/ acc > ${RESULTDIR}/elegans-anno2signallabel-acc.log
# Alternative: output in human readable format (much slower)
#${BINDIR}/mgene_anno2signallabel ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-acc-label.spf acc

echo 2b. train the signal detector \[Log file in ${RESULTDIR}/elegans-signal_train-acc.log\]
rm -rf ${RESULTDIR}/elegans-acc.tsp
${BINDIR}/signal_train ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans-acc-label.bspf/ ${RESULTDIR}/elegans-acc.tsp > ${RESULTDIR}/elegans-signal_train-acc.log

echo 2c. use signal detector to predict on genomic DNA \[Log file in ${RESULTDIR}/elegans-signal_predict-acc.log\]
${BINDIR}/signal_predict ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans-acc.tsp ${RESULTDIR}/elegans-acc-prediction.bspf/ > ${RESULTDIR}/elegans-signal_predict-acc.log

echo 2d. perform evaluation
${BINDIR}/signal_eval ${RESULTDIR}/elegans-acc-label.bspf/ ${RESULTDIR}/elegans-acc-prediction.bspf/ | tail -6 


echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Content Prediction %
echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo 3a. generate labels usable for training the content detector for cds_exon \[Log file in ${RESULTDIR}/elegans-anno2contentlabel-cds_exon.log\]
${BINDIR}/anno2contentlabel ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-cds_exon-label.bcpf/ cds_exon > ${RESULTDIR}/elegans-anno2contentlabel-cds_exon.log

echo 3b. train the content detector \[Log file in ${RESULTDIR}/elegans-content_train-cds_exon.log\]
rm -rf ${RESULTDIR}/elegans-cds_exon.tcp/
${BINDIR}/content_train ${RESULTDIR}/elegans.gio/ ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-cds_exon-label.bcpf/ ${RESULTDIR}/elegans-cds_exon.tcp/ > ${RESULTDIR}/elegans-content_train-cds_exon.log

echo 3c. use content detector to predict on genomic DNA \[Log file in ${RESULTDIR}/elegans-content_predict-cds_exon.log\]
${BINDIR}/content_predict ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans-cds_exon.tcp ${RESULTDIR}/elegans-cds_exon-prediction.bspf/ > ${RESULTDIR}/elegans-content_predict-cds_exon.log

echo 3d. perform evaluation
${BINDIR}/content_eval ${RESULTDIR}/elegans-cds_exon-label.bcpf/ ${RESULTDIR}/elegans-cds_exon-prediction.bspf/ | tail -6 

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo


