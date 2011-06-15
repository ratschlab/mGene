#/bin/bash

set -e 

PROG=`basename $0`

echo ${PROG}: This program is part of mGene version 0.1 beta.
echo Please make sure you read the DISCLAIMER.
echo
echo This script trains all six signal, five content predictors, and the gene predictor.
echo Moreover, it performs gene prediction on the DNA and performs an evaluation
echo \(This version uses the individual modules to achieve this goal\)
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
RESULTDIR=./results-2-$1
mkdir -p $RESULTDIR

echo All results can be found in $RESULTDIR
echo

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Genome and Annotation %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo
 
echo 1a. read DNA fasta file and create a genome information object \(GIO\) \[Log file in ${RESULTDIR}/elegans-genometool.log\]
${BINDIR}/mgene_genometool ${FASTA_INPUT} ${RESULTDIR}/elegans.gio > ${RESULTDIR}/elegans-genometool.log

echo 1b. load the genome annotation in GFF3 format \(format version = wormbase\), create an annotation object \[Log file in ${RESULTDIR}/elegans-gff2anno.log\]
${BINDIR}/mgene_gff2anno ${RESULTDIR}/elegans.gio ${GFF3_INPUT} ${RESULTDIR}/elegans.anno wormbase > ${RESULTDIR}/elegans-gff2anno.log

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. Signal Prediction %
echo %%%%%%%%%%%%%%%%%%%%%%%%
echo

for signal in tss tis acc don cdsStop cleave 
do
  echo 2a. generate labels usable for training the signal detector for signal ${signal} \[Log file in ${RESULTDIR}/elegans-anno2signallabel-${signal}.log\]
  ${BINDIR}/mgene_anno2signallabel ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-${signal}-label.bspf/ ${signal} > ${RESULTDIR}/elegans-anno2signallabel-${signal}.log

  echo 2b. train the signal detector \[Log file in ${RESULTDIR}/elegans-signal_train-${signal}.log\]
  rm -rf ${RESULTDIR}/elegans-${signal}.tsp
  ${BINDIR}/mgene_signal_train ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans-${signal}-label.bspf/ ${RESULTDIR}/elegans-${signal}.tsp > ${RESULTDIR}/elegans-signal_train-${signal}.log
  
  echo 2c. use signal detector to predict on genomic DNA \[Log file in ${RESULTDIR}/elegans-signal_predict-${signal}.log\]
  ${BINDIR}/mgene_signal_predict ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans-${signal}.tsp ${RESULTDIR}/elegans-${signal}-prediction.bspf/ > ${RESULTDIR}/elegans-signal_predict-${signal}.log

  echo 2d. perform evaluation
  ${BINDIR}/mgene_signal_eval ${RESULTDIR}/elegans-${signal}-label.bspf/ ${RESULTDIR}/elegans-${signal}-prediction.bspf/ | tail -6 

  echo
done

echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Content Prediction %
echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo

for content in intergenic utr5exon cds_exon utr3exon intron
do
  echo 3a. generate labels usable for training the content detector for ${content} \[Log file in ${RESULTDIR}/elegans-anno2contentlabel-${content}.log\]
  ${BINDIR}/mgene_anno2contentlabel ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-${content}-label.bcpf/ ${content} > ${RESULTDIR}/elegans-anno2contentlabel-${content}.log
  
  echo 3b. train the content detector \[Log file in ${RESULTDIR}/elegans-content_train-${content}.log\]
  rm -rf ${RESULTDIR}/elegans-${content}.tcp/
  ${BINDIR}/mgene_content_train ${RESULTDIR}/elegans.gio/ ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-${content}-label.bcpf/ ${RESULTDIR}/elegans-${content}.tcp/ > ${RESULTDIR}/elegans-content_train-${content}.log
  
  echo 3c. use content detector to predict on genomic DNA \[Log file in ${RESULTDIR}/elegans-content_predict-${content}.log\]
  ${BINDIR}/mgene_content_predict ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans-${content}.tcp ${RESULTDIR}/elegans-${content}-prediction.bspf/ > ${RESULTDIR}/elegans-content_predict-${content}.log

  echo 3d. perform evaluation
  ${BINDIR}/mgene_content_eval ${RESULTDIR}/elegans-${content}-label.bcpf/ ${RESULTDIR}/elegans-${content}-prediction.bspf/ | tail -6 

  echo
done

echo %%%%%%%%%%%%%%%%%%%
echo % Gene Prediction %
echo %%%%%%%%%%%%%%%%%%%

echo 4a. training the gene predictor \[Log file in ${RESULTDIR}/elegans-gene_train.log\]
${BINDIR}/mgene_gene_train ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.anno \
	${RESULTDIR}/elegans-tss-prediction.bspf ${RESULTDIR}/elegans-tis-prediction.bspf ${RESULTDIR}/elegans-acc-prediction.bspf \
	${RESULTDIR}/elegans-don-prediction.bspf ${RESULTDIR}/elegans-cdsStop-prediction.bspf ${RESULTDIR}/elegans-cleave-prediction.bspf \
	${RESULTDIR}/elegans-intergenic-prediction.bspf ${RESULTDIR}/elegans-utr5exon-prediction.bspf ${RESULTDIR}/elegans-cds_exon-prediction.bspf \
	${RESULTDIR}/elegans-utr3exon-prediction.bspf ${RESULTDIR}/elegans-intron-prediction.bspf ${RESULTDIR}/elegans.tgp > ${RESULTDIR}/elegans-gene_train.log

echo 4b. using the trained gene predictor for prediction \[Log file in ${RESULTDIR}/elegans-gene_predict.log\]
${BINDIR}/mgene_gene_predict ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.tgp \
	${RESULTDIR}/elegans-tss-prediction.bspf ${RESULTDIR}/elegans-tis-prediction.bspf ${RESULTDIR}/elegans-acc-prediction.bspf \
	${RESULTDIR}/elegans-don-prediction.bspf ${RESULTDIR}/elegans-cdsStop-prediction.bspf ${RESULTDIR}/elegans-cleave-prediction.bspf \
	${RESULTDIR}/elegans-intergenic-prediction.bspf ${RESULTDIR}/elegans-utr5exon-prediction.bspf ${RESULTDIR}/elegans-cds_exon-prediction.bspf \
	${RESULTDIR}/elegans-utr3exon-prediction.bspf ${RESULTDIR}/elegans-intron-prediction.bspf ${RESULTDIR}/elegans-prediction.gff3 > ${RESULTDIR}/elegans-gene_predict.log

echo 4c. comparing the reference annotation and the gene prediction \[Log file in ${RESULTDIR}/elegans-gene_eval.log\]
${BINDIR}/mgene_gene_eval ${RESULTDIR}/elegans.gio ${RESULTDIR}/elegans.anno ${RESULTDIR}/elegans-prediction.gff3 > ${RESULTDIR}/elegans-gene_eval.log
tail -30 ${RESULTDIR}/elegans-gene_eval.log

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo

