#/bin/bash

set -e 

PROG=`basename $0`

echo ${PROG}: This program is part of mGene version 0.2.
echo Please make sure you read the DISCLAIMER.
echo
echo This script trains one signal \(acceptor splice sites\) 
echo and one content predictor \(coding exon sequences\).
echo 
OPT=$1
if [ -z "$OPT" ]; then OPT="--help"; fi 

case "$OPT" in
	"--help")
  		echo Usage: $0 small\|nGASP
  		echo "   or:" $0 --help
  		false
	;;
	"small")
  		echo Note: Running this script takes about 5 minutes \(on a single CPU\).
  		FASTA_INPUT=data/elegans.fasta 
  		FASTA_DIR=data/ 
  		GFF3_INPUT=data/elegans.gff3 
	;;
	"nGASP")
		echo Note: Running this script takes 60 minutes \(on a single CPU\).
		FASTA_INPUT=data/nGASP-Train.fasta 
		FASTA_DIR=data/ 
		GFF3_INPUT=data/nGASP-Train.gff3 
	;;
	*)
  		echo invalid parameter
  		false
	;;
esac



BINDIR=../bin
RESULTDIR=./results-1-$OPT
mkdir -p $RESULTDIR

echo All results can be found in $RESULTDIR
echo

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Genome and Annotation %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

GENOME_DIR=${RESULTDIR}/elegans.gio/
mkdir -p ${GENOME_DIR}
if ! test -f ${GENOME_DIR}/genome.config; then
	echo 1a. read DNA fasta file and create a genome information object \(GIO\) \[Log file in ${GENOME_DIR}/genometool.log\]
	${BINDIR}/genometool ${FASTA_INPUT} ${FASTA_DIR} ${GENOME_DIR}/log ${GENOME_DIR} > ${GENOME_DIR}/genometool.log
fi

ANNO_DIR=${RESULTDIR}/elegans.anno
mkdir -p ${ANNO_DIR}
if ! test -f ${ANNO_DIR}/genes.mat; then
	echo 1b. load the genome annotation in GFF3 format, create an annotation object \[Log file in ${ANNO_DIR}/GFFParser.log\]
	../src/parsegff/GFFParser.sh ${GFF3_INPUT} ${ANNO_DIR}/genes.mat ${GENOME_DIR} > ${ANNO_DIR}/GFFParser.log
fi


echo
echo %%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. Signal Prediction %
echo %%%%%%%%%%%%%%%%%%%%%%%%
echo

sig=acc

SIG_DIR=${RESULTDIR}/elegans-$sig-signal
mkdir -p ${SIG_DIR}/label
mkdir -p ${SIG_DIR}/trained_predictor
mkdir -p ${SIG_DIR}/prediction

echo 2a. generate labels usable for training the signal detector for $sig sites \[Log file in ${SIG_DIR}/anno2signallabel.log\]
${BINDIR}/anno2signallabel ${SIG_DIR}/anno2signallabel.log.tmp ${ANNO_DIR} ${GENOME_DIR} ${SIG_DIR}/label/tmp.txt  ${SIG_DIR}/label/ $sig > ${SIG_DIR}/anno2signallabel.log

echo 2b. train the signal detector \[Log file in ${SIG_DIR}/signal_train.log\]
rm -rf ${SIG_DIR}/trained_predictor
${BINDIR}/signal_train  '-' ${SIG_DIR}/label/ $sig ${GENOME_DIR} ${SIG_DIR}/trained_predictor/log ${SIG_DIR}/trained_predictor 1 > ${SIG_DIR}/signal_train.log

echo 2c. use signal detector to predict on genomic DNA \[Log file in ${SIG_DIR}/signal_predict.log\]
${BINDIR}/signal_predict ${SIG_DIR}/trained_predictor ${GENOME_DIR} ${SIG_DIR}/prediction/log ${SIG_DIR}/prediction/ 1 > ${SIG_DIR}/signal_predict.log



echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Content Prediction %
echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo

content=cds_exon
CONT_DIR=${RESULTDIR}/elegans-$content-content
mkdir -p ${CONT_DIR}/label
mkdir -p ${CONT_DIR}/trained_predictor
mkdir -p ${CONT_DIR}/prediction

echo 3a. generate labels usable for training the content detector for $content \[Log file in ${CONT_DIR}/anno2contentlabel.log\]
${BINDIR}/anno2contentlabel ${ANNO_DIR}/log  ${ANNO_DIR} ${GENOME_DIR} ${CONT_DIR}/label/log ${CONT_DIR}/label $content all_other default #> ${CONT_DIR}/anno2contentlabel.log

echo 3b. train the content detector \[Log file in ${CONT_DIR}/content_train.log\]
rm -rf ${RESULTDIR}/elegans-$content.tcp/
${BINDIR}/content_train ${CONT_DIR}/label/log ${CONT_DIR}/label $content ${GENOME_DIR} ${ANNO_DIR} ${CONT_DIR}/trained_predictor/log ${CONT_DIR}/trained_predictor 1 #> ${CONT_DIR}/content_train.log

echo 3c. use content detector to predict on genomic DNA \[Log file in ${CONT_DIR}/content_predict.log\]
${BINDIR}/content_predict ${CONT_DIR}/trained_predictor ${GENOME_DIR} ${CONT_DIR}/prediction/log ${CONT_DIR}/prediction 1  #> ${CONT_DIR}/content_predict.log

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo

