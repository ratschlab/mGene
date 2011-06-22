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
OPT=$1
if [ -z "$OPT" ]; then OPT="--help"; fi 

case "$OPT" in
	"--help")
  		echo Usage: $0 small\|nGASP
  		echo "   or:" $0 --help
  		false
	;;
	"small")
  		echo Note: Running this script takes about 60 minutes \(on a single CPU\).
  		FASTA_INPUT=data/elegans.fasta 
  		FASTA_DIR=data/ 
  		GFF3_INPUT=data/elegans.gff3 
	;;
	"nGASP")
		echo Note: Running this script takes 10-20 hours \(on a single CPU\).
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
RESULTDIR=./results-2-$OPT
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

for sig in tss tis acc don cdsStop cleave 
do
	SIG_DIR=${RESULTDIR}/elegans-$sig-signal
	FN_DONE=${SIG_DIR}/prediction/DONE 
	if test -f $FN_DONE; then echo $sig already done; continue; fi
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

	touch $FN_DONE
  echo
done

echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Content Prediction %
echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo

for content in intergenic utr5exon cds_exon utr3exon intron
do
	CONT_DIR=${RESULTDIR}/elegans-$content-content
	FN_DONE=${CONT_DIR}/prediction/DONE 
	if test -f $FN_DONE; then echo $content already done; continue; fi
	mkdir -p ${CONT_DIR}/label
	mkdir -p ${CONT_DIR}/trained_predictor
	mkdir -p ${CONT_DIR}/prediction

	echo 3a. generate labels usable for training the content detector for $content \[Log file in ${CONT_DIR}/anno2contentlabel.log\]
	${BINDIR}/anno2contentlabel ${ANNO_DIR}/log  ${ANNO_DIR} ${GENOME_DIR} ${CONT_DIR}/label/log ${CONT_DIR}/label $content all_other default > ${CONT_DIR}/anno2contentlabel.log
	
	echo 3b. train the content detector \[Log file in ${CONT_DIR}/content_train.log\]
	rm -rf ${RESULTDIR}/elegans-$content.tcp/
	${BINDIR}/content_train ${CONT_DIR}/label/log ${CONT_DIR}/label $content ${GENOME_DIR} ${ANNO_DIR} ${CONT_DIR}/trained_predictor/log ${CONT_DIR}/trained_predictor 1 > ${CONT_DIR}/content_train.log
	
	echo 3c. use content detector to predict on genomic DNA \[Log file in ${CONT_DIR}/content_predict.log\]
	${BINDIR}/content_predict ${CONT_DIR}/trained_predictor ${GENOME_DIR} ${CONT_DIR}/prediction/log ${CONT_DIR}/prediction 1  > ${CONT_DIR}/content_predict.log

	touch $FN_DONE
  echo
done

echo %%%%%%%%%%%%%%%%%%%
echo % Gene Prediction %
echo %%%%%%%%%%%%%%%%%%%

LSL_DIR=${RESULTDIR}/elegans-lsl
GENE_TRAIN_OPTS="use_rna_seq_for_label_gen=-1"
FN_DONE=$LSL_DIR/TRAIN_DONE
if test -f $FN_DONE; then 
	echo gene_train already done
else
	mkdir -p $LSL_DIR
	echo 4a. training the gene predictor \[Log file in ${LSL_DIR}/gene_train.log\]
	${BINDIR}/gene_train ${GENOME_DIR} ${ANNO_DIR} ${LSL_DIR}/gene_train.log2  \
		${RESULTDIR}/elegans-tss-signal/prediction/log ${RESULTDIR}/elegans-tss-signal/prediction \
		${RESULTDIR}/elegans-tis-signal/prediction/log ${RESULTDIR}/elegans-tis-signal/prediction \
		${RESULTDIR}/elegans-acc-signal/prediction/log ${RESULTDIR}/elegans-acc-signal/prediction \
		${RESULTDIR}/elegans-don-signal/prediction/log ${RESULTDIR}/elegans-don-signal/prediction \
		${RESULTDIR}/elegans-cdsStop-signal/prediction/log ${RESULTDIR}/elegans-cdsStop-signal/prediction \
		${RESULTDIR}/elegans-cleave-signal/prediction/log ${RESULTDIR}/elegans-cleave-signal/prediction \
		${RESULTDIR}/elegans-intergenic-content/prediction/log ${RESULTDIR}/elegans-intergenic-content/prediction \
		${RESULTDIR}/elegans-utr5exon-content/prediction/log ${RESULTDIR}/elegans-utr5exon-content/prediction \
		${RESULTDIR}/elegans-cds_exon-content/prediction/log ${RESULTDIR}/elegans-cds_exon-content/prediction \
		${RESULTDIR}/elegans-intron-content/prediction/log ${RESULTDIR}/elegans-intron-content/prediction \
		${RESULTDIR}/elegans-utr3exon-content/prediction/log ${RESULTDIR}/elegans-utr3exon-content/prediction \
		${LSL_DIR}/log ${LSL_DIR} $GENE_TRAIN_OPTS #> ${LSL_DIR}/gene_train.log
	
	${BINDIR}/start_training '-' '-' ${LSL_DIR}
	
	touch $FN_DONE
fi

FN_DONE=$LSL_DIR/PRED_DONE
if test -f $FN_DONE; then 
	echo gene_predict already done
else
	echo 4b. using the trained gene predictor for prediction \[Log file in ${LSL_DIR}/gene_predict.log\]
	${BINDIR}/gene_predict ${GENOME_DIR} ${LSL_DIR}/log ${LSL_DIR} ${LSL_DIR}/gene_predict.log2 \
		${RESULTDIR}/elegans-tss-signal/prediction/log ${RESULTDIR}/elegans-tss-signal/prediction \
		${RESULTDIR}/elegans-tis-signal/prediction/log ${RESULTDIR}/elegans-tis-signal/prediction \
		${RESULTDIR}/elegans-acc-signal/prediction/log ${RESULTDIR}/elegans-acc-signal/prediction \
		${RESULTDIR}/elegans-don-signal/prediction/log ${RESULTDIR}/elegans-don-signal/prediction \
		${RESULTDIR}/elegans-cdsStop-signal/prediction/log ${RESULTDIR}/elegans-cdsStop-signal/prediction \
		${RESULTDIR}/elegans-cleave-signal/prediction/log ${RESULTDIR}/elegans-cleave-signal/prediction \
		${RESULTDIR}/elegans-intergenic-content/prediction/log ${RESULTDIR}/elegans-intergenic-content/prediction \
		${RESULTDIR}/elegans-utr5exon-content/prediction/log ${RESULTDIR}/elegans-utr5exon-content/prediction \
		${RESULTDIR}/elegans-cds_exon-content/prediction/log ${RESULTDIR}/elegans-cds_exon-content/prediction \
		${RESULTDIR}/elegans-intron-content/prediction/log ${RESULTDIR}/elegans-intron-content/prediction \
		${RESULTDIR}/elegans-utr3exon-content/prediction/log ${RESULTDIR}/elegans-utr3exon-content/prediction \
		${LSL_DIR}/prediction.gff3 ${LSL_DIR} 1 $GENE_TRAIN_OPTS #> ${LSL_DIR}/gene_predict.log

	touch $FN_DONE
fi
FN_DONE=$LSL_DIR/GENES_DONE
if test -f $FN_DONE; then 
	echo genes_from_prediction already done
else
	${BINDIR}/genes_from_predictions ${LSL_DIR} ${GENOME_DIR} 0
	${BINDIR}/write_gff3_wrapper ${LSL_DIR} '-' '-' genes

	touch $FN_DONE
fi

echo 4c. comparing the reference annotation and the gene prediction
${BINDIR}/gene_eval ${ANNO_DIR} ${LSL_DIR}/genome_wide_predictions | grep "level"

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo

