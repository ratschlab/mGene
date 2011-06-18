#!/bin/bash
## USAGE: GFFParser.sh <GFF file> <Result file name in .mat extension> <gio>

sh `dirname $0`/../mgene_config.sh
export PYTHONPATH=~/lib/lib/python2.6/dist-packages/:~/lib/python2.5/site-packages/:~/lib/python/:/fml/ag-raetsch/home/vipin/svn/projects/genefinding/parsegff/
export PYTHONPATH=/fml/ag-raetsch/home/vipin/lib/python2.5/site-packages/:/fml/ag-raetsch/home/vipin/lib/python/:$PYTHONPATH
echo 'GFF file parsing started...'
## %%% Python GFF parsing started%%%
python `dirname $0`/GFFParser.py $1 $2
## %%% GFF Parsing done.%%%  
oct=octave
## %%% Rearranging gene models%%% 
echo "addpath $MGENE_SRC_PATH; paths; arrange_genes('$2');quit;"
$oct --eval "addpath $MGENE_SRC_PATH; paths; arrange_genes('$2');quit;"

## USAGE: CurateGenes <GIO-genome_info_path> <gene_models_file-Created by GFFParser.sh program>
oct=octave
$oct --eval "addpath $MGENE_SRC_PATH; paths; CurateGenes('$3', '$2'); quit;"

