#!/bin/bash


### access the last argument
#eval lfile=\$$#
#
#echo "---------------------------------"
#echo $lfile
#echo "---------------------------------"

#for arg in $* ; do
#	echo "$arg" 
#done

# this script calls a matlab/octave function to set the LD_LIBRARY_PATH to the same value 
# as the matlab path is set later
`/home/galaxy/svn/projects/mGene_core/webserver/set_library_path.sh`

export MOSEKLM_LICENSE_FILE=~/svn/tools/utils_octave/mosek_interface/src/mosek.dat; 
export LD_LIBRARY_PATH=~/svn/tools/utils_octave/mosek_interface/src/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/galaxy/svn/tools/utils_octave/mosek_interface/src/linux64x86/bin/:$LD_LIBRARY_PATH


tmpfile=`tempfile`

oct=/home/galaxy/software/octave-3.0.4/bin/octave 

addpaths="rmpath('/home/galaxy/mgene_galaxy'); addpath('/home/galaxy/svn/projects/mGene_core/'); paths; addpath('/home/galaxy/svn/tools/utils_octave/mosek_interface/inst/'); addpath('/home/galaxy/svn/tools/utils_octave/mosek_interface/src');"
gio=$1; shift
anno=$1; shift
annod=""

#signals
tss=$1; shift
tssd=$1; shift
tis=$1; shift
tisd=$1; shift
acc=$1; shift
accd=$1; shift
don=$1; shift
dond=$1; shift
cdsStop=$1; shift
cdsStopd=$1; shift
cleave=$1; shift
cleaved=$1; shift

# contents
intergenic=$1; shift
intergenicd=$1; shift
utr5exon=$1; shift
utr5exond=$1; shift
cds_exon=$1; shift
cds_exond=$1; shift
intron=$1; shift
intrond=$1; shift
utr3exon=$1; shift
utr3exond=$1; shift

# output
Trained_Gene_Predictor=$1; shift
Trained_Gene_Predictord=$1; shift



# concat agruments for train options to a single string
train_opts=""
if [ $# != 1 ]; then
	train_opts="$1"
	shift
fi

until [ -z $1 ] ; do
	if [ $# != 1 ]; then
		train_opts="$train_opts;$1"
		shift
	else
		# last argument is the log file
		Log_File=$1
		shift
	fi
done

echo
echo "---------------------------------"
echo annotation: $anno
echo output dir: $Trained_Gene_Predictord
echo logfile:    $Log_File
echo "---------------------------------"
echo

# label generation
$oct --eval "$addpaths gene_train('$gio', '$anno', '$annod', '$tss', '$tssd', '$tis', '$tisd', '$acc', '$accd', '$don', '$dond', '$cdsStop', '$cdsStopd', '$cleave', '$cleaved', '$intergenic', '$intergenicd','$utr5exon', '$utr5exond','$cds_exon', '$cds_exond','$intron', '$intrond','$utr3exon', '$utr3exond','$Trained_Gene_Predictor', '$Trained_Gene_Predictord', '$train_opts');" 2>$tmpfile >>$Log_File; 

# training
PAR_file_name=$Trained_Gene_Predictord/lsl/data/training_PAR.mat
echo
echo "--------------------"
echo training parameters
echo $PAR_file_name
echo "--------------------"
echo
$oct --eval "$addpaths load('$PAR_file_name', 'PAR'); [fn_predictor, iter] = train_path_caller(PAR);" 2>$tmpfile >>$Log_File; 

cat $tmpfile | ~/bin/output_filter.sh | grep error 1>&2
cat $tmpfile | ~/bin/output_filter.sh | grep -v error >> $Log_File
rm $tmpfile

~/bin/cleanup_output.sh $Log_File
