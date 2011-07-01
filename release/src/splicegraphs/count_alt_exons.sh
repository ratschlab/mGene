#!/bin/tcsh 
#
#$ -cwd
#$ -j y
#$ -N stats
#$ -l arch=lx24-amd64
#$ -l h_vmem=3.5G
#$ -l matlab=1
####$ -R y

matlab -nojvm -nodisplay -r "count_alt_exons('altsplicecount.txt');exit"
