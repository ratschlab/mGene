#!/bin/bash
samtoolsdr="/fml/ag-raetsch/share/software/samtools/"

for m in `ls *.bam`
do
	echo "starting with file $m"
	# BAM must be sorted by start position to use random access
	$samtoolsdr./samtools sort $m ${m}_sorted
	new=`echo ${m}_sorted.bam | sed 's/bam_sorted/sorted/'`
	mv ${m}_sorted.bam $new
	# index for BAM file
	$samtoolsdr./samtools index $new
done
