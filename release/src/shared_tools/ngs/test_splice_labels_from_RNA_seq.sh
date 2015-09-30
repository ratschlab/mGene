#!/bin/bash

#bamfile1=/cbio/grlab/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_left_sam_stranded.mapped.2.bam
#bamfile2=/cbio/grlab/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_right_sam_stranded.mapped.2.bam

bamfile1=/cbio/grlab/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans.dummy/polyA_left_trim_new.bam
bamfile2=/cbio/grlab/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans.dummy/polyA_right_trim_new.bam


genome_config_dir=/cbio/grlab/nobackup/projects/rgasp/mgene_predictions/elegans/genome_dir/
genome_config=/cbio/grlab/nobackup/projects/rgasp/mgene_predictions/elegans/genome_dir/genome.config

output_dir=/cbio/grlab/home/jonas/tmp/RNA_seq_splice_label_unbiased

time ./splice_labels_from_RNA_seq $genome_config $output_dir $bamfile1 $bamfile2


matlab -r "addpath ~/svn/projects/mGene_core/; paths; addpath ~/svn/projects/mGene_core/RNA_seq_label; train_splice_signals($genome_config_dir, $output_dir)"
