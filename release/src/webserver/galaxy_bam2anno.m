function galaxy_bam2anno(dir_config, bam_dir, output_dir)
% galaxy_bam2anno(dir_config, bam_files, output_dir)

addpath('~/svn/tools/genomes/');
addpath('~/svn/tools/utils/');
addpath('~/svn/projects/mGene_core/webserver/');
addpath('~/svn/projects/mGene_core/RNA_seq_label/');
addpath('~/svn/projects/mGene_core/auxiliary_data/');
addpath('~/svn/projects/mGene_core/data_processing_signals/');
addpath('~/svn/projects/mGene_core/signals/');
addpath('~/svn/projects/mGene_core/utils/');
addpath('~/svn/projects/mGene_core/tools/');
addpath('~/svn/projects/RNASeq_galaxy/rdiff.web/');


[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf('-------------------------------------------------- \n');
fprintf('BAM2Anno started at %s\n', timedate) ;
fprintf('-------------------------------------------------- \n\n');

label_gen(sprintf('%s/genome.config', dir_config), sprintf('%s/alignments.bam', bam_dir), output_dir);

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf('\n-------------------------------------------------- \n');
fprintf('BAM2Anno finished at %s\n', timedate) ; 
fprintf('--------------------------------------------------- \n');