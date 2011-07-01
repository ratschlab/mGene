function galaxy_anno2gff(gene_dir, gff_fname)
% galaxy_anno2gff(gene_dir, gff_fname)

addpath('~/svn/projects/mGene_core/webserver/');
addpath('~/svn/tools/utils/');

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf('-------------------------------------------------- \n');
fprintf('Anno2GFF started at %s\n', timedate) ;
fprintf('-------------------------------------------------- \n\n');

load(sprintf('%s/genes.mat', gene_dir), 'genes');
write_gff3(genes, gff_fname, 'galaxy', 'galaxy');

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf('\n-------------------------------------------------- \n');
fprintf('Anno2GFF finished at %s\n', timedate) ; 
fprintf('--------------------------------------------------- \n');