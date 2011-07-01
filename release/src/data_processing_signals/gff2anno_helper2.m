function [genes, FAULTS1, FAULTS2] = gff2anno_helper2(gff_fname, save_fname, par_fname, log_fname, genome_config)
% [genes,FAULTS1,FAULTS2] = gff2anno_helper2(gff_fname, save_fname, par_fname, log_fname, genome_config)
paths

load(par_fname)

genome_info = init_genome(genome_config);

fd=fopen(log_fname, 'w+') ;
disp('-----------------------------') ;
disp('Step 3: Parsing gff3 file ...') ;
disp('-----------------------------') ;

genes = init_genes(0) ;
assert(isempty(genes)) ;
for i=1:length(cname)
  [g,F] = parse_gff3(gff_fname, '', source, cname{i}, TYPE, log_fname);
  genes = [genes g];
  FAULTS1(i) = F;
  clear g F
end

fprintf('Read %i gene structures from GFF file\n', length(genes)) ;
fprintf(fd, 'Read %i gene structures from GFF file\n\n', length(genes)) ;
fclose(fd) ;

disp('-----------------------------') ;
disp('Step 4: Processing parsed genes ...') ;
disp('-----------------------------') ;

% save_struct(sprintf('%s_genes_parse.mat', save_fname), genes, 'genes') ;

[genes,FAULTS2] = process_gff3_genes(genes, genome_info, level3_coding_names) ;

if do_correct_tis_stop,
  disp('-----------------------------') ;
  disp('Step 5a: correcting TIS and Stop positions ...') ;
  disp('-----------------------------') ;
  
  genes = correct_tis_stop(genes, genome_info);
end ;

disp('-----------------------------') ;
disp('Step 5b: checking processed genes ...') ;
disp('-----------------------------') ;

which_checks.exons_sorted = 1;
which_checks.intron_length = 1;
which_checks.splicesites = 1;
which_checks.orf = 1;
which_checks.gene_length = 1;
which_checks.graph = 0;
which_checks.transacc = 0;
which_checks.complete = 1;

genes = check_genes(genes, genome_info,which_checks);

% save_struct(sprintf('%s_genes_check.mat',save_fname), genes, 'genes') ;

all_genes = genes;
save_struct(sprintf('%sgenes_all.mat', save_fname),all_genes, 'all_genes') ;
save_struct(sprintf('%sFAULTS1.mat', save_fname),FAULTS1, 'FAULTS1') ;
save_struct(sprintf('%sFAULTS2.mat', save_fname),FAULTS2, 'FAULTS2') ;

return

disp('-------------------------------------') ;
disp('Step 6: Filter out invalid genes');
disp('-------------------------------------') ;

genes = filter_invalid_genes(genes, genome_info,log_fname);

disp('-------------------------------------') ;
disp('Step 7: curate genes');
disp('-------------------------------------') ;

genes = curate_anno(genes, genome_info, log_fname);
fprintf('Done.\n\n') ;

genes = check_genes(genes, genome_info,which_checks);

disp('-----------------------------') ;
disp('Step 8: build splicegraph ...') ;
disp('-----------------------------') ;

% build splicegraphs
genes = build_splice_graph_caller(genes) ;
genes = infer_splice_graph_caller(genes) ;
fprintf('Sorting exons...\n');
if length(genes)==1 && isempty(genes(1).exons)
  genes(1)=[];
else
  for i=1:length(genes)
    [dummy,exon_order] = sort(genes(i).splicegraph{1}(1,:),2,'ascend');
    genes(i).splicegraph{1} = genes(i).splicegraph{1}(:,exon_order);
    genes(i).splicegraph{2} = genes(i).splicegraph{2}(exon_order,exon_order);
  end
  genes = alt_const(genes);
end

fprintf('Done.\n\n') ;
%annotation_fname

save_struct(sprintf('%sgenes.mat', save_fname), genes, 'genes') ;


