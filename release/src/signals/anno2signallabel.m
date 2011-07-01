function anno2signallabel(anno_fname, anno_dir,genome_dir, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;
% anno2label(anno_fname, anno_dir,genome_dir, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;

if ~exist('covered_genic_positions_only', 'var'),
  covered_genic_positions_only = 1 ;
end ;

disp('------------------------------------') ;
disp('Step 1: Loading genome annotation...') ;
disp('------------------------------------') ;

genes = load_genes(sprintf('%s/genes.mat', anno_dir)) ;
genes = update_splicegraph(genes);


fprintf('Loaded annotation with %i genes\n', length(genes)) ;

disp('-------------------------------------') ;
disp('Step 2: Loading genome information...') ;
disp('-------------------------------------') ;

genome_config = [genome_dir '/genome.config'] ;
genome_info = init_genome(genome_config);

if length(genes)==0 
  error('There are no genes in your annotation gene structure (AGS). Probably something went wrong during parsing');
end ;

if isempty(genes(1).chr_num)
	genes = get_chr_num(genes, genome_info)
end
% genes = filter_invalid_genes(genes,genome_info);

genes2label_helper(genes, genome_info, genome_config, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;
