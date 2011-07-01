function anno2contentlabel(anno_fn, anno_dir, genome_dir, cpf_fname, cpf_dir, content_name, vs_content_names,intergenic_dist_to_genes) ;
% anno2label( anno_dir, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;

if ~exist('covered_genic_positions_only', 'var'),
  covered_genic_positions_only = 1 ;
end ;

disp('------------------------------------') ;
disp('Step 1: Loading genome annotation...') ;
disp('------------------------------------') ;

genes = load_genes(sprintf('%s/genes.mat', anno_dir)) ;

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
if strcmp(vs_content_names, 'all_other')
	vs_content_names = '';
end
if strcmp(intergenic_dist_to_genes, 'default')
	intergenic_dist_to_genes = '';
end

genes2content_helper(genes, genome_config, cpf_fname, cpf_dir, content_name, vs_content_names,intergenic_dist_to_genes) ;
