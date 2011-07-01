function genes2label(genes_fname, spf_fname, spf_dir, genome_config_dir, signal_name, covered_genic_positions_only, fix_donor_offset) ;
% genes2label(genes_fname, spf_fname, spd_dir, genome_config_dir, signal_name, covered_genic_positions_only, fix_donor_offset) ;


if ~exist('fix_donor_offset'),
  fix_donor_offset = 1 ;
end ;
if ~exist('covered_genic_positions_only'),
  covered_genic_positions_only = 1 ;
end ;

disp('------------------------------------') ;
disp('Step 1: Setting up datastructures...') ;
disp('------------------------------------') ;

genome_config = [genome_config_dir '/genome.config'] ;
genome_info = init_genome(genome_config) ;

% load genes
load(genes_fname, 'genes') ;

% set strand, if missing
for i=1:length(genes),
  if ~isfield(genes, 'strand') || isempty(genes(i).strand),
    genes(i).strand = genes(i).strands(1) ;
  end ;
  % fix donor offset
  if fix_donor_offset, 
    for j=1:length(genes(i).exons)
      if genes(i).strand == '+',
        genes(i).exons{j}(:,2) = genes(i).exons{j}(:,2) + 1 ;
      else
        genes(i).exons{j}(:,1) = genes(i).exons{j}(:,1) - 1 ;
      end ;
    end 
  end ;
end ;
% add missing fields
genes = init_genes(length(genes)+1, [], genes) ;
genes = genes(1:end-1) ;

if 1,
  genes = process_gff3_genes(genes) ;
end ;
fprintf('Done.\n\n') ;


genes = filter_invalid_genes(genes,genome_info);
genes2label_helper(genes, genome_info, genome_config, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;

