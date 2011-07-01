function gff2label(gff_fname, spf_fname, spf_dir, genome_config_dir, coding_flag, signal_name, ignore_missing_referenced_entries, covered_genic_positions_only) ;
% gff2label(gff_fname, spf_fname, spf_dir, genome_config_dir, coding_flag, signal_name, ignore_missing_referenced_entries, covered_genic_positions_only) ;

if ~exist('ignore_missing_referenced_entries'),
  ignore_missing_referenced_entries = 1 ;
end ;
if ~exist('covered_genic_positions_only'),
  covered_genic_positions_only = 1 ;
end ;

disp('------------------------------------') ;
disp('Step 1: Setting up datastructures...') ;
disp('------------------------------------') ;

genome_config = [genome_config_dir '/genome.config'] ;

needs_coding = ~isempty(strmatch(signal_name, {'tis', 'stop'})) ;

while (1),

  % preparing data structures
  source = '' ;
  level1_names = {'gene'} ;
  level2_names = {'mRNA'} ;
  if coding_flag,
    level3_names = {'five_prime_UTR', 'three_prime_UTR', 'CDS', 'coding_exon'} ;
  else
    level3_names = {'exon'} ;
  end ;
  TYPE.short_names = {level1_names{:} level2_names{:} level3_names{:}} ;
  TYPE.IDs = TYPE.short_names ;
  TYPE.hierachies = [ones(1,length(level1_names)) 2*ones(1,length(level2_names)) 3*ones(1,length(level3_names))] ;
  
  genome_info = init_genome(genome_config) ;
  fprintf('Done.\n\n') ;
  
  disp('-----------------------------') ;
  disp('Step 2: Parsing gff3 file ...') ;
  disp('-----------------------------') ;
  
  
  [ret, nofgenes] = unix(['grep gene ' gff_fname '| wc -l']);
  assert(ret==0)
  nofgenes=str2num(nofgenes);
  
  genes = init_genes(max(nofgenes,1));
  len = 0 ;
  % split by contig
  for contig=1:length(genome_info.contig_names)
    cname = genome_info.contig_names{contig};
    new_genes = parse_gff3(gff_fname, '', '', cname, TYPE) ;
    genes(len+1:len+length(new_genes))=new_genes;
    len = len + length(new_genes);
    clear new_genes
  end 
  genes = genes(1:len) ;
  fprintf(1, 'Read %i gene structures from gff file\n', length(genes)) ;

  take_map = ones(1,length(genes)) ;
  for i=1:length(genes),
    empty = 1 ;
    for j=1:length(genes(i).exons),
      if ~isempty(genes(i).exons{j})
        empty=0 ;
      end ;
    end ;
    if empty,
      take_map(i) = 0 ;
    end ;
  end ;
  if sum(take_map==0)>0,
    fprintf(1, 'Ignoring %i genes with no exon seqments.\n', sum(take_map==0)) ;
    genes = genes(take_map==1) ;
  end ;

  if length(genes)==0 && coding_flag && ~needs_coding,
    warning('Did not find any genes with coding information in file. \nSwitching to ignore the coding regions for this signal (%s).', signal_name);
    coding_flag = 0 ;
    pause(5) ;
  else
    break ;
  end ;
end ;

if ignore_missing_referenced_entries,
  if sum(~[genes.is_correctly_gff3_referenced])>0,
    fprintf(1, 'Ignoring %i genes with missing references in the GFF3 file.\n', sum(~[genes.is_correctly_gff3_referenced]) ) ;
  end ;
  assert(length([genes.is_correctly_gff3_referenced]) == length(genes)) ;
  genes = genes(find([genes.is_correctly_gff3_referenced])) ;

  take_map = zeros(1,length(genes)) ;
  for i=1:length(genes),
    take_map(i) = ~isempty(genes(i).transcripts) ;
  end  ;
  if sum(take_map==0)>0,
    fprintf(1, 'Ignoring %i genes with no transcripts.\n', sum(take_map==0)) ;
    genes=genes(find(take_map)) ;
  end ;
end ;
genes = process_gff3_genes(genes) ;
if isequal(level3_names, {'exon'}),
  for i=1:length(genes),
    genes(i).tis = [] ;
    genes(i).cdsStop = [] ;
  end ;
end ;
fprintf('Done.\n\n') ;

genes = filter_invalid_genes(genes,genome_info);

genes2label_helper(genes, genome_info, genome_config, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;
