function genes = filter_invalid_genes(genes, genome_info,log_fname)
%  genes = filter_invalid_genes(genes, genome_info,log_fname)
  
if exist('log_fname')
  fd=fopen(log_fname, 'a') ;
end

ignore_missing_referenced_entries = 1;  
ignore_genes_without_transcripts = 1;
ignore_genes_without_exons = 1;
ignore_genes_without_chr = 1;

take_map.is_valid = ones(1,length(genes)) ;
take_map.no_parents = ones(1,length(genes)) ;
take_map.no_transcripts = zeros(1,length(genes)) ;
take_map.no_exons = ones(1,length(genes)) ;
take_map.no_chr = ones(1,length(genes)) ;

%--------------

if (length([genes.is_valid]) ~= length(genes)) 
	which_checks = struct; % use default values
	genes = check_genes(genes, genome_info, which_checks);
end
take_map.is_valid(find(~[genes.is_valid]))=0 ;

if isfield(genes, 'is_correctly_gff3_referenced')&&~isempty(genes(1).is_correctly_gff3_referenced)
	assert(length([genes.is_correctly_gff3_referenced]) == length(genes)) ;
	take_map.no_parents(find(~[genes.is_correctly_gff3_referenced]))=0 ;
end

num_invalid_transcripts = 0 ;

for i=1:length(genes),
  % find genes with no transcript
  take_map.no_transcripts(i) = ~isempty(genes(i).transcripts) ;
  % find genes with no exons 
  empty = 1;
  for j=1:length(genes(i).exons),
    if ~isempty(genes(i).exons{j})
      empty=0 ;
      break
    end ;
  end ;
  if empty,
    take_map.no_exons(i) = 0 ;
  end ;
  % find genes of unknown contig
  chr_num = find(ismember(lower(genome_info.contig_names), lower(genes(i).chr))) ;
  take_map.no_chr(i) = length(chr_num)==1 ;
  
  % remove invalid transcripts
  num_invalid_transcripts = num_invalid_transcripts + sum(~genes(i).transcript_valid) ;
  idx = find(genes(i).transcript_valid) ;
  genes(i).exons = genes(i).exons(idx) ;
  if isfield(genes, 'exons_confirmed') && ~isempty(idx) && length(genes(i).exons_confirmed)>=max(idx)
  	genes(i).exons_confirmed = genes(i).exons_confirmed(idx);
  end
  genes(i).cds_exons = genes(i).cds_exons(idx) ;
  genes(i).transcripts = genes(i).transcripts(idx) ;
  genes(i).transcript_valid = genes(i).transcript_valid(idx) ;
  try
    genes(i).transcript_complete = genes(i).transcript_complete(idx) ; 
  catch
  end
  if isfield(genes, 'transcript_status') && ~isempty(genes(i).transcript_status)
    genes(i).transcript_status = genes(i).transcript_status(idx) ;
  end  
  genes(i).utr5_exons = genes(i).utr5_exons(idx) ;
  genes(i).utr3_exons = genes(i).utr3_exons(idx) ;
  if ~isempty(genes(i).tis),
    genes(i).tis = genes(i).tis(idx) ;
  end ;
  if ~isempty(genes(i).cdsStop),
    genes(i).cdsStop = genes(i).cdsStop(idx) ;
  end ;
  try,
    genes(i).tss = genes(i).tss(idx) ;
  catch
    % warning('tss elements missing') ;
  end ;
  try,
    genes(i).cleave = genes(i).cleave(idx) ;
  catch
    % warning('cleave elements missing') ;
  end ;
  try,
    genes(i).polya = genes(i).polya(idx) ;
  catch
    % warning('polya elements missing') ;
  end ;
end ;

fprintf(1, ' %i valid genes.\n', sum(take_map.is_valid==1)) ;
fprintf(1, ' %i invalid genes.\n', sum(take_map.is_valid==0)) ;
fprintf(1, '      %i genes with no exon segments.\n', sum(take_map.no_exons==0)) ;
fprintf(1, '      %i genes with no transcripts.\n', sum(take_map.no_transcripts==0)) ;
fprintf(1, '      %i genes incorrectly referenced (parent gene missing).\n', sum(take_map.no_parents==0)) ;
fprintf(1, '      %i genes on unknown contigs/chromosomes.\n', sum(take_map.no_chr==0)) ;
fprintf(1, '      %i invalid transcripts found.\n',num_invalid_transcripts);
if exist('fd')
fprintf(fd,'\n')
fprintf(fd, ' %i valid genes.\n', sum(take_map.is_valid==1)) ;
fprintf(fd, ' %i invalid genes.\n', sum(take_map.is_valid==0)) ;
fprintf(fd, '      %i genes with no exon segments.\n', sum(take_map.no_exons==0)) ;
fprintf(fd, '      %i genes with no transcripts.\n', sum(take_map.no_transcripts==0)) ;
fprintf(fd, '      %i genes incorrectly referenced (parent gene missing).\n', sum(take_map.no_parents==0)) ;
fprintf(fd, '      %i genes on unknown contigs/chromosomes.\n', sum(take_map.no_chr==0)) ;
fprintf(fd, '      %i invalid transcripts found.\n',num_invalid_transcripts);
fclose(fd);
end

take_idx = ones(1,length(genes));
take_idx = take_idx&take_map.is_valid;

if ignore_missing_referenced_entries,
  take_idx = take_idx&take_map.no_parents;
end ;
if ignore_genes_without_transcripts
  take_idx = take_idx&take_map.no_transcripts;
end ;
if ignore_genes_without_exons
  take_idx = take_idx&take_map.no_exons;
end ;
if ignore_genes_without_chr
  take_idx = take_idx&take_map.no_chr;
end ;
  
genes = genes(take_idx) ;
  
      
fprintf('Done.\n\n') ;

