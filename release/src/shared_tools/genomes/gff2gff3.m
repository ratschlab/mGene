function [] = gff2gff3(infile, genome_config_fname, outfile, confirm_status, exclude_regions, rename_double_genes)
% function [] = gff2gff3(infile_gff2,genome_info,outfile_gff3,confirm_status, exclude_regions, rename_double_genes)
  if nargin < 4
    error('to few parameters');
  end
  if nargin < 5
    exclude_regions=[];
  end
  if nargin < 6
    rename = false;
  end
  % test if infile is existent
  if isempty(infile)
    error('no infile specified');
  elseif ~fexist(infile)
    error('infile not existent or maybe wrong path');
  end
  % test if file 'save_fname' already exists
  if isempty(outfile)
    error('no outfile specified');
  elseif fexist(outfile)
    warning('file %s.mat already exists. Will be overwritten!\n',outfile);
  end
  genome_info=init_genome(genome_config_fname);
  genes = parse_transcripts_gff2(infile, [], genome_info);
  if ~isempty(exclude_regions)
    % filter out some regions
    genes=filter_regions_from_genes(genes, exclude_regions);
  end
  % handle needed fields for merge
  for i=1:length(genes)
    genes(i).exons_confirmed{1}=[];
    genes(i).transcript_info(1)=2;
  end
  contigs=unique([genes.chr_num]);
  merged_genes=[];
  for i=1:length(contigs)
    part=genes(find([genes.chr_num]==contigs(i)));
    genes=genes(find([genes.chr_num]~=contigs(i)));
    part=merge_transcripts_by_colocation(part,genome_info);
    merged_genes=[merged_genes part];
  end
  if rename_double_genes && length(merged_genes) ~= length(unique({merged_genes.name}))
    % rename double gene identifiers
    sorted=unique({merged_genes.name});
    for i=1:length(sorted)
      ids=find(strcmp({merged_genes.name},sorted{i}));
      if length(ids)==1
        continue;
      end
      for j=1:length(ids)
	merged_genes(ids(j)).name=[merged_genes(ids(j)).name char(64+j)];
      end
    end 
  end
  [tmp,source]=unix(sprintf('head -1 %s | cut -f 2',infile));
  write_gff3(merged_genes,outfile,deblank(source),confirm_status);
end
