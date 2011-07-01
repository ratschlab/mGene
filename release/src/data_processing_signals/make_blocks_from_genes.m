function blocks = make_blocks_from_genes(genes, fn_config, margin, merge_blocks, min_dist)
% blocks = make_blocks_from_genes(genes, fn_config, margin=2000, merge_blocks=1, min_dist=100)

%fn_config = '/fml/ag-raetsch/share/databases/genomes//F_graminearum/F_graminearumOct2007/genebuild_noproteins/genome.config';
genome_info = init_genome(fn_config);

if ~exist('margin', 'var')
  margin = 2000;
end ;
if ~exist('merge_blocks', 'var')
  merge_blocks = 1 ;
end ;
if ~exist('min_dist', 'var')
  min_dist = 100;
end ;

%blocks = genes;
%blocks = struct('strand',{genes.strand},'chr',{genes.chr},'chr_num',{genes.chr_num},'start',{genes.start},...
%		'stop',{genes.stop},'gene_ids',{genes.id},'paralogs',{genes.paralogs},'is_alt',{genes.is_alt},...
%		'gene_complete',{genes.is_complete}, 'config' ,fn_config );
blocks = struct('strand',{genes.strand},'chr',{genes.chr},'chr_num',{genes.chr_num},'start',{genes.start},...
		'stop',{genes.stop},'gene_ids',{genes.id},'is_alt',{genes.is_alt},...
		'gene_complete',{genes.is_complete}, 'config' ,fn_config );

for i=1:length(blocks)
  blocks(i).chr_num = strmatch(blocks(i).chr, genome_info.contig_names, 'exact') ;
  assert(length(blocks(i).chr_num)==1) ;
end ;

% ADJUST BOUNDARIES
for j=1:length(blocks)
  margin_left = margin ;
  margin_right = margin ;

  if isfield(genes(j), 'next_gene_left') && ~isempty(genes(j).next_gene_left),
    if margin_left > round((genes(j).start-genes(j).next_gene_left)/2),
      margin_left = round((genes(j).start-genes(j).next_gene_left)/2) ;
    end ;
  end ;
  if isfield(genes(j), 'next_n_left') && ~isempty(genes(j).next_n_left),
    if margin_left > genes(j).start-genes(j).next_n_left,
      margin_left = genes(j).start-genes(j).next_n_left ;
    end ;
  end ;
  if isfield(genes(j), 'next_gene_right') && ~isempty(genes(j).next_gene_right),
    if margin_right > round((genes(j).next_gene_right-genes(j).stop)/2),
      margin_right = round((genes(j).next_gene_right-genes(j).stop)/2) ;
    end ;
  end ;
  if isfield(genes(j), 'next_n_right') && ~isempty(genes(j).next_n_right),
    if margin_right > genes(j).next_n_right-genes(j).stop,
      margin_right = genes(j).next_n_right-genes(j).stop ;
    end ;
  end ;

  blocks(j).gene_ids = genes(j).id;
  blocks(j).paralogs = [];
  blocks(j).gene_status = max(genes(j).transcript_status); 
  blocks(j).id = j; 
  if blocks(j).start<1
    if keyboard_allowed()
      keyboard
    end
  end
  flat_fn = genome_info.flat_fnames{blocks(j).chr_num};
  d = dir(flat_fn);
  len = d.bytes;
 
  if blocks(j).stop+margin_right<=len
    blocks(j).stop=blocks(j).stop+margin_right;
  else
    assert(blocks(j).stop<=len)
    blocks(j).stop=len;
  end
  if blocks(j).start-margin_left>1
    blocks(j).start=blocks(j).start-margin_left;
  else
    assert(blocks(j).start>=1)
    blocks(j).start=1;
  end
end

if merge_blocks,
  % MERGE BLOCKS
  len = length(blocks)
  for s='+-'
    strands = [blocks.strand];
    str_blocks = blocks(strands==s);
    x = unique({str_blocks.chr});
    for c = 1:length(x)
      c_idx = find(strcmp({str_blocks.chr},x(c)));
      contig_blocks = str_blocks(c_idx);
      %keyboard
      starts = [contig_blocks.start];
      stops  = [contig_blocks.stop];
      [starts sort_idx] = sort(starts);
      stops = stops(sort_idx);
      contig_blocks = contig_blocks(sort_idx);
      dist = starts(2:end)-stops(1:end-1);
      merge_idx = find(dist<min_dist);
      
      rm_idx = [];
      for j=merge_idx;
        contig_blocks(j+1)=merge(contig_blocks(j),contig_blocks(j+1));
        %keyboard
        rm_idx = [rm_idx j];
      end 
      contig_blocks(rm_idx) = []; 
      blocks(c_idx) = [];
      blocks = [blocks contig_blocks];
    end
  end
  len = length(blocks)
end ;


for j=1:length(blocks)
  blocks(j).id = j;
  off = [blocks(j).start-1 blocks(j).stop+1];
  blocks(j).offset = off((blocks(j).strand=='-')+1); 
end

return

function block = merge(block1,block2)
  block = block1;
  block.gene_ids = [block1.gene_ids block2.gene_ids];
  block.paralogs = unique([block1.paralogs block2.paralogs]);
  block.gene_status = [block1.gene_status block2.gene_status];
  block.start = min([block1.start block2.start]);
  block.stop = max([block1.stop block2.stop]);
return


