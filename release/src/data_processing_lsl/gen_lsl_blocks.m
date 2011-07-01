function [blocks genes split] = gen_lsl_blocks(genes, blocks, fn_genome_config, model, fn_split)
% [blocks genes split] = gen_lsl_blocks(genes, blocks, fn_genome_config, model, fn_split)


% remove blocks that contain no or incomplete genes
%--------------------------------------------------------------
block_rm_idx = [];
for j=1:length(blocks)
  gene_ids = blocks(j).gene_ids;
  if isempty(gene_ids)
    block_rm_idx = [block_rm_idx j];
  end
end
if 1
	fprintf('remove %i (%i) blocks containing no genes\n',length(block_rm_idx),length(blocks))
	blocks(block_rm_idx)=[];
else
	warning('do not remove blocks containing no genes');
end

block_rm_idx = [];
for j=1:length(blocks)
  gene_ids = blocks(j).gene_ids;
  for id=gene_ids
	assert(genes(id).id==id)
    if ~genes(id).is_complete
      block_rm_idx = [block_rm_idx j];
      break;
    end
  end
end
fprintf('remove %i (%i) blocks containing genes with is_complete flag set to zero\n',length(block_rm_idx),length(blocks))
blocks(block_rm_idx)=[];


%block_rm_idx = [];
%for j=1:length(blocks)
%  gene_ids = blocks(j).gene_ids;
%  for id=gene_ids
%	assert(genes(id).id==id)
%    if genes(id).start<=blocks(j).start + 10
%      	block_rm_idx = [block_rm_idx j];
%      	break;
%    end
%    if genes(id).stop>=blocks(j).stop - 10
%      	block_rm_idx = [block_rm_idx j];
%      	break;
%    end
%  end
%end
%
%fprintf('remove %i (%i) blocks where genes reach block boundaries aftern UTR inference \n',length(block_rm_idx),length(blocks))
%blocks(block_rm_idx)=[];

clear curr_genes curr_blocks s_blocks_idx s_genes_idx c_blocks_idx c_genes_idx s j k g block_rm_idx gene_idxs



% filter out blocks that contain overlapping genes
%------------------------------------------------------------
fprintf('filter out blocks that contain overlapping genes\n')
rm_idx = [];
for j=1:length(blocks)
  if length(blocks(j).gene_ids) <= 1, continue; end
  starts = [genes(blocks(j).gene_ids).start];
  stops = [genes(blocks(j).gene_ids).stop];
  [starts sort_idx] = sort(starts);
  stops = stops(sort_idx);
  for k=2:length(starts)
    if ~(starts(k)>stops(k-1))
      rm_idx = [rm_idx j];
      break; 
    end
  end
end
fprintf('remove %i (%i) blocks containing overlapping genes\n',length(rm_idx),length(blocks))
blocks(rm_idx) = [];
clear rm_idx sort_idx starts stops

% add segmentation to genes
%-----------------------------------------------------------
genes(1).segments{1} = [];
tmp.genes = genes;
tmp.segment_ids = model.segments;
genes = gen_gene_segmentation(tmp);
clear tmp

%rm_idx = [];
%for j = 1:length(genes)
%  tr_keep_idx = [];
%  for k = 1:length(genes(j).transcripts)
%    if isempty(genes(j).segments)
%      rm_idx = [rm_idx j];
%      continue; 
%    end
%    if ~isempty(genes(j).segments{k})
%      tr_keep_idx = [tr_keep_idx k];
%    end
%  end
%  if length(tr_keep_idx) == 0
%     rm_idx = [rm_idx j];
%  else
%    % remove those transcripts
%    genes(j).transcripts = genes(j).transcripts{tr_keep_idx};
%    genes(j).transcript_status = genes(j).transcript_status(tr_keep_idx);
%    genes(j).transcript_valid = genes(j).transcript_valid(tr_keep_idx);
%    genes(j).transcript_complete = genes(j).transcript_complete(tr_keep_idx);
%    genes(j).exons = genes(j).exons{tr_keep_idx};
%    genes(j).cds_exons = genes(j).cds_exons{tr_keep_idx};
%    genes(j).utr5_exons = genes(j).utr5_exons{tr_keep_idx};
%    genes(j).utr3_exons = genes(j).utr3_exons{tr_keep_idx};
%    genes(j).tis = genes(j).tis{tr_keep_idx};
%    genes(j).cdsStop = genes(j).cdsStop{tr_keep_idx};
%    genes(j).cleave = genes(j).cleave{tr_keep_idx};
%    genes(j).tss = genes(j).tss{tr_keep_idx};
%  end
%end
%fprintf('remove %i (%i) blocks where segmentation failed\n',length(rm_idx),length(blocks)) 
%blocks(rm_idx) = [];


% add segmentation to blocks
%-------------------------------------------------------------
[blocks.config] = deal(fn_genome_config);
blocks = load_sequence(blocks);
blocks = gen_block_segmentation(genes,blocks,model);

% get split from training blocks
%-------------------------------------------------------------
split = [];
if ~isfield(blocks,'split')||isempty(blocks(1).split)
  l = load(fn_split,'split');
  split = -1*ones(1,length(blocks));
  for j=1:length(blocks), 
    for k=1:length(fieldnames(l.split)), 
      idx = find(l.split.(['split_' num2str(k)])==blocks(j).id); 
      if ~isempty(idx), 
        split(j)=k;
        blocks(j).split=k; 
      end, 
    end, 
  end
  assert(isequal(unique(split), 1:5))
  clear l j k idx
end
% remove blocks with empty truth 
%-------------------------------------------------------------
fprintf(1,'filter out blocks with empty truth\n')
rm_idx = [];
val_idx = [];
for b = 1:length(blocks)
  if isempty(blocks(b).truth)
    rm_idx = [rm_idx b];
    continue;
  end
end
blocks(rm_idx) = [];
fprintf('remove %i (%i) blocks where generation of truth has failed\n',length(rm_idx), length(blocks))


