function blocks_new=merge_blocks(blocks, model, num_blocks_per_blocks) ;
% blocks=merge_blocks(blocks, model, num_blocks_per_blocks) ;

%idx=randperm(length(blocks)) ;
%idx = reshape(idx(1:(floor(length(blocks)/num_blocks_per_blocks)*num_blocks_per_blocks)), floor(length(blocks)/num_blocks_per_blocks), num_blocks_per_blocks) ;

rand('seed', 1234);

% merge blocks, that have similar intergenic segment lengths. this tries not to spoil the intergenic length  distribution
left_len = zeros(1,length(blocks)) ;
right_len = zeros(1,length(blocks)) ;
for i=1:length(blocks)
  left_len(i) = blocks(i).truth(1).segments(1,2)-blocks(i).truth(1).segments(1,1) ;
  right_len(i) = blocks(i).truth(1).segments(end,2)-blocks(i).truth(1).segments(end,1) ;
end 

taken=zeros(1,length(blocks)) ;
for i=1:floor(length(blocks)/num_blocks_per_blocks),
  iidx=find(taken==0) ;
  iidx=iidx(randperm(length(iidx))) ;
  idx(i,1) = iidx(1) ;
  taken(iidx(1))=1 ;
  for j=2:num_blocks_per_blocks,
    iidx=find(taken==0) ;
    [tmp,iidx2]=min(abs(left_len(iidx)-right_len(idx(i,j-1)))) ;
    idx(i,j) = iidx(iidx2) ;
    taken(iidx(iidx2))=1 ;
    %[tmp abs(left_len(idx(i,j))-right_len(idx(i,j-1)))]
  end ;
end ;
    
empty_block = blocks([]);
empty_block(1).split = [];
signal_names=fieldnames(blocks(1).Signals) ;

blocks_new=blocks([]) ;
for i=1:size(idx,1),
  blocks_new(i).id = i ;
  blocks_new(i).seq = '' ;
  blocks_new(i).orig_block_ids = idx(i,:) ;
  blocks_new(i).orig_offsets = [] ;
  blocks_new(i).truth.segments = [] ;
  blocks_new(i).is_alt = 0 ;
  blocks_new(i).gene_complete = [] ;
  blocks_new(i).gene_ids = [] ;

  for s=1:length(signal_names),
    blocks_new(i).Signals.(signal_names{s})=struct ;
    blocks_new(i).Signals.(signal_names{s}).pos = [] ;
    blocks_new(i).Signals.(signal_names{s}).contig_pos = [] ;
    blocks_new(i).Signals.(signal_names{s}).Conf_cum = [] ;
  end 
  cnt=0 ;
  for j=idx(i,:),
    cnt = cnt + 1 ;
    %if length(blocks(j).truth)>1,
    %  warning('merge_blocks:multiple_truths', 'ignoring multiple truths') ;
    %end ;

    blocks_new(i).orig_offsets(end+1) = length(blocks_new(i).seq) ;
    if cnt==1,
      blocks_new(i).start = blocks(j).start ;
    end ;
    if cnt == num_blocks_per_blocks,
      blocks_new(i).stop = blocks(j).stop + blocks_new(i).orig_offsets(end) ;
    end ;

    % content predictions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if isfield(blocks, 'content_pred') && cnt==1
      blocks_new(i).content_pred = [ blocks(j).content_pred ] ;
    elseif isfield(blocks, 'content_pred')
      blocks_new(i).content_pred = [ blocks_new(i).content_pred blocks(j).content_pred+blocks_new(i).content_pred(end) ] ;
    end

    % tracks of additional features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if isfield(blocks, 'tracks') && cnt==1
      blocks_new(i).tracks = blocks(j).tracks;
    elseif isfield(blocks, 'tracks')
      blocks_new(i).tracks = [blocks_new(i).tracks blocks(j).tracks];
    end

    % lists of segment information and scores
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if isfield(blocks, 'segment_lists') && cnt==1
      blocks_new(i).segment_lists = blocks(j).segment_lists;
      blocks_new(i).segment_scores = blocks(j).segment_scores;
    elseif isfield(blocks, 'segment_lists')
      for k = 1:length(blocks_new(i).segment_lists)
        blocks_new(i).segment_lists{k} = [ blocks_new(i).segment_lists{k}; blocks(j).segment_lists{k}+blocks_new(i).orig_offsets(end) ] ;
        blocks_new(i).segment_scores{k} = [blocks_new(i).segment_scores{k}; blocks(j).segment_scores{k}];
	assert(size(blocks_new(i).segment_lists{k}, 1) == size(blocks_new(i).segment_scores{k}, 1))
      end
    end

	% mask to switch of loss in alternative regions	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	if isfield(blocks, 'loss_mask') && cnt==1
		blocks_new(i).loss_mask = blocks(j).loss_mask; 
	elseif isfield(blocks, 'loss_mask')
		blocks_new(i).loss_mask = [blocks_new(i).loss_mask blocks(j).loss_mask];
	end
		

    blocks_new(i).seq = [ blocks_new(i).seq blocks(j).seq ] ;
    blocks_new(i).is_alt = blocks_new(i).is_alt | blocks(j).is_alt ;

    blocks_new(i).gene_status = [blocks_new(i).gene_status blocks(j).gene_status] ;
    blocks_new(i).gene_ids = [blocks_new(i).gene_ids blocks(j).gene_ids] ;
    blocks_new(i).gene_complete = [blocks_new(i).gene_complete blocks(j).gene_complete] ;
    blocks_new(i).truth.segments = [  blocks_new(i).truth.segments;
                        [blocks(j).truth(1).segments(:,1:2)+blocks_new(i).orig_offsets(end) blocks(j).truth(1).segments(:,3) blocks(j).truth(1).segments(:,4)*cnt] ] ;


    %signal_names=fieldnames(blocks(j).Signals) ;
    for s=1:length(signal_names),
      blocks_new(i).Signals.(signal_names{s}).pos = [blocks_new(i).Signals.(signal_names{s}).pos;
                          blocks(j).Signals.(signal_names{s}).pos+blocks_new(i).orig_offsets(end)] ;
      blocks_new(i).Signals.(signal_names{s}).Conf_cum = [blocks_new(i).Signals.(signal_names{s}).Conf_cum;
                          blocks(j).Signals.(signal_names{s}).Conf_cum] ;
    end 
	blocks(j) = empty_block;
  end ;

  % remove multiple intergenic segments between two genes
  while (1),
    truth_type = blocks_new(i).truth.segments(:,3) ;
    idx_t = find(truth_type(1:end-1)==model.segments.intergenic & truth_type(2:end)==model.segments.intergenic, 1, 'first') ;
    if isempty(idx_t), 
      break ;
    end ;
    S= blocks_new(i).truth.segments ;
    blocks_new(i).truth.segments = [S(1:idx_t-1, :);
                        S(idx_t,1) S(idx_t+1, 2:end);
                        S(idx_t+2:end, :) ] ;
  end ;
  % fix column 4 of segmentation if there was more than one gene in each block 
  gene_count=0;
  for line=1:size(blocks_new(i).truth.segments,1)
    if blocks_new(i).truth.segments(line,4)==0
      gene_count = gene_count+1;
    else
      blocks_new(i).truth.segments(line,4) = gene_count;
    end
  end

end ;


