function blocks = initialize_blocks(blocks, model, max_neg, FN)

% blocks = initialize_blocks(blocks, model)
% - transforms conservation to DNA_alphabet 
% - generates weights and loss for truth
% - adds decoy paths 

%%%%%%%

if isfield(model, 'use_conservation') && model.use_conservation 
  disp('Converting conservation levels') ;
  for id = 1:length(blocks),
    blocks(id).genestr_conserv = conservation2DNA_alphabet(blocks(id).conservation) ;
  end ;  
end ;


%%%%%%
rm_idx = [];
blocks(1).pred = struct;
fprintf('Processing %i blocks: ', length(blocks)) ;
for block_id = 1:length(blocks)
  if mod(block_id,50)==0, 
    fprintf('%i ... ', block_id) ;
  end ;
  block = blocks(block_id);
  block.truth(1).seg_path = [] ;
  block.seq=upper(block.seq);
  notacgt_idx = find(block.seq~='A'&block.seq~='C'&block.seq~='G'&block.seq~='T');
  if length(notacgt_idx)>0
    fprintf('Change %i non-[ACGT] symbold to ''A''\n',length(notacgt_idx))
    block.seq(notacgt_idx)='A';
  end
  for t=1:length(block.truth)
    %block.truth(t) = gen_loss_matrix(block.truth(t), block.all_pos, block.gene_confirmed, model) ;
    block.truth(t) = gen_loss_matrix(block.truth(t), block.all_pos, block.gene_status, model) ;
	try
    	[weights,loss,score,losses,scores] = compute_weights(block.truth(t), block, model, FN,[],1);
	catch
		fprintf('Discard block %i becase compute_weights returned -inf\n',block_id);
		rm_idx = [rm_idx, block_id];
		continue
	end
    if t==1, assert(loss==0) ; end
    if isempty(weights)
      %block.truth(t).weights = blocks(1).truth(1).weights ; % hack
      error('empty weights') ;
    else
      block.truth(t).weights = weights ;
    end
  end

  if ~isempty(rm_idx)&&rm_idx(end)==block_id
	continue
  end

  % generate wrong paths
  block.pred = gen_decoys(block, model, max_neg) ;
  % draw_pred(block.truth,block.pred,model)
  max_path = length(block.pred) ;
  good_idx = [] ;
  for i = 1:max_path,
    [weights,loss] = compute_weights(block.pred(i), block, model,FN,[],0);
    if isempty(weights)
      block.pred(i) = [];
      warning('invalid decoy path')
      max_path = max_path-1;
    else
      if loss>0,
        good_idx(end+1) = i ;
        assert(loss~=0) ;
        block.pred(i).weights = weights ;
        block.pred(i).loss = loss ;
      end ;
    end
  end
  block.pred = block.pred(good_idx) ;
  blocks(block_id) = block ;
  clear weights ;
end ;
blocks(rm_idx) = [];
fprintf('Done.\n') ;
