function [weights, loss, score, losses, scores, ok] = compute_weights(SEG, block, model, FN, parameters, check_result) ;
% weights=compute_weights(segmentation, block, model) ;

if nargin<5
  parameters = [];
end
if nargin<6
  check_result = 1 ;
end
ok=true ;
[penalty_array, content_weights, state_signals, a_trans] = gen_penalty_array(model, block.split, FN, parameters) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i=1:length(segmentation.genes)
%  segmentation.genes{i}(1,3)=model.segments.utr5exon ;
%  segmentation.genes{i}(end,3)=model.segments.utr3exon ;
%end ;

% simple state sequence -> generate this from segmentation
% if ~isfield(SEG,'path') | isfield(SEG,'pos_idx') 
%   [my_path,pos_idx,ok]=segmentation2path(SEG.segmentation, block, model);
% else
%   my_path = SEG.path ;
%   pos_idx = SEG.pos_idx
% end
my_path = SEG.path ;
pos_idx = SEG.pos_idx ;

% keyboard
% sg('send_command', 'loglevel WARNING') ;
% sg('send_command', 'loglevel ALL') ;
%model.loss.segments(:)=1 ;
%model.loss.segments(:,9:16)=2 ;
%keyboard

%[my_path,pos_idx,pos,ok]=segmentation2path(SEG, block, model) ;

seg_path = block.truth(1).seg_path;
% if the alternatives disagree, then don't incur a loss
for i=2:length(block.truth)
  if isempty(block.truth(i).seg_path), continue ; end ;
  idx = find(block.truth(i).seg_path(1,:) ~= block.truth(1).seg_path(1,:)) ;
  seg_path(2,idx) = 0 ;
end ;
if isfield(block, 'loss_mask')
	seg_path(2,:) = seg_path(2,:).*block.loss_mask(block.all_pos);
end
if iscell(block.features)
  features = full_features(block.features);
else 
  features = block.features;
end


if isfield(model, 'use_conservation') && model.use_conservation
  error('not yet implemented')
else
  use_orf=1;
  num_svms=8;
  use_long_transitions = 1;
  threshold = model.long_trans_thresh;
  long_transition_max_len = 1000000;
  sg_set_plif(penalty_array)
  sg('init_dyn_prog', num_svms)
  if isfield(block, 'content_pred')
    block.content_pred(end+1:num_svms,:) = deal(0);
    sg('set_lin_feat', block.seq, int32(block.all_pos-1), block.content_pred);
  else
    sg('precompute_content_svms',block.seq,int32(block.all_pos-1), content_weights');
  end
  sg('set_model', model.transition_pointers, use_orf, int32(model.mod_words), int32(state_signals),int32(model.orf_info))
  sg('set_feature_matrix',features)
  for j = 1:length(model.track_names)
    local_pos = find(~isnan(block.tracks(j,:)));
    plif_ids_name = sprintf('track_%i_plif_ids', j);
    sg('precompute_tiling_features',block.tracks(j,local_pos), int32(local_pos-1), int32(model.(plif_ids_name)-1))
  end
  for j = 1:length(model.segment_feature_names);
   sg('init_intron_list', int32(block.segment_lists{j}(:,1)-1)', int32(block.segment_lists{j}(:,2)-1)', int32(block.segment_scores{j})', int32(block.all_pos-1));
  end
  sg('long_transition_settings', use_long_transitions, threshold, long_transition_max_len);

  [p_deriv, q_deriv, a_deriv, penalty_deriv, scores, losses]= sg('best_path_trans_deriv', model.p', model.q', seg_path, a_trans, model.loss.segments, int32(my_path-1), int32(pos_idx-1));

end ;

loss = sum(losses) ;
score = sum(scores) ;

if isempty(parameters) && check_result
  if ~sum(scores)==0
    idx = find(scores)
    pos_idx(idx(1)+[-1:1]) 
    block.all_pos(pos_idx(idx(1)+[-1:1]))
    my_path(idx(1)+[-1:1])
    block.id
    error('invalid path: sum(scores) ~= 0')
    %weights = [];
    %ok=false ;
    %return
  end
else
  if check_result,
    assert(~isinf(sum(scores))) ;
    if isinf(sum(scores)), 
      ok=false; 
    end ;
  end ;
end

weights = init_weights(model) ;
for i=1:size(penalty_deriv,1)
  pen = penalty_array{i} ;
  weights = setfield(weights, pen.name, penalty_deriv(pen.id,:)) ;
end ;

%%% transitions
weights.transitions(:) = a_deriv(model.active_transitions) ;



%%%% for debugging
if 0 
  idx = find(isinf(scores))
  
  states = SEG.path(idx:idx+1)
  
  %%%check state

  features(SEG.path(idx+1),SEG.pos_idx(idx+1))
  %%%check transition 
  
  trans_lengths = SEG.pos(idx+1)-SEG.pos(idx)
  idx2 = find(ismember(SEG.segments(:,1:2),SEG.pos(idx:idx+1),'rows'));
  seg_type = SEG.segments(idx2,3)
 
  model.transition_pointers(states(2),states(1),:)
  model.lengths
 
  range = model.lengths_range

end
