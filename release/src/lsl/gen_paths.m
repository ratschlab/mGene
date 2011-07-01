function [pred, dummy] = gen_paths(PAR)
% pred = gen_paths(PAR)
  
%sg_fn = which('sg')
%unix(sprintf('ldd %s', sg_fn))

% dummy return argument
dummy=[] ;

%%%%%%%%%%%%%%%%%

paths

% avoid multiple loading of the big constant parameter file
global GLOBAL_cPAR ;%GLOBAL_blocks;
if isempty(GLOBAL_cPAR),
  L = load(PAR.cPAR_fname, 'cPAR') ;
  GLOBAL_cPAR = L ;
  %blocks = load_struct(L.cPAR.blocks_fname) ; 
  %GLOBAL_blocks = blocks ;
else
  L = GLOBAL_cPAR ;
  %blocks = GLOBAL_blocks ;
end ;

%block = blocks(PAR.block_id) ;
if isfield(PAR, 'block')
  block = PAR.block ;
else
  load(PAR.block_name, 'block')
end



% fix some field ordering
L.cPAR.model = fix_model(L.cPAR.model) ;
L.cPAR.parameters = reorder_fields(L.cPAR.parameters, L.cPAR.parameters_fields) ;
block = fix_blocks(block, L.cPAR.model) ;

model = L.cPAR.model;
parameters = L.cPAR.parameters ;

%%%%%%%%%%%%%%%%%

if ~isfield(PAR, 'max_len')
  PAR.max_len=[] ;
end ;
if ~isfield(PAR, 'prediction')
  PAR.prediction=0 ;
end ;
block.seq(find(block.seq~='A'&block.seq~='C'&block.seq~='G'&block.seq~='T'))='A';


% subset the positions to approximate the full dynamic program
if isfield(PAR, 'approximate') && isequal(PAR.approximate, 1),

  %if keyboard_allowed(),
  %  keyboard
  %end ;

  % make sure all previously predicted positions, true positions and strong positions are in the list
  used_positions = [] ;
  for i=1:length(block.pred),
	  used_positions = [used_positions block.pred(i).pos] ;
  end ;
  true_positions = [] ;
  for i=1:length(block.truth),
    true_positions = [true_positions block.truth(i).pos] ;
  end ;
  strong_positions_acc = block.Signals.acc.pos(find(block.Signals.acc.Conf_cum>0.75)) ;
  strong_positions_don = block.Signals.don.pos(find(block.Signals.don.Conf_cum>0.75)) ;
  strong_positions = union(strong_positions_acc, strong_positions_don) ;
  grid_pos = block.all_pos(1:10:length(block.all_pos));
  strong_positions = union(strong_positions, grid_pos);
  good_positions = union(union(used_positions, strong_positions), true_positions) ;

  [tmp,idx_pos]=intersect(block.all_pos, good_positions) ;

  block.all_pos = block.all_pos(idx_pos) ;
  if isfield(block, 'pos')
	  block.pos = block.pos(idx_pos) ;
  end

  if iscell(block.features)
    for i=1:length(block.features)
      block.features{i} = block.features{i}(:, idx_pos) ;
    end ;
  else
    block.features = block.features(:, idx_pos, :) ;
  end ;
  block.content_pred = block.content_pred(:, idx_pos) ;

  % update seg_path in truth (this is all_pos based)
  for t=1:length(block.truth)
    block.truth(t) = gen_loss_matrix(block.truth(t), block.all_pos, block.gene_status, model) ;
  end ;
else
  idx_pos = 1:length(block.all_pos) ;
end ;

% generate full feature matrix, if still in sparse format
if iscell(block.features)
  block.features = full_features(block.features);
end

[penalty_array,content_weights,state_signals,a_trans] = gen_penalty_array(model, block.split, PAR.FN, parameters, PAR.max_len) ;

loss = model.loss.segments ;
if PAR.prediction,
  loss(:)=0 ;
  seg_path = zeros(2,length(block.all_pos)) ;
else
  %block.truth(1) = gen_loss_matrix(block.truth(1), block.pos, ...
  %                                 block.gene_confirmed, model) ;
  
  seg_path = block.truth(1).seg_path;
  % if the alternatives disagree, then don't incur a loss
  for i=2:length(block.truth)
    idx = find(block.truth(i).seg_path(1,:) ~= block.truth(1).seg_path(1,:)) ;
    seg_path(2,idx) = 0 ;
  end ;
end ;
if isfield(block, 'loss_mask')
	seg_path(2,:) = seg_path(2,:).*block.loss_mask(block.all_pos);
end

% SHOGUN interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t = cputime;
  use_orf = 1;
  num_svms = 8;
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
  sg('set_feature_matrix', block.features)
  for j = 1:length(model.track_names)
    local_pos = find(~isnan(block.tracks(j,:)));
    plif_ids_name = sprintf('track_%i_plif_ids', j);
    sg('precompute_tiling_features',block.tracks(j,local_pos), int32(local_pos-1), int32(model.(plif_ids_name)-1))
  end
  for j = 1:length(model.segment_feature_names);
   sg('init_intron_list', int32(block.segment_lists{j}(:,1)-1)', int32(block.segment_lists{j}(:,2)-1)', int32(block.segment_scores{j})', int32(block.all_pos-1));
  end
  if length(PAR.viterbi_nbest)<2,
    PAR.viterbi_nbest(2) = 0 ;
  end ;
  sg('long_transition_settings', use_long_transitions, threshold, long_transition_max_len) 
  [path_scores, path, ppos]= sg('best_path_trans', model.p', model.q', int32(PAR.viterbi_nbest), seg_path, a_trans, loss);
decode_time = cputime-t;

% very important check to assert that feature generation and veterbi decoding are equivalent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(isfield(PAR, 'skip_checks') && isequal(PAR.skip_checks,1)), 
  path_ = path(1,ppos(1,:)~=-1);
  ppos_ = ppos(1,ppos(1,:)~=-1);
  [p_deriv, q_deriv, a_deriv, penalty_deriv, scores, losses]= sg('best_path_trans_deriv', model.p', model.q', seg_path, a_trans, loss, int32(path_), int32(ppos_));
  %%% test
  A = model.A ;
  A(~isinf(A)) = parameters.transitions ;
  A(isinf(A)) = 0 ;
  score = sum(sum(a_deriv.*A))+sum(losses);
  for i=1:size(penalty_deriv,1)
    score = score + sum(penalty_deriv(i,:).*penalty_array{i}.penalties);
  end ;
  assert(all(losses>=0))
  if ~(abs(score-path_scores(1))<1e-4)
    fprintf('unexpected score differences: %1.3f %1.3f (diff=%1.3f)\n', ...
            score, path_scores, score-path_scores) ;
    if keyboard_allowed(),
      keyboard ;
    end ;
  end ;
end ;

%% clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sg('clean_up_dyn_prog');


true_path_weights = struct;
pred = []; loss=[] ;
for i=1:sum(PAR.viterbi_nbest)
  pred(i).loss = nan ;
  pred(i).path = path(i,ppos(i,:)~=-1)+1;
  pred(i).pos_idx = ppos(i,ppos(i,:)~=-1)+1 ;
  pred(i).pos = block.all_pos(pred(i).pos_idx);
  [pred(i).segments, pred(i).genes] = path2segmentation(pred(i).path, pred(i).pos_idx, block, model);
  if ~PAR.prediction
    [pred(i).weights, pred(i).loss, pred_x(i).score, pred_x(i).losses, pred_x(i).scores, pred_x(i).ok] = compute_weights(pred(i), block, model, PAR.FN, parameters) ;   
    if isfield(PAR, 'approximate') && isequal(PAR.approximate, 1),
      % the reason is that the predictions with long_transitions differ a little if the candidate lists are not the same. 
      % therefor the features of the true path has to be recomputed with this candidate list
      [tmp tmp2 idxx] = intersect(block.truth(1).pos, block.all_pos);
      block.truth(1).pos_idx = idxx;
      true_path_weights = compute_weights(block.truth(1), block, model, PAR.FN, parameters, 0) ;
    else
      true_path_weights = block.truth(1).weights;
    end ;
  end ;
  % undo position subseting in index to positions (done for approximation)
  pred(i).pos_idx = idx_pos(pred(i).pos_idx) ;
end
pred = orderfields(pred) ;

%disp(['saving pred to ' PAR.fn_save]) ;
cPAR = L.cPAR;
if isfield(PAR, 'block')
  % save mem
  PAR = rmfield(PAR, 'block');
end
save('-v7', PAR.fn_save,'decode_time', 'pred', 'true_path_weights', 'PAR')

return

