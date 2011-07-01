function blocks = gen_features(blocks, model)
% blocks = gen_features(block, model)

%blocks=blocks(1);

signal_names = fieldnames(model.signals) ;
signal_ids = zeros(1,length(signal_names)) ;
for i=1:length(signal_names) 
  signal_ids(i) = model.signals.(signal_names{i}) ;
end ;

%for i=1:length(model.states)
%  motifs{i}=[];
%  if ~isempty(model.states(i).consensus)
%    for j=1:length(model.states(i).consensus)
%      motifs{i}{j}=lower(model.states(i).consensus{j});
%    end
%  end
%end

blocks(1).features = {[]};
blocks(1).all_pos = [];
fprintf('Processing %i blocks: ', length(blocks)) ;
for b_idx = 1:length(blocks)
  if mod(b_idx, 50)==0,
    fprintf('%i ... ', b_idx);
  end ;
  block = blocks(b_idx);
  all_pos= [] ;
  for i=1:length(signal_names) 
    new_pos = block.Signals.(signal_names{i}).pos ;
    all_pos = [all_pos new_pos(:)'] ;
  end ;
  all_pos = unique(all_pos) ;

  thirddim =0;
  for i=1:length(model.states)
    thirddim=max(thirddim,length(model.states(i).signal));
  end ;

  features = -inf(model.cnt_states, length(all_pos), thirddim) ;
  genestr  = lower(block.seq) ;
  
  for idx_states=1:length(model.states) 
    state = model.states(idx_states);
    if isempty(state.signal)
      continue ; 
    end ;
    for j=1:length(state.signal)
      sig_name = signal_names{find(state.signal(j)==signal_ids)} ;
      sig_pos_temp{j} = block.Signals.(sig_name).pos ;
      sig_scores_temp{j} = block.Signals.(sig_name).Conf_cum ;
      assert(isequal(sort(sig_pos_temp{j}), sig_pos_temp{j}));
      assert(length(sig_pos_temp{j})==length(unique(sig_pos_temp{j})));
      % [sig_pos_temp{j},iidx] = unique(sig_pos_temp{j}) ;
      % sig_scores_temp{j} = sig_scores_temp{j}(iidx) ;
    end
    sig_pos=sig_pos_temp{1};
    for j=2:length(state.signal)
      sig_pos=intersect(sig_pos, sig_pos_temp{j}) ;
    end ;
    sig_scores =[];
    for j=1:length(state.signal)
      sig_scores{j} = -inf(size(sig_pos));
      [tmp,idx1,idx2] = intersect(sig_pos,sig_pos_temp{j});
      sig_scores{j}(idx1) = sig_scores_temp{j}(idx2) ;
      assert(length(idx1)==length(sig_scores{j}));
    end ;clear sig_pos_temp sig_scores_temp

    if ~isempty(state.consensus) 
      for k=1:length(state.consensus) 
		motifs{k}=lower(state.consensus{k});
		%motifl(k)=length(motifs{k});
      end
      cons_found = check_consensus(genestr, motifs, state.consensus_pos, sig_pos'-1);
      clear motifs;
       %if sum(cons_found==0)~=0
       %  keyboard
       %end
      %assert(all(cons_found)) ;
      sig_pos(~cons_found)=[] ;
      for j=1:length(sig_scores)
        sig_scores{j}(~cons_found)=[] ;
      end
    end ;
  
    if ~isempty(state.non_consensus) 
      for k=1:length(state.non_consensus)
        motifs{k}=lower(state.non_consensus{k});
      end
      ncons_found = check_consensus(genestr, motifs, state.non_consensus_pos, sig_pos'-1);
      clear motifs
      sig_pos(ncons_found~=0)=[] ;
      for j=1:length(sig_scores)
        sig_scores{j}(ncons_found~=0)=[] ;
      end
    end ;

    [tmp,idx1,idx2]=intersect(all_pos,sig_pos) ;
    assert(length(idx2)==length(unique(sig_pos))) ;
    for j=1:length(sig_scores)
      features(idx_states, idx1, j) = sig_scores{j}(idx2) ;
    end
  end ; %%loop over states

  % first position only start-state allowed
  features(:, 1, :) = -inf ;
  features(model.state_ids.seq_start, 1, 1) = 0 ;

  % last position only end-state allowed
  features(:, end, :) = -inf ;
  features(model.state_ids.seq_end, end, 1) = 0 ;

  %% every position can be used as rna_seq state
  if isfield(model.state_ids,'rna_seq_polya')
    features(model.state_ids.rna_seq_polya, :, 1) = 0 ;
  end
  block.all_pos = all_pos;
  block.features = sparse_features(features); % teuer
  blocks(b_idx)=block;
  
  %%%check don
  state_ids = model.state_ids.don;
  num=0;
  for j=state_ids
    num=num+1;
    cons{num} = model.states(j).consensus;
    nocons{num} = model.states(j).non_consensus;
  end
  for j=2:length(cons)
    for k=1:j-1
      if isequal(cons{j},cons{k})&&isequal(nocons{j},nocons{k})
        %%!!!
        assert(isequal(features(state_ids(j),:,:),features(state_ids(k),:,:)))
      end
    end
  end
  
  %%%check acc
  state_ids = model.state_ids.acc;
  num=0;
  for j=state_ids
    num=num+1;
    cons{num} = model.states(j).consensus;
    nocons{num} = model.states(j).non_consensus;
  end
  for j=2:length(cons)
    for k=1:j-1
      if isequal(cons{j},cons{k})&&isequal(nocons{j},nocons{k})
        assert(isequal(features(state_ids(j),:,:),features(state_ids(k),:,:)))
      end
    end
  end
  %assert(all(any(~isinf(features(:,:,1)))))
end
fprintf('Done.\n') ;

return
%features(features>20)=-1 ;
%figure;imagesc(features(:,:,1)) ;
%figure;imagesc(features(:,:,2)) ;



