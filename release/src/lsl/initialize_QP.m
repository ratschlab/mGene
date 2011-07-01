function QP = initialize_QP(blocks, method, model, ...
                            generate_neg_constraints, weight_names)
% QP= initialize_QP(blocks , method, model, generate_neg_constraints=1)
  

INF = method.INF;
C_regul = method.par_ms.C_regul; 
bins = method.plif_bins;
cnt_plif_detector = model.cnt_signals + model.cnt_contents  ;
cnt_plif_lengths = model.cnt_lengths  ;
cnt_transitions = model.cnt_transitions ;
max_neg = method.max_neg ;
margin = method.margin;

% just for older versions 
if ~isfield(model, 'track_names')
  model.track_names = {};
end
if ~isfield(model, 'segment_feature_names')
  model.segment_feature_names = {};
end

%%

num_pos_tracks = length(model.track_names);
num_seg_tracks = length(model.segment_feature_names);

%%%%%%%


content_names = fieldnames(model.contents) ;
id_contents = [] ;

for i=1:length(content_names),
  m = getfield(model.contents, content_names{i}) ;
  id_contents(end+1) = m;
  C_regul.id2C(m) = C_regul.contents;
end ;
signal_names = fieldnames(model.signals) ;
id_signals = [] ;
for i=1:length(signal_names),
  m = getfield(model.signals, signal_names{i}) ;
  id_signals(end+1) = m ;
  C_regul.id2C(m) = C_regul.signals;
end ;
length_names = fieldnames(model.lengths) ;
id_length = [] ;
for i=1:length(length_names),
  m = getfield(model.lengths, length_names{i}) ;
  id_length(end+1) = m ;
  C_regul.id2C(m) = C_regul.lengths;
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect plif ids for different tracks of position associated features and store 
% them in an cell array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:num_pos_tracks
  fieldname = sprintf('track_%i',j);
  track_fields = fieldnames(model.(fieldname));
  id_tracks{j} = [] ;
  for i=1:length(track_fields),
    m = getfield(model.(fieldname), track_fields{i}) ;
    id_tracks{j}(end+1) = m ;
    C_regul.id2C(m) = C_regul.tracks;
  end ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect plif ids for different segment associated features and store 
% them in an cell array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:num_seg_tracks
  fieldname = sprintf('segment_feature_%i',j);
  seg_fields = fieldnames(model.(fieldname));
  id_segment_features{j} = [] ;
  for i=1:length(seg_fields),
    m = getfield(model.(fieldname), seg_fields{i}) ;
    id_segment_features{j}(end+1) = m ;
    C_regul.id2C(m) = C_regul.segment_features;
  end ;
  fieldname = sprintf('segment_score_%i',j);
  seg_fields = fieldnames(model.(fieldname));
  id_segment_scores{j} = [] ;
  for i=1:length(seg_fields),
    m = getfield(model.(fieldname), seg_fields{i}) ;
    id_segment_scores{j}(end+1) = m ;
    C_regul.id2C(m) = C_regul.segment_features;
  end ;

end

%%%%%%
monoton_id = [id_contents id_signals] ;
nof_monoton = 0;
for m = monoton_id,
  param_ids = model.plif2param{m};
  nof_monoton = nof_monoton + length(param_ids) - 1;
end ;
smooth_id = [monoton_id id_length] ;

for j=1:num_pos_tracks
  nof_tracks(j) = 0;
  for m = id_tracks{j},
    param_ids = model.plif2param{m};
    nof_tracks(j) = nof_tracks(j) + length(param_ids) - 1;
  end ;
  smooth_id = [smooth_id id_tracks{j}] ;
end

for j=1:num_seg_tracks
  nof_seg_feat(j) = 0;
  for m = id_segment_features{j},
    param_ids = model.plif2param{m};
    nof_seg_feat(j) = nof_seg_feat(j) + length(param_ids) - 1;
  end ;
  smooth_id = [smooth_id id_segment_features{j}] ;
end

for j=1:num_seg_tracks
  nof_seg_score(j) = 0;
  for m = id_segment_scores{j},
    param_ids = model.plif2param{m};
    nof_seg_score(j) = nof_seg_score(j) + length(param_ids) - 1;
  end ;
  smooth_id = [smooth_id id_segment_scores{j}] ;
end

nof_smooth = 0;
for m = smooth_id,
  param_ids = model.plif2param{m};
  nof_smooth = nof_smooth + length(param_ids) - 1;
end ;

nof_trans_smooth = model.cnt_transitions ;
nof_smooth = nof_smooth + nof_trans_smooth ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define initial set of margin constraints 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear A A1 A1_sparse
A1 = zeros(10000, model.cnt_parameters) ;
A1_sparse = [] ;
A2 = sparse([],[],[],length(blocks)*max_neg, length(blocks)+nof_smooth,...
            length(blocks)*max_neg) ; 

q = 0 ; r = 0 ; q2 = 0 ;
num_complete = 0 ;


fprintf('initialize_QP: ') ;
for block_id = 1:length(blocks)
  if mod(block_id,50)==0, 
    fprintf('%i ... ', block_id) ;
  end ;
  block = blocks(block_id);
  if generate_neg_constraints,
    % change: max_path was not initialized before
    max_path = length(block.pred);
    for i=1:max_path,
      q = q+1 ; 
      q2 = q2+1 ; 
      A1(q2,1:model.cnt_parameters) = weights2vector(block.pred(i).weights, weight_names)-weights2vector(block.truth(1).weights, weight_names) ;    
      if all(abs(A1(q2,:))<1e-3)
        error('haeh')
      end ;
      if q2==size(A1,1)
        A1_sparse = [A1_sparse; sparse(A1)] ;
        q2 = 0 ;
        A1(:) = 0 ;
      end ;
      A2(q, block_id)= -1 ; 
    end ;
  end ;
end ;

A1 = [A1_sparse; sparse(A1(1:q2,:))] ;
assert(q==size(A1,1));
clear A1_sparse
A2 = A2(1:q,:) ;
A = [A1 A2] ;
clear A1 A2
b = -margin*ones(q,1) ;

% penalization of 1st derivatives
B = spzeros(2*nof_smooth, model.cnt_parameters+length(blocks)+nof_smooth, ...
            6*nof_smooth) ;
QP.f = zeros(model.cnt_parameters+length(blocks)+nof_smooth,1) ;
QP.f( model.cnt_parameters+(1:length(blocks)) ) = 1;
q = 0 ;
for m = smooth_id
  param_ids = model.plif2param{m};
  for i = param_ids(1:end-1) 
    q = q + 1;
    B(q+q-1,i) = 1 ;
    B(q+q-1,i+1)=-1 ;
    B(q+q-1,model.cnt_parameters+length(blocks)+q)=-1 ;
    B(q+q,i)   = -1 ;
    B(q+q,i+1) = +1 ;
    B(q+q,model.cnt_parameters+length(blocks)+q)=-1 ;
    QP.f(model.cnt_parameters+length(blocks)+q) = C_regul.id2C(m) ;
  end;
end ;

% add L1 constraints and penalty for transitions
param_ids = model.transitions2param ;
for i = param_ids
  q = q + 1;
  B(q+q-1,i) = 1 ;
  B(q+q-1,model.cnt_parameters+length(blocks)+q)=-1 ;
  B(q+q,i)   = -1 ;
  B(q+q,model.cnt_parameters+length(blocks)+q)=-1 ;
  QP.f(model.cnt_parameters+length(blocks)+q) = C_regul.transitions ;
end;
% figure; imagesc(B)
assert( q == nof_smooth );
A = [A;B]; 
b = [b; zeros(2*nof_smooth,1)] ;
clear B


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% monotonicity for all PLIFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = spzeros(nof_monoton, model.cnt_parameters+length(blocks)+nof_smooth, 2*nof_monoton) ;
q = 0 ;
for m = monoton_id,
  param_ids = model.plif2param{m};
  for i = param_ids(1:end-1)
    q = q+1 ;
    C(q,i) = 1 ;
    C(q,i+1) = -1 ;
  end ;
end ;
% figure;imagesc(C)
QP.A = [A;C]; 
QP.b = [b; zeros(nof_monoton,1)] ;

for j=1:num_pos_tracks
  T = spzeros(nof_tracks(j), model.cnt_parameters+length(blocks)+nof_smooth, 2*nof_tracks(j)) ;
  q = 0 ;
  for m = id_tracks{j},
    param_ids = model.plif2param{m};
    fieldname = sprintf('track_%i',j);
    names = fieldnames(model.(fieldname));
    for k=1:length(names), if model.(fieldname).(names{k})==m, break;end,end
    name = names{k};
    fieldname_monoton = sprintf('track_%i_monoton',j);
    if model.(fieldname_monoton).(name)==-1 %monotonical decreasing
      for i = param_ids(1:end-1)
        q = q+1 ;
        T(q,i)   = -1 ;
        T(q,i+1) = 1 ;
      end ;
    elseif model.(fieldname_monoton).(name)==1 %monotonical increasing
      for i = param_ids(1:end-1)
        q = q+1 ;
        T(q,i)   = 1 ;
        T(q,i+1) = -1 ;
      end ;
    end
  end ;
%  figure;imagesc(T)
  QP.A = [QP.A;T]; 
  QP.b = [QP.b; zeros(nof_tracks(j),1)] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segment lists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:num_seg_tracks
  T = spzeros(nof_seg_feat(j), model.cnt_parameters+length(blocks)+nof_smooth, 2*nof_seg_feat(j)) ;
  q = 0 ;
  for m = id_segment_features{j},
    param_ids = model.plif2param{m};
    fieldname = sprintf('segment_feature_%i',j);
    names = fieldnames(model.(fieldname));
    for k=1:length(names), if model.(fieldname).(names{k})==m, break;end,end
    name = names{k};
    fieldname_monoton = sprintf('segment_feature_%i_monoton',j);
    if model.(fieldname_monoton).(name)==-1 %monotonical decreasing
      for i = param_ids(1:end-1)
        q = q+1 ;
        T(q,i)   = -1 ;
        T(q,i+1) = 1 ;
      end ;
    elseif model.(fieldname_monoton).(name)==1 %monotonical increasing
      for i = param_ids(1:end-1)
        q = q+1 ;
        T(q,i)   = 1 ;
        T(q,i+1) = -1 ;
      end ;
    end
  end ;
%  figure;imagesc(T)
  QP.A = [QP.A;T]; 
  QP.b = [QP.b; zeros(nof_seg_feat(j),1)] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segment quality scores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:num_seg_tracks
  T = spzeros(nof_seg_score(j), model.cnt_parameters+length(blocks)+nof_smooth, 2*nof_seg_score(j)) ;
  q = 0 ;
  for m = id_segment_scores{j},
    param_ids = model.plif2param{m};
    fieldname = sprintf('segment_score_%i',j);
    names = fieldnames(model.(fieldname));
    for k=1:length(names), if model.(fieldname).(names{k})==m, break;end,end
    name = names{k};

    % assume monotonicity is equal to segment list
    fieldname_monoton = sprintf('segment_feature_%i_monoton',j);
    if model.(fieldname_monoton).(name)==-1 %monotonical decreasing
      for i = param_ids(1:end-1)
        q = q+1 ;
        T(q,i)   = -1 ;
        T(q,i+1) = 1 ;
      end ;
    elseif model.(fieldname_monoton).(name)==1 %monotonical increasing
      for i = param_ids(1:end-1)
        q = q+1 ;
        T(q,i)   = 1 ;
        T(q,i+1) = -1 ;
      end ;
    end
  end ;
%  figure;imagesc(T)
  QP.A = [QP.A;T]; 
  QP.b = [QP.b; zeros(nof_seg_score(j),1)] ;
end

QP.lb = [-INF*ones(model.cnt_parameters,1) ; zeros(length(blocks)+nof_smooth,1)] ;
QP.ub = [+INF*ones(model.cnt_parameters+length(blocks)+nof_smooth,1)]  ;
% QP.ub = [];

Q = spzeros(model.cnt_parameters+length(blocks)+nof_smooth) ;
for i=1:model.cnt_parameters,
  if( i < model.cnt_parameters-model.cnt_transitions )
    Q(i,i) = C_regul.plif_ys_sq;
  else
    Q(i,i) = C_regul.transitions_sq ;
  end;
end ;
for i=model.cnt_parameters+length(blocks)+(1:nof_smooth-nof_trans_smooth)
  Q(i,i) = Q(i,i)+C_regul.smoothness_sq ;
end ;
QP.Q = Q;


if isfield(method, 'fn_regul_solution') ,
  L=load(method.fn_regul_solution) ;
  fprintf('\nregularizing with C=%1.2f against solution (norm=%1.2f)\nfrom file %s \n', ...
          C_regul.solution, norm(L.QP.res(1:model.cnt_parameters)), ...
          method.fn_regul_solution) ;
  assert(L.TRAIN_PAR.model.cnt_parameters==model.cnt_parameters) ;
  for k=1:L.TRAIN_PAR.model.cnt_parameters
    QP.Q(k,k) = QP.Q(k,k) + C_regul.solution ;
    QP.f(k) = QP.f(k) - C_regul.solution*L.QP.res(k) ;
  end ;
end ;

%%   checks
r = zeros(size(QP.f)) ;
r(model.cnt_parameters+(1:length(blocks))) = 1;
assert(all(QP.A*r-QP.b==0)) ;

fprintf('Done \n') ;
