function [penalty_array,content_weights,state_signals,a_trans] = gen_penalty_array(model,cont_subset_id, FN,  parameters, max_len) ;
% [penalty_array,content_weights,state_signals,a_trans] = gen_penalty_array(model,cont_subset_id, parameters, max_len) ;
  
if nargin<4 || isempty(parameters)
  use_param = 0 ;
else
  use_param = 1 ;
end ;

content_names = fieldnames(model.contents) ;
%idx=sort_numeric_fields(model.contents) ;
%content_names = content_names(idx) ;

length_names = fieldnames(model.lengths) ;
%idx=sort_numeric_fields(model.lengths) ;
%length_names = length_names(idx) ;

signal_names = fieldnames(model.signals) ;
%idx=sort_numeric_fields(model.signals) ;
%signal_names = signal_names(idx) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate penalty_array. set plif penalties to zero 

if nargin<5 || isempty(max_len)
  max_len = 0 ;
  for s = 1:model.cnt_lengths,
    len = getfield(model.lengths_range, length_names{s}) ;
    if len(2)>max_len,
      max_len = len(2) ;
    end ;
  end ;
end ;

penalty_array = {} ;
content_weights = [] ;
num_svms = 0 ;
for s = 1:model.cnt_contents,
  pen.id        = getfield(model.contents, content_names{s}) ;
  assert(pen.id == length(penalty_array)+1) ;
  pen.name      = ['contents_' content_names{s}] ;
  pen.limits    = getfield(model.boundaries.contents, content_names{s}) ;
  pen.limits    = pen.limits(1:end-1) ;
  if ~use_param
    pen.penalties = zeros(1,length(pen.limits)) ;
  else
    pen.penalties = getfield(parameters, pen.name) ;
    pen.penalties = pen.penalties(:)' ;
  end ;
  % keyboard
  pen.transform = '' ;
  pen.min_value   = -1e6 ;
  pen.max_value   = 1e6 ;
  pen.use_cache = 0 ;
  num_svms = num_svms + 1 ;
  pen.use_svm   = num_svms ;
  if strcmp(content_names{s},'dummy')
    content_weights(num_svms, :) = deal(0);
  else
    %fn_content_all_models = [FN.output_cont.(content_names{s}).fn_SVMs '_partition=' num2str(cont_subset_id) '.mat'];
    if ~isempty(FN) && isfield(FN, 'output_cont') && isfield(FN.output_cont, content_names{s})
      fn_content_all_models = FN.output_cont.(content_names{s}).fn_SVMs;
    else
      fn_content_all_models = '';
    end
    if ~fexist(fn_content_all_models) || exist(fn_content_all_models, 'dir') 
      %fprintf('file %s with content svms not found, use zeros as weights\n', fn_content_all_models);
    elseif ~ismember('w', who_file(fn_content_all_models)), 
      %fprintf('no weight vector found in file %s, use zeros as weights\n', fn_content_all_models);
    else
      l = load(fn_content_all_models, 'w');
      content_weights(num_svms, :) = l.w;
    end
    clear l
  end
  penalty_array{end+1} = pen ;
end ;

if size(content_weights,1)<8
  % current shogun requirement
  content_weights(end+1:8,:)=0 ;
end ;

for s = 1:model.cnt_lengths,
  pen.id        = getfield(model.lengths, length_names{s}) ;
  assert(pen.id == length(penalty_array)+1) ;
  pen.name      = ['lengths_' length_names{s}] ;
  pen.limits    = getfield(model.boundaries.lengths, length_names{s}) ;
  pen.limits    = pen.limits(1:end-1) ;
  if ~use_param
    pen.penalties = zeros(1,length(pen.limits)) ;
  else
    pen.penalties = getfield(parameters, pen.name) ;
    pen.penalties = pen.penalties(:)' ;
  end ;
  pen.transform = 'log(+1)' ;
  len=getfield(model.lengths_range, length_names{s}) ;
  pen.min_value   = len(1) ;
  pen.max_value   = min(len(2), max_len) ;
  pen.use_cache = 1 ;
  pen.use_svm   = 0 ;
  penalty_array{end+1} = pen ;
end ;

for s = 1:model.cnt_signals,
  pen.id        = getfield(model.signals, signal_names{s}) ;
  assert(pen.id == length(penalty_array)+1) ;
  pen.name      = ['signals_' signal_names{s}] ;
  pen.limits    = getfield(model.boundaries.signals, signal_names{s}) ;
  pen.limits    = pen.limits(1:end-1) ;
  pen.penalties = zeros(1,length(pen.limits)) ;
  if ~use_param
    pen.penalties = zeros(1,length(pen.limits)) ;
  else
    pen.penalties = getfield(parameters, pen.name) ;
    pen.penalties = pen.penalties(:)' ;
  end ;
  pen.transform = '' ;
  pen.min_value   = -1e6 ;
  pen.max_value   = 1e6 ;
  pen.use_cache = 0 ;
  pen.use_svm   = 0 ;
  penalty_array{end+1} = pen ;
end ;

for j = 1:length(model.track_names)
  fieldname = sprintf('track_%i',j);
  names = fieldnames_sorted(model.(fieldname));
  for s = 1:model.cnt_tracks(j),
    pen.id        = getfield(model.(fieldname), names{s}) ;
    if ~(pen.id == length(penalty_array)+1)
	fprintf('name:%s, pen.id:%i, length(penalty_array)+1:%i',[fieldname '_' names{s}], pen.id, length(penalty_array)+1);
	error('pen.id ~= length(penalty_array)+1');
    end
    pen.name      = [fieldname '_' names{s}] ;
    pen.limits    = getfield(model.boundaries.(fieldname), names{s}) ;
    pen.limits    = pen.limits(1:end-1) ;
    if ~use_param
      pen.penalties = zeros(1,length(pen.limits)) ;
    else
      pen.penalties = getfield(parameters, pen.name) ;
      pen.penalties = pen.penalties(:)' ;
    end ;
    % keyboard
    pen.transform = '' ;
    pen.min_value   = -1e6 ;
    pen.max_value   = 1e6 ;
    pen.use_cache = 0 ;
    num_svms = num_svms + 1 ;
    pen.use_svm   = num_svms ;
    penalty_array{end+1} = pen ;
  end ;
end

for j = 1:length(model.segment_feature_names)
  fieldname = sprintf('segment_feature_%i',j);
  names = fieldnames_sorted(model.(fieldname));
  for s = 1:model.cnt_segment_features(j),
    pen.id        = getfield(model.(fieldname), names{s}) ;
    if ~(pen.id == length(penalty_array)+1)
	fprintf('name:%s, pen.id:%i, length(penalty_array)+1:%i',[fieldname '_' names{s}], pen.id, length(penalty_array)+1);
	error('pen.id ~= length(penalty_array)+1');
    end
    pen.name      = [fieldname '_' names{s}] ;
    pen.limits    = getfield(model.boundaries.(fieldname), names{s}) ;
    pen.limits    = pen.limits(1:end-1) ;
    if ~use_param
      pen.penalties = zeros(1,length(pen.limits)) ;
    else
      pen.penalties = getfield(parameters, pen.name) ;
      pen.penalties = pen.penalties(:)' ;
    end ;
    % keyboard
    pen.transform = '' ;
    pen.min_value   = -1e6 ;
    pen.max_value   = 1e6 ;
    pen.use_cache = 0 ;
    num_svms = num_svms + 1 ;
    pen.use_svm   = num_svms ;
    penalty_array{end+1} = pen ;
  end ;

  fieldname = sprintf('segment_score_%i',j);
  names = fieldnames(model.(fieldname));
  for s = 1:model.cnt_segment_scores(j),
    pen.id        = getfield(model.(fieldname), names{s}) ;
    if ~(pen.id == length(penalty_array)+1)
	fprintf('name:%s, pen.id:%i, length(penalty_array)+1:%i',[fieldname '_' names{s}], pen.id, length(penalty_array)+1);
	error('pen.id ~= length(penalty_array)+1');
    end
    pen.name      = [fieldname '_' names{s}] ;
    pen.limits    = getfield(model.boundaries.(fieldname), names{s}) ;
    pen.limits    = pen.limits(1:end-1) ;
    if ~use_param
      pen.penalties = zeros(1,length(pen.limits)) ;
    else
      pen.penalties = getfield(parameters, pen.name) ;
      pen.penalties = pen.penalties(:)' ;
    end ;
    % keyboard
    pen.transform = '' ;
    pen.min_value   = -1e6 ;
    pen.max_value   = 1e6 ;
    pen.use_cache = 0 ;
    num_svms = num_svms + 1 ;
    pen.use_svm   = num_svms ;
    penalty_array{end+1} = pen ;
  end ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_signals=0;
for i=1:model.cnt_states
  num_signals=max([num_signals length(model.states(i).signal)]);
end

state_signals=zeros(model.cnt_states, num_signals) ;
for i=1:model.cnt_states,
  state_signals(i,1:length(model.states(i).signal))=model.states(i).signal ;
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A = model.A ;
if use_param
  A(~isinf(A)) = parameters.transitions ;
  model.A = A;
end

a_trans = model.a_trans ;
k=0 ;
for i=1:size(A,1)
  idx = find(~isinf(A(i,:))) ;
  val = A(i,idx) ;
  assert(all(a_trans(1,k+1:k+length(idx))==i-1)) ;
  assert(all(a_trans(2,k+1:k+length(idx))==idx-1)) ;  
  a_trans(3,k+1:k+length(idx))=val ;    
  k=k+length(idx) ;
end ;
a_trans=a_trans' ;
[tmp,idx]=sort(a_trans(:,2)) ;
a_trans = a_trans(idx,:) ;

function names = fieldnames_sorted(structure)
  names = fieldnames(structure);
  for j = 1:length(names)
    assert(isnumeric(structure.(names{j})))
    id(j) = structure.(names{j});
  end
  [tmp idx] = sort(id);
  names = names(idx);
return

function idx=sort_numeric_fields(s)
% idx=sort_numeric_fields(s)

names=fieldnames(s) ;
values=[] ;
for i=1:length(names),
  values(i) = s.(names{i}) ;
end ;
[tmp,idx]=sort(values) ;

