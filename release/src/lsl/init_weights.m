function weights=init_weights(model) ;
% weights=init_weights(model) ;

weights = struct;

%%% content detectors
content_names = fieldnames(model.contents) ;
for s = 1:model.cnt_contents,
  weights = setfield(weights, ['contents_' content_names{s}], zeros(1,model.bins));
end

%%% length detectors
length_names = fieldnames(model.lengths) ;
for s = 1:model.cnt_lengths,
  weights = setfield(weights, ['lengths_' length_names{s}], zeros(1,model.bins));
end

%%% signal detectors
signal_names = fieldnames(model.signals) ;
for s = 1:model.cnt_signals,
  weights = setfield(weights, ['signals_' signal_names{s}], zeros(1,model.bins));
end

%% positons associated features
for j=1:length(model.track_names);
  fieldname = sprintf('track_%i',j);
  names = fieldnames(model.(fieldname)) ;
  for s = 1:model.cnt_tracks(j),
    weights = setfield(weights, [fieldname '_' names{s}], zeros(1,model.bins));
  end
end

%%% segment associated features
for j=1:length(model.segment_feature_names);
  fieldname = sprintf('segment_feature_%i',j);
  names = fieldnames(model.(fieldname)) ;
  for s = 1:model.cnt_segment_features(j),
    weights = setfield(weights, [fieldname '_' names{s}], zeros(1,model.bins));
  end
end

%%% segment associated scores
for j=1:length(model.segment_feature_names);
  fieldname = sprintf('segment_score_%i',j);
  names = fieldnames(model.(fieldname)) ;
  for s = 1:model.cnt_segment_scores(j),
    weights = setfield(weights, [fieldname '_' names{s}], zeros(1,model.bins));
  end
end

%%% transitions
weights.transitions = zeros(1,model.cnt_transitions) ;
