function weight_names=init_weight_names(model) ;
%  weight_names=weight_names(model) ;

%%% content detectors
content_names = fieldnames(model.contents) ;
for s = 1:model.cnt_contents,
  weight_names{model.contents.(content_names{s})} = ['contents_' content_names{s}] ;
end

%%% length detectors
length_names = fieldnames(model.lengths) ;
for s = 1:model.cnt_lengths,
  weight_names{model.lengths.(length_names{s})} = ['lengths_' length_names{s}] ;
end

%%% signal detectors
signal_names = fieldnames(model.signals) ;
for s = 1:model.cnt_signals,
  weight_names{model.signals.(signal_names{s})} = ['signals_' signal_names{s}] ;
end

%% positons associated features
for j=1:length(model.track_names)
  fieldname = sprintf('track_%i',j);
  names = fieldnames(model.(fieldname)) ;
  for s = 1:model.cnt_tracks,
    weight_names{model.(fieldname).(names{s})} = [ fieldname '_' names{s}] ;
  end
end
%%% segment associated features
for j=1:length(model.segment_feature_names)
  fieldname = sprintf('segment_feature_%i',j);
  names = fieldnames(model.(fieldname)) ;
  for s = 1:model.cnt_segment_features,
    weight_names{model.(fieldname).(names{s})} = [ fieldname '_' names{s}] ;
  end
end

%%% segment associated features
for j=1:length(model.segment_feature_names)
  fieldname = sprintf('segment_score_%i',j);
  names = fieldnames(model.(fieldname)) ;
  for s = 1:model.cnt_segment_scores,
    weight_names{model.(fieldname).(names{s})} = [ fieldname '_' names{s}] ;
  end
end

%%% transitions
weight_names{end+1}='transitions' ;
