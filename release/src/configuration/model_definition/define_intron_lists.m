function model = define_intron_lists(model)
  


for nn = 1:length(model.segment_feature_names)
  name = sprintf('segment_feature_%i', nn);
  start = model.cnt_plifs+1 ;
  cnt=0;

  model.(name).intron     = start+cnt ;cnt=cnt+1;
  
  %cnt_name = ['cnt_' name]; 
  model.cnt_segment_features(nn) = length(fieldnames(model.(name)));
  plif_ids_name = [name '_plif_ids'];
  model.(plif_ids_name) = [start:start+model.cnt_segment_features(nn)-1];

  % declare range
  range_name = [name '_range'];
  model.(range_name) = struct;
  field_names = fieldnames(model.(name)) ;
  for s = 1:length(field_names)
    model.(range_name) = setfield(model.(range_name), field_names{s},[0,19]);
  end
  fields_name = [name '_fields'];
  model.(fields_name) = fieldnames(model.(name));
  model.cnt_plifs = model.cnt_plifs + model.cnt_segment_features(nn);

  %% QUALITY track
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  name = sprintf('segment_score_%i',nn);
  start = model.cnt_plifs+1 ;
  cnt=0;

  model.(name).intron     = start+cnt ;cnt=cnt+1;
  
  %cnt_name = ['cnt_' name]; 
  model.cnt_segment_scores(nn) = length(fieldnames(model.(name)));
  plif_ids_name = [name '_plif_ids'];
  model.(plif_ids_name) = [start:start+model.cnt_segment_scores(nn)-1];

  % declare range
  range_name = [name '_range'];
  model.(range_name) = struct;
  field_names = fieldnames(model.(name)) ;
  for s = 1:length(field_names)
    model.(range_name) = setfield(model.(range_name), field_names{s},[0,19]);
  end
  fields_name = [name '_fields'];
  model.(fields_name) = fieldnames(model.(name));
  model.cnt_plifs = model.cnt_plifs + model.cnt_segment_scores(nn);

end
