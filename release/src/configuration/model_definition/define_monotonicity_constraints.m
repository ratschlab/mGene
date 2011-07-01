function model = define_monotonicity_constraints(model)


for nn = 1:length(model.track_names)
  fieldname = sprintf('track_%i_monoton', nn);
  model.(fieldname) = feval(model.track_monoton_functions{nn});
end

for nn = 1:length(model.segment_feature_names)
  fieldname = sprintf('segment_feature_%i_monoton', nn);
  model.(fieldname) = feval(model.segment_feature_monoton_functions{nn});
end
