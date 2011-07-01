function model = plif_to_param_mapping(model)

map_plif2param = init_weights(model);
map_plif2param = rmfield(map_plif2param,'transitions');
weights_name = fieldnames(map_plif2param);
start=0;
for w=1:length(weights_name)
  len=length(getfield(map_plif2param, weights_name{w}));
  model.plif2param{w} = start+(1:len);
  start = start + len;
end
model.transitions2param = start+1:model.cnt_parameters ;

assert(model.cnt_parameters-model.cnt_transitions==start)
