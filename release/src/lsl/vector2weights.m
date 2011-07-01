function weights = vector2weights(vector, weights, weight_names)
% weights = vector2weights(vector, weights)
  
%weights_name = fieldnames(weights);
last_pos = 0 ;
for w = 1:length(weight_names)
  old_w = getfield(weights, weight_names{w}) ;
  weights = setfield(weights, weight_names{w}, vector(last_pos+1:last_pos+length(old_w)));
  last_pos = last_pos + length(old_w) ;
end
assert(last_pos==length(vector)) ;