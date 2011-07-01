function vector = weights2vector(weights, weight_names)
% vector = weights2vector(weights)
  
vector = [];

%weights_name = fieldnames(weights);
for w=1:length(weight_names)
  vector=[vector, getfield(weights, weight_names{w})];
end
