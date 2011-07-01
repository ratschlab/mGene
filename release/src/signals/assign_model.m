function method = assign_model(method, models, par_names,  model_number)

%function method = assign_model(method, models, par_names, model_number)  

  
for num_param = 1:size(models,1)
  method.par_ms.(par_names{num_param}) = models(num_param,model_number) ;
end

if isfield(method, 'kernels') && ~isempty(method.kernels{1})
  nofKernels = length(method.kernels);
  for k =1:nofKernels
    kernel = method.kernels{k};
    for f =fieldnames(kernel)'
      if isnumeric(kernel.(f{1}))&~isempty(kernel.(f{1}))
        kernel.(f{1})= method.par_ms.(sprintf('kernel_%i_%s',k,f{1}));
      end
    end
    method.kernels{k}=kernel;
  end
end
