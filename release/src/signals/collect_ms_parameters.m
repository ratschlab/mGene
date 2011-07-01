function par_ms= collect_ms_parameters(kernels, par_ms)
%par_ms= collect_model_parameters(kernels, par_ms)

nofKernels = length(kernels);
assert(nofKernels>=1)
for k = 1:nofKernels
  kernel = kernels{k};
  for f =fieldnames(kernel)'
    if isnumeric(kernel.(f{1}))&~isempty(kernel.(f{1}))
      par_ms.(sprintf('kernel_%i_%s',k,f{1})) = kernel.(f{1});
    end
  end
end

