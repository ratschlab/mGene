function method = SVM_COMBINED(Signal)

% Signals = SVM_COMBINED(Signals,signal_name)
% INPUT  - Signals is a struct with a field for each Signal. In this function only
%          the field (signal_name) is modified.
%          Particuarly a field method is added.
%        - signal_name can be something in {'acc','don','tss','cleave','polya','tis','cdsStop'};
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: SVM_COMBINMED -> SVM -> method_template


method = SVM(Signal) ;

nofKernels = length(Signal.kernel_name);

assert(length(Signal.lwin)==nofKernels)
assert(length(Signal.rwin)==nofKernels)


for k=1:nofKernels
  for n = Signal.kernel_par_names 
    method.kernels{k}.(n{1}) = Signal.(n{1}){k};
  end
end
if isfield(Signal,'subkernel_weights')
  method.subkernel_weights = Signal.subkernel_weights;
end

