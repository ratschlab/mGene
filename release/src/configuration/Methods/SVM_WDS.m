function method = SVM_WDS(Signal)

% Signals = SVM_WDS(Signals,signal_name)
% INPUT  - Signals is a struct with a field for each Signal. In this function only
%          the field (signal_name) is modified.
%          Particuarly a field method is added.
%        - signal_name can be something in {'acc','don','tss','cleave','polya','tis','cdsStop'};
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: SVM_WDS -> SVM -> method_template

  
method = SVM(Signal) ;
method.feature_type = 'char';


kernel.name = 'WDS';
kernel.order = Signal.order;
kernel.shift = Signal.shift;
kernel.lwin =  Signal.lwin{1};
kernel.rwin = Signal.rwin{1};

method.kernels{1} = kernel; 
