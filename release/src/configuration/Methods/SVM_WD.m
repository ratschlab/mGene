function method = SVM_WD(Signal)

method = SVM(Signal);
method.feature_type = 'char';

kernel.name = 'WD';
kernel.order = Signal.order;
kernel.shift = 0 ;
kernel.lwin = Signal.lwin{1};
kernel.rwin = Signal.rwin{1};

method.kernels{1} = kernel; 

