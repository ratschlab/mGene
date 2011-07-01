function method = SVM(Signal)
  
method = method_template;

% set parameters, model selection 

method.name = 'SVM';
method.train_fct = 'svm_train';
method.predict_fct = 'svm_predict';
method.kernels = {[]};
method.poim_maxOrder = 8;
method.svm_name = 'LIGHT';
method.qpsize = 50;
method.threads = 16;
method.center = Signal.lwin_big;
method.par_ms.C = Signal.C;

%PAR.method.svm_name = 'GPBTSVM';PAR.method.qpsize=500;

% use multi kernel learning
method.use_mkl = 0;




