function method = method_template

% Signals = method_template(Signal,signal_name)
% Signals is a struct with a field for each Signal. In this function only
% the field Signal.(signal_name).method is added. Default values are set.

method.name = [];
method.train_fct = [];
method.predict_fct = [];

% PAR.method.svm_name = 'GPBTSVM';PAR.method.qpsize=500;
% use conservation information
% method.use_cons = 0;
% load training examples only once and calculate all models in one job
% method.submit_batch=0;

method.plif.bins = 50;
method.plif.min_incr = 1e-6;

method.use_cons = 0;
method.feature_type = [];
%assert(max(abs([PAR.Signals.(PAR.Signal_name).lwin{:}]))<PAR.Signals.(PAR.Signal_name).lwin_big)

%assert(max(abs([PAR.Signals.(PAR.Signal_name).lwin{:}]))<PAR.Signals.(PAR.Signal_name).rwin_big)



