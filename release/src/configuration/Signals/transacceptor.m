function Signal = transacceptor(label_source)

% Signals = transacceptor(Signals,label_source)
% INPUT  - Signals is a struct with a field for each Signal. In this function only
%          the field transcc is modified.  
%        - label_source is a struct that should contain a field
%          transacc. this is used in the called function splicesite_template.
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: tansacceptor -> acceptor -> splicesite_template -> signal_template

INF = 1e6;  
signal_name = 'transacc';

Signal = acceptor(label_source,signal_name);

%%%% SET transacc STUFF

Signal.lwin = {[-190]} ;
Signal.rwin = {[80]};


Signal.filter_label.train.conf.anno_SL1 = 0;
Signal.filter_label.train.conf.anno_SL2 = [0 INF];
Signal.filter_label.train.conf.pred = [0 INF];

Signal.filter_label.eval.conf.anno_SL1 = 0;
Signal.filter_label.eval.conf.anno_SL2 = [0 INF];
Signal.filter_label.eval.conf.pred = [0 INF];


Signal.filter_label.test.conf.anno_SL1 = 0;
Signal.filter_label.test.conf.anno_SL2 = [0 INF];
Signal.filter_label.test.conf.pred = [0 INF];


Signal.filter_label.train.use_intergenic = 1;
Signal.filter_label.train.use_alt = 0;
Signal.filter_label.eval.use_intergenic = 1;
Signal.filter_label.eval.use_alt = 0;
Signal.filter_label.test.use_intergenic = 1;
Signal.filter_label.test.use_alt = 0;
 
Signal.filter_label.train.use_pred = 1;
Signal.filter_label.eval.use_pred = 1;
Signal.filter_label.test.use_pred = 1;

Signal.filter_label.train.only_acc = 1;
Signal.filter_label.eval.only_acc = 1;
Signal.filter_label.test.only_acc = 1;

Signal.filter_label.train.pos_dist = 20;
Signal.filter_label.eval.pos_dist = 20;
Signal.filter_label.test.pos_dist = 20;

Signal.export_settings.conf_cum_thresh = 0.05 ;
Signal.export_settings.resolution = 1 ;

Signal.Conf_names = {Signal.Conf_names{:},'isacc','has_tss','has_trans','dist'};


Signal.label_source_par.transacc.win_anno_SL1 = 50;
Signal.label_source_par.transacc.win_anno_SL2 = 50;
Signal.label_source_par.transacc.win_pred = 5;
Signal.label_source_par.transacc.threshold_pred = 0.5;

