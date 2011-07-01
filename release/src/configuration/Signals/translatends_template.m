function Signal = translatends_template(label_source,signal_name)

% Signals = translatends_template(Signals,label_source,signal_name)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field signal_name is added or
%          modified. 
%          Default values are set.
%        - label_source is a struct that should contain a field (signal_name). 
%          this is used in the called function signal_template.
%        - signal_name:'tis', 'cdsStop'
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: tis/cdsStop -> translatends_template -> signal_template
%
% see also splicesite_template/transcriptends_template
  
INF=1e6;
  
Signal = signal_template(label_source,signal_name);

%%%% SET TRANSLAT END STUFF

Signal.type = 'translatends' ;
Signal.label_fct = 'get_cand_others';

Signal.offset = 0;
Signal.lwin_big = 210;
Signal.rwin_big = 200;



Signal.filter_label.train.conf_pos.anno = [0 1] ;
Signal.filter_label.train.conf_pos.maxorf = INF;
Signal.filter_label.test.conf_pos.anno = [0 1] ;
Signal.filter_label.test.conf_pos.maxorf =INF;

Signal.filter_label.train.use_intergenic=0;
Signal.filter_label.train.use_alt=0;
Signal.filter_label.test.use_intergenic=0;
Signal.filter_label.test.use_alt=0;

Signal.filter_label.eval = Signal.filter_label.test;

% PAR.label_source.cds.par.max_nof_paths = 100;
% PAR.label_source.cds.par.nof_best_kept = 2;

