function Signal = splicesites_template(label_source,signal_name)

% Signals = splicesites_template(Signals,label_source,signal_name)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field signal_name is added or
%          modified. 
%          Default values are set.
%        - label_source is a struct that should contain a field (signal_name). 
%          this is used in the called function signal_template.
%        - signal_name: 'don', 'acc','transacc'
%  
% OUTPUT Signals with modified field signal_name; 
%  
% HIERACHY: acceptor/donor -> splicesite_template -> signal_template
%
% see also transcriptends_template/translat_template
  
  
INF=1e6;

Signal = signal_template(label_source,signal_name);

%%%% SET SPLICE SITE STUFF

Signal.type = 'splicesites'; 
Signal.label_fct = 'get_cand_splicesites';

Signal.lwin_big = 199;
Signal.rwin_big = 199;

Signal.filter_label.train.use_intergenic=0;
Signal.filter_label.train.use_alt=0;
Signal.filter_label.train.use_alt_neg=0;
Signal.filter_label.test.use_intergenic=0;
Signal.filter_label.test.use_alt=0;
Signal.filter_label.test.use_alt_neg=0;

Signal.filter_label.train.use_label = 'label';
Signal.filter_label.test.use_label = 'label';


Signal.filter_label.train.use_decoys = 1;
Signal.filter_label.test.use_decoys = 1;


Signal.filter_label.train.conf_pos.est_clusters = [0 INF];
Signal.filter_label.train.conf_pos.cDNA = [0 INF];
Signal.filter_label.train.conf_pos.fulllength = [0 INF];
Signal.filter_label.train.conf_pos.anno = [0:1] ;

Signal.filter_label.train.conf_pos_skip.est_clusters = [0 INF];
Signal.filter_label.train.conf_pos_skip.cDNA = [0 INF];
Signal.filter_label.train.conf_pos_skip.fulllength = [0 INF];
Signal.filter_label.train.conf_pos_skip.anno = [-inf] ;

Signal.filter_label.train.conf_neg.est_clusters = [0 INF];
Signal.filter_label.train.conf_neg.cDNA = [0 INF];
Signal.filter_label.train.conf_neg.fulllength = [0 INF];
Signal.filter_label.train.conf_neg.anno = [-1:1] ;

Signal.filter_label.train.conf_neg_skip.est_clusters = [0 INF];
Signal.filter_label.train.conf_neg_skip.cDNA = [0 INF];
Signal.filter_label.train.conf_neg_skip.fulllength = [0 INF];
Signal.filter_label.train.conf_neg_skip.anno = [-1:1] ;

Signal.filter_label.test = Signal.filter_label.train;
Signal.filter_label.eval = Signal.filter_label.train;

% Signal.save_confs = 1;
% Signal.save_alts = 1;


Signal.export_settings.conf_cum_thresh = 0.5 ;
Signal.export_settings.resolution = 1 ;


Signal.Conf_names = {Signal.Conf_names{:},'alt_exon','alt_intron'};




