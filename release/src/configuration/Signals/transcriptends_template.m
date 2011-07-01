function Signal = transcriptends_template(label_source,signal_name)

% Signals = transcriptends_template(Signals,label_source,signal_name)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field signal_name is added or
%          modified. 
%          Default values are set.
%        - label_source is a struct that should contain a field (signal_name). 
%          this is used in the called function signal_template.
%        - signal_name:'tss', 'cleave'
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: tss/cleave -> transcriptends_template -> signal_template
%
% see also splicesite_template/translat_template
  
  
  
Signal = signal_template(label_source,signal_name);

%%%% SET TRANSCRIPT END STUFF

Signal.type = 'transcriptends'; 
Signal.label_fct = 'get_cand_others';


Signal.consensus = [];
Signal.offset = 0;

Signal.sampling.downsample_neg = 1/10; %% take randomly every 10th porsitions 
Signal.sampling.mindist = 10 ;


Signal.resolution   = 20; % for lsl

Signal.filter_label.train.conf_pos.anno = [0 1];
Signal.filter_label.train.conf_pos.est_clusters = [2:3];
Signal.filter_label.train.conf_pos.db = 0;
Signal.filter_label.train.conf_pos.fulllength = 0;
Signal.filter_label.train.onlywithpos=0; 

Signal.filter_label.train.use_intergenic=1;
Signal.filter_label.train.use_alt=0;

Signal.filter_label.test = Signal.filter_label.train;
Signal.filter_label.eval = Signal.filter_label.train;

Signal.export_settings.conf_cum_thresh = 0.01 ;
Signal.export_settings.resolution = 10 ;

Signal.Conf_names = {Signal.Conf_names{:},['has_' signal_name]};


