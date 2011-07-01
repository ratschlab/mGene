function Signal = signal_template(label_source,signal_name)

% Signals = signal_template(Signal,label_source,signal_name)
% Signals is a struct with a field for each Signal. In this function only
% the field signal_name is added. Default values are set.


if nargin<2
  signal_name = 'dummy';
end

%%%%%

Signal.type = []; 
Signal.name = signal_name;
Signal.consensus = {[]};
Signal.offset = [];

Signal.lwin_big = [];
Signal.rwin_big = [];

Signal.label_fct = [];
Signal.filter_fct = 'filter_cands';
Signal.load_fct = 'load_data';
Signal.score_fct = [];

Signal.Conf_names = { 'altgenic','gene_id' };

Signal.sampling.downsample_neg = 1; % <=1 (for no downsampling: 1)  
Signal.sampling.mindist = 1; % >=1 (for all: 1)  
Signal.sampling.pos_neg_ratio = 1/100;

Signal.resolution   = 1; % for lsl

conf.USE_ALL = 0;
conf.conf_pos = struct;
conf.conf_pos_skip = struct;
conf.any_conf_pos = [];
conf.any_conf_pos_skip = [];
conf.conf_neg = struct;
conf.conf_neg_skip = struct;
conf.any_conf_neg = [];
conf.any_conf_neg_skip = [];
conf.use_intergenic=[];
conf.use_alt = [];
conf.use_label = 'label';

Signal.filter_label.train = conf ; 
Signal.filter_label.test = conf ;
Signal.filter_label.eval = conf ;

Signal.export_settings.conf_cum_thresh = 0 ;
Signal.export_settings.resolution = 1 ;

Signal.domain_adapt.name = 'target_only';
Signal.domain_adapt.Source_organisms = [];

Signal.method = struct;

label_source = label_source.(signal_name);
sources = fieldnames(label_source);
for i=1:length(sources)
  if label_source.(sources{i})
    Signal.Conf_names = {Signal.Conf_names{:},[sources{i} '_conf']};
    if isequal(signal_name,'acc')|isequal(signal_name,'don')
      Signal.Conf_names = {Signal.Conf_names{:},[sources{i} '_conf_skip']};
    end
  end
end



