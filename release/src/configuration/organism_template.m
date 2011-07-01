function [label_source] = organism_template(label_source_orig,signal_names)


% PAR = set_label_source(PAR)
% 
% sets defaults for all label_sources (don't use anno, use EST confirmed only.)
% values that are already assigned will NOT be overriden.
  
  
%% default values for general use of annotation and EST confrimed sequences
label_source = label_source_orig;

if isempty(label_source.use_anno)
  label_source.use_anno = 0;
end
if isempty(label_source.use_confirmed)
  label_source.use_confirmed = 1;
end


%% default values for individual signals


for j=1:length(signal_names)
  signal_name = signal_names{j};
  if label_source.use_anno
    label_source.(signal_name).anno = 1;
  else
    label_source.(signal_name).anno = 0;
  end
  if label_source.use_confirmed & ~isequal(signal_name,'transacc')
    label_source.(signal_name).est_clusters = 1;
    label_source.(signal_name).cDNA = 1;
    label_source.(signal_name).fulllength = 1;
  elseif ~label_source.use_confirmed & ~isequal(signal_name,'transacc')
    label_source.(signal_name).est_clusters = 0;
    label_source.(signal_name).cDNA = 0;
    label_source.(signal_name).fulllength = 0;
  end
  label_source.(signal_name).db = 0;
  label_source.(signal_name).from_fn_candsites = 1 ;
end

label_source.tis.maxorf = 0;
label_source.cdsStop.maxorf = 0;

label_source.transacc.anno_SL1 = 0;
label_source.transacc.anno_SL2 = 0;
label_source.transacc.pred = 1;
label_source.transacc.db = 0;


