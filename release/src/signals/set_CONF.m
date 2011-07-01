function CONF = set_CONF(Signal,ALL,sig,info_names,signal_name)
% CONF = set_CONF(Signal,ALL,sig,info_names)

%%% GENERAL 
idx = find(strcmp(Signal.Conf_names,'altgenic'));
if ~isempty(idx)
  CONF(idx,:) = ALL.ALTGENIC;
end
idx = find(strcmp(Signal.Conf_names,'gene_id'));
if ~isempty(idx)
  CONF(idx,:) = ALL.GENE_ID;
end

for i=1:length(info_names)
  idx = find(strcmp(Signal.Conf_names,[info_names{i} '_conf']));
  CONF(idx,:) = ALL.COV1(i,:);
  if ~isequal(sig,'other')&~isequal(Signal.name,'transacc')
    idx = find(strcmp(Signal.Conf_names,[info_names{i} '_conf_skip']));
    CONF(idx,:) = ALL.COV2(i,:);
  end
end

%%% ALT SPLICE
idx = find(strcmp(Signal.Conf_names,'alt_exon'));
if ~isempty(idx)
  CONF(idx,:) = ALL.LABEL_2(1,:);
end
idx = find(strcmp(Signal.Conf_names,'alt_intron'));
if ~isempty(idx)
  CONF(idx,:) = ALL.LABEL_2(2,:);
end
idx = find(strcmp(Signal.Conf_names,'alt_3prime'));
if ~isempty(idx)
  CONF(idx,:) = ALL.LABEL_2(3,:);
end
idx = find(strcmp(Signal.Conf_names,'alt_5prime'));
if ~isempty(idx)
  CONF(idx,:) = ALL.LABEL_2(3,:);
end

idx = find(strcmp(Signal.Conf_names,['is' sig]));
if ~isempty(idx)
  CONF(idx,:) = ALL.ISSPLICE;
end

%%% SPECIAL THINGS
idx = find(strcmp(Signal.Conf_names,['has_' signal_name]));
if ~isempty(idx)
  CONF(idx,:) = ALL.HAS_SIG;
end
idx = find(strcmp(Signal.Conf_names,'has_trans'));
if ~isempty(idx)
  CONF(idx,:) = ALL.TRANS;
end


idx = find(strcmp(Signal.Conf_names,'dist'));
if ~isempty(idx)
  CONF(idx,:) = ALL.DIST;
end

if exist('CONF', 'var') && ~(size(CONF,1)==length(Signal.Conf_names))
  warning('this was an assert') ;
end ;
