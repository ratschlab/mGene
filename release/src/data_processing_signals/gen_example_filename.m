function [fn, exists]=gen_example_filename(Signal,fn_filter_settings,fn_examples,which,split)

% [fn, exists]=gen_example_filename(Signal,fn_filter_settings,fn_examples,which,split)
% [fn, exists]=gen_example_filename(PAR,which,split)
% which is 'train','eval','test'

% sig = PAR.Signal_name;
% Signal = PAR.Signals.(sig); 
% fn_filter_settings =  PAR.FN.input_sig.(sig).input_sig;
% fn_examples = PAR.FN.input_sig.(sig).fn_examples;

whichs = fieldnames(Signal.filter_label);


if isnumeric(which)
  assert(which<=length(whichs))
  which = whichs{which};
else
  assert(ismember(which,whichs))
end
clear whichs

current_settings = Signal.filter_label.(which);

if fexist(fn_filter_settings)
  load(fn_filter_settings, 'filter_settings')
  idx = find([filter_settings.split]==split);
else
  filter_settings = [];
  idx = [];
end
idx_setting = [];is_equal=0;
for i=idx
  temp = rmfield(filter_settings(i),'split');
  temp = rmfield(temp,'filename');
  if isequal(temp, current_settings)
    idx_setting = [idx_setting i];
    is_equal=1;
  else
    assert(length(fieldnames(temp))==length(fieldnames(current_settings))&&...
           all(strcmp(sort(fieldnames(temp)),sort(fieldnames(current_settings)))))
  end
end
assert(length(idx_setting)<=1);
if isempty(idx_setting)
  temp = struct;
end
assert(length(setdiff(fieldnames(temp),fieldnames(current_settings)))==0)

if ~isempty(idx_setting)
  fn = filter_settings(idx_setting).filename;
  if fexist(fn)
    exists = 1;
  else
    exists = 0;
  end
else 
  settings_names = fieldnames(current_settings);
  s = length(filter_settings)+1;
  for j = 1:length(settings_names)
    filter_settings(s).(settings_names{j}) = current_settings.(settings_names{j});
  end
  filter_settings(s).split = split;
  filter_settings(s).filename = sprintf('%s_settings=%i_split=%i.mat',fn_examples, s, split);
  fn = filter_settings(s).filename;
  if fexist(fn)
    error('I do not know where that file comes from!!!')
  else
    exists = 0;
  end 
  save(fn_filter_settings, 'filter_settings')
end



       
%fn = sprintf('%s_%s_pc=%i_pa=%i_nc=%i_na=%i_split=%i.mat',...
%	input_sig.fn_examples,PAR.Signal_name,pos_conf,...
%	pos_anno,neg_conf,neg_anno, s);
