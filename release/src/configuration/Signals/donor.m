function Signal = donor(label_source)

% Signals = donor(Signals,signal_name,label_source)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field don is modified.
%        - label_source is a struct that should contain a field don. 
%          this is used in the called function signal_template.
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: acceptor/donor -> splicesite_template -> signal_template

 
signal_name = 'don';

Signal = splicesites_template(label_source,signal_name);

%%%% SET DONOR STUFF

Signal.consensus = {'GC' 'GT'};
Signal.offset = 0;

Signal.lwin = {[-80]};
Signal.rwin = {[60]};

Signal.order = [22];
Signal.shift = 0.5;
Signal.C = [3]; 

Signal.Conf_names = {Signal.Conf_names{:},'alt_5prime'};



