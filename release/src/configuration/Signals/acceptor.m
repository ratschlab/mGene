function Signal = acceptor(label_source,signal_name)

% Signal = acceptor(label_source,signal_name)
% INPUT  - Signals is a struct with a field for each Signal. In this function only
%          the field acc (transcc) is modified.  
%        - label_source is a struct that should contain a field
%          (signal_name). this is used in the called function splicesite_template.
%        - name can either be 'acc' or 'transacc'  
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: acceptor/donor -> splicesite_template -> signal_template
% (function acceptor is also called from transacc)
  
  
if nargin<3 
  signal_name = 'acc';
end

Signal = splicesites_template(label_source,signal_name);

%%%% SET ACCEPTOR STUFF

Signal.consensus = {'AG'};
Signal.offset = 2;

Signal.lwin = {[-60]} ;
Signal.rwin = {[80]} ;

Signal.order = [ 22];
Signal.shift = 0.5;
Signal.C = [3]; 

Signal.Conf_names = {Signal.Conf_names{:},'alt_3prime'};

Signals.(signal_name) = Signal;


