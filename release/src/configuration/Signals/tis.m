function Signal = tis(label_source)
  
% Signals = tis(Signals,label_source)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field tis is modified.
%        - label_source is a struct that should contain a field tis. 
%          this is used in the called function signal_template.
%       
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: tis/cdsStop -> translatends_template -> signal_template

    
signal_name = 'tis';

Signal = translatends_template(label_source,signal_name);

%%%% SET TIS STUFF

Signal.consensus = {'ATG'};


Signal.kernel_name= {'SP','WDS','SP'};
Signal.kernel_par_names={'kernel_name','order','shift','shift_const','wordlen','lwin','rwin'};
Signal.lwin = {[-200] [-30] [0]} ; 
Signal.rwin = {[0] [110] [180]}; 
Signal.shift={[],[0],[]};
Signal.shift_const={[],[0],[]};
Signal.order={[],[24],[]};  
Signal.wordlen = {[3] [] [5]}; 
Signal.C = 1; % log2(C)=1.0









