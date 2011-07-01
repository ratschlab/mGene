function Signal = cdsStop(label_source)
  
% Signals = cdsStop(Signals,label_source)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field cdsStop is modified.
%        - label_source is a struct that should contain a field cdsStop. 
%          this is used in the called function signal_template.
%       
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: tis/cdsStop -> translatends_template -> signal_template
   
signal_name = 'cdsStop';

Signal = translatends_template(label_source,signal_name);

%%%% SET cdsStop STUFF

Signal.consensus = {'TAA','TAG','TGA'};


Signal.kernel_name= {'SP','WDS','SP'};
Signal.kernel_par_names={'kernel_name','order','shift','shift_const','wordlen','lwin','rwin'};
Signal.lwin = {[-200] [-110] [0]} ; 
Signal.rwin = {[0] [10] [180]};
Signal.shift={[],[0],[]};
Signal.shift_const={[],[0],[]};
Signal.order={[],[20],[]};  
Signal.wordlen = {[2] [] [4]}; 
Signal.C = 2.8284; % logC=1.5


