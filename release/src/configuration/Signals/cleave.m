function Signal = cleave(label_source)

% Signals = cleave(Signals,signal_name,label_source)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field cleave is modified.
%        - label_source is a struct that should contain a field cleave. 
%          this is used in the called function signal_template.
%
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: tss/cleave -> transcriptends_template -> signal_template

  
signal_name = 'cleave';

Signal = transcriptends_template(label_source,signal_name);

%%%% SET CLEAVE STUFF

Signal.lwin_big = 410;
Signal.rwin_big = 810;

Signal.kernel_name= {'SP','WDS','SP'};
Signal.kernel_par_names={'kernel_name','order','shift','shift_const','wordlen','lwin','rwin'};

Signal.lwin = {[-400], [-60], [0]};
Signal.rwin = {[100], [60], [800]};
Signal.wordlen = {[5], [], [3]}; 

Signal.order = {[],[18],[]};
Signal.shift = {[],[0.1],[]};
Signal.shift_const = {[],[0],[]};
Signal.C = 1; % log2(C)=0


