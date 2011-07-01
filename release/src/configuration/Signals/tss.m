function Signal = tss(label_source,has_transacc)
  
% Signals = tss(Signals,label_source,has_transacc)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field tss is modified.
%        - label_source is a struct that should contain a field tss. 
%          this is used in the called function signal_template.
%        - has_transacc: is per default set to zero, only 1 for nematodes. 
%       
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: tss/cleave -> transcriptends_template -> signal_template

  
if nargin<3
  has_transacc=0; 
end
   
signal_name = 'tss';

Signal = transcriptends_template(label_source,signal_name);

%%%% SET TSS STUFF


Signal.lwin_big = 601;
Signal.rwin_big = 901;

Signal.kernel_name= {'SP','WDS','SP'};
Signal.kernel_par_names={'kernel_name','order','shift','shift_const','wordlen','lwin','rwin'};

Signal.lwin = { [-600] [-70] [   0] };
Signal.rwin = { [+100] [+70] [+900] };
Signal.wordlen = { [4] [] [4] };

Signal.order ={[] [16] []};
Signal.shift ={[] [0] []};
Signal.shift_const ={[] [28] []};
Signal.C = 2; % log2(C)=1.0

if has_transacc
  Signal.Conf_names = {Signal.Conf_names{:},'has_trans'};
  Signal.filter_label.train.rm_trans = 1;
  Signal.filter_label.eval.rm_trans = 1;
  Signal.filter_label.test.rm_trans = 1;
else
  Signal.filter_label.train.rm_trans = 0;
  Signal.filter_label.eval.rm_trans = 0;
  Signal.filter_label.test.rm_trans = 0;
end

%%%%%% for human
%Signal.sampling.downsample_neg = 1/160; % <=1 (for no downsampling: 1)


