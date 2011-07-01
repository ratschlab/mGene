function Signal = polya(label_source,fn_polya)
  
% Signals = polya(Signals,label_source)
% INPUT  - Signals is a struct with a field for each Signal. 
%          In this function only the field polya is modified.
%        - label_source is a struct that should contain a field polya. 
%          this is used in the called function signal_template.
%       
% OUTPUT Signals with modified field signal_name; 
%
% HIERACHY: polya -> signal_template

signal_name = 'polya' ; 
Signal = signal_template(label_source,signal_name);

%%%% SET POLYA STUFF

Signal.name = signal_name;
Signal.type = signal_name; 
Signal.label_fct = 'get_cand_others';
Signal.filter_fct = 'filter_cands';

Signal.offset = 0;

Signal.lwin_big = 310;
Signal.rwin_big = 310;

Signal.kernel_name= {'WDS','SP','SP','SP','SP','SP'};
Signal.kernel_par_names={'kernel_name','order','shift','shift_const','wordlen','lwin','rwin'};
Signal.subkernel_weights = [1 0.2 0.2 0.2 0.2 0.2];

Signal.lwin = {-300 -300 -149 1 151 -300} ;
Signal.rwin = {306 -150 0 150 306 306};
Signal.shift={[0 0.1 0.4],[],[],[],[],[]};
%Signal.shift={[0.1],[],[],[],[],[]};
Signal.shift_const={[0 10 20],[],[],[],[],[]};
%Signal.shift_const={[0],[],[],[],[],[]};
Signal.wordlen={[],[4],[4],[4],[4],[4]};
Signal.order={[18],[],[],[],[],[]};
Signal.C = [0.01 0.5 10];

Signal.consensus = {''};
if exist('fn_polya')
  Signal.consensus_fn = fn_polya;
  if fexist(Signal.consensus_fn)
    load(Signal.consensus_fn,'consensus')
    Signal.consensus = consensus;
  end
end

Signal.filter_label.train.use_intergenic=0;
Signal.filter_label.train.use_alt=0;
%Signal.filter_label.eval.use_intergenic=0;
%Signal.filter_label.eval.use_alt=0;
Signal.filter_label.test.use_intergenic=0;
Signal.filter_label.test.use_alt=0;

Signal.filter_label.train.conf_pos.anno = 0;
Signal.filter_label.train.conf_pos.from_cleave = [201 301 401];
Signal.filter_label.train.conf_pos.db = 0;
Signal.filter_label.train.conf_pos.fulllength = 0;
Signal.filter_label.train.onlywithpos = 0; 

%Signal.filter_label.eval.conf_pos.anno = 0;
%Signal.filter_label.eval.conf_pos.from_cleave = [201 301 401];
%Signal.filter_label.eval.conf_pos.db = 0;
%Signal.filter_label.eval.conf_pos.fulllength = 0;
%Signal.filter_label.eval.onlywithpos = 0; 

Signal.filter_label.test.conf_pos.anno = 0;
Signal.filter_label.test.conf_pos.from_cleave = [201 301 401];
Signal.filter_label.test.conf_pos.db = 0;
Signal.filter_label.test.conf_pos.fulllength = 0;
Signal.filter_label.test.onlywithpos = 0; 


Signal.Conf_names = {Signal.Conf_names{:},['has_' signal_name]};



Signal.label_source_par.polya.wl = 1;
Signal.label_source_par.polya.POSconstr =1;
Signal.label_source_par.polya.maxpcCov = 95;
Signal.label_source_par.polya.logdir = [];
