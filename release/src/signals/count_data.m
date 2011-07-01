function [num_data, num_pos, num_neg] = count_data(Signal, fn_filter_settings, fn_examples, which, filter_mode, subsample_neg, subsample_pos)
% [num_data, num_pos, num_neg] = count_data(Signal, fn_filter_settings, fn_examples, which, filter_mode,subsample_neg,subsample_pos)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% count data of training set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = Signal.method;

LT = [] ; POS.pos=[]; 
if nargin<4
  subsample_neg=1;
end
if nargin<5
  subsample_pos=1;
end

for w = which
  % create filename of input file
  filename = gen_example_filename(Signal,fn_filter_settings,fn_examples,filter_mode,w);
  
  fprintf('loading LT and POS from file %s \n',filename)
  L = load(filename, 'LT');
  P = load(filename, 'POS');
  
  have_pos = isfield(P, 'POS') ;
  have_lt = isfield(L, 'LT') ;
  if ~have_lt, L.LT=[] ; end ;

  [tmp, L.LT, P.POS] = subsample(L.LT, L.LT, P.POS, [], subsample_neg, subsample_pos) ;

  if ~have_lt, L=rmfield(L, 'LT') ; end ;
  
  if isfield(L,'LT')
    if size(L.LT,1)>size(L.LT,2)
      L.LT = L.LT';
    end
    LT = [LT L.LT] ;
  end
  
  POS.pos = [POS.pos P.POS.pos];
end ;
%%% just for compatibility with older versions 
  
clear L
if ~isempty(LT)
  fprintf('%i positive, %i negative examples (%.2f %% positives)\n',sum(LT==1),sum(LT==-1),sum(LT==1)/length(LT)*100);
  num_data = length(LT) ;
  num_pos = sum(LT==1) ;
  num_neg = sum(LT==-1) ;
else
  num_data = length(POS.pos) ;
  num_pos = nan ;
  num_neg = nan ;
end
