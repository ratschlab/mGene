function [Data,par] = load_data(Signal,fn_filter_settings,fn_examples, which, filter_mode,subsample_neg,subsample_pos)
% [Data,par] = load_data(Signal,fn_filter_settings,fn_examples, which, filter_mode,subsample_neg,subsample_pos)
% [Data,par] = load_data(PAR, which, filter_mode,subsample_neg,subsample_pos)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load data of training set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method = PAR.method;
% Signal = PAR.Signals.(PAR.Signal_name);
%  FN = PAR.FN;
% input_sig = FN.input_sig.(Signal_name);
% fn_filter_settings = input_sig.fn_filter_settings;
% fn_examples = input_sig.fn_examples;

method = Signal.method;

XT = []; LT = [] ; POS.region_id = []; POS.pos=[]; POS_tmp = []; 
SC = []; use_sc = [];
if nargin<4
  subsample_neg=1;
end
if nargin<5
  subsample_pos=1;
end

use_cons = Signal.method.use_cons;
% if nargin<6
%   use_cons = 0;
% end
for w = which
  % create filename of input file
  filename = gen_example_filename(Signal,fn_filter_settings,fn_examples,filter_mode,w);
  
  fprintf('loading XT, LT and POS from file %s \n',filename)
  L = load(filename,'LT');
  L.XT = load_matrix(filename, 'XT');
  P = load(filename, 'POS');
  %if isfield(signal,'center')
  %  par.center = signal.center;
  %else
  %  par.center = size(L.XT,1)./2;
  %end
  
  have_pos = isfield(P, 'POS') ;
  have_lt = isfield(L, 'LT') ;
  if ~have_lt, L.LT=[] ; end ;
  have_sc = isfield(L, 'SC_') ;
  if ~have_sc, L.SC_=[] ; end ;  

  [L.XT, L.LT, P.POS, L.SC_] = subsample(L.XT, L.LT, P.POS, L.SC_, subsample_neg, subsample_pos) ;

  if ~have_lt, L=rmfield(L, 'LT') ; end ;
  if ~have_sc, L=rmfield(L, 'SC_') ; end ;
  
  if ~isfield(method,'feature_type') || isempty(method.feature_type) ||isequal(method.feature_type, 'char')
    L.XT_char = char(L.XT) ;
    for i=1:size(L.XT,1)
      L.XT_char(i,:) = char(double(L.XT(i,:))) ;
    end ;
    XT = [XT L.XT_char];
  elseif isequal(method.feature_type, 'uint8')
    L.XT_char = uint8(L.XT) ;
    for i=1:size(L.XT,1)
      L.XT_char(i,:) = uint8(double(L.XT(i,:))) ;
    end ;
    XT = [XT L.XT_char];
  elseif isequal(method.feature_type, 'double')
    XT = [XT L.XT];
  else
    error('unknown feature type; should be "char" or "double"')
  end
  if isfield(L,'LT')
    if size(L.LT,1)>size(L.LT,2)
      L.LT = L.LT';
    end
    LT = [LT L.LT] ;
  end
  
  if (isfield(L,'SC_')&& isfield(L.use_sc)) || use_cons
    fprintf('using consensus information\n');
    SC = [SC L.SC_] ;
    use_sc = [use_sc L.use_sc];
  end

  POS.region_id =[POS.region_id P.POS.region_id];
  POS.pos = [POS.pos P.POS.pos];
end ;
%%% just for compatibility with older versions 
  
% check whether everything is ok (matlab bug successfully avoided?)
for i=1:size(XT,1)
  if all(XT(i,:)==0)
    % i
    % keyboard
  end
  if  ~isfield(method,'feature_type') || isequal(method.feature_type, 'char')
    assert(all(XT(i,:)~=0)) ;
  end
end ;

clear L
if ~isempty(LT)
  fprintf('%i positive, %i negative examples (%.2f %% positives), sequences of length %i \n',sum(LT==1),sum(LT==-1),sum(LT==1)/length(LT)*100, size(XT,1));
end

if ~isempty(SC)
  XT(:,~use_sc) = [] ;
  POS(~use_sc) =[];
  LT(~use_sc) = [];
  SC(:,~use_sc) = [];
  if ~isempty(LT)
    fprintf('%i positive, %i negative examples (%.2f %% positives) \n',sum(LT==1),sum(LT==-1),sum(LT==1)/sum(LT==1|LT==-1)*100);
  end
end

Data.LT = LT;
% Data.XT = char(XT) ;
Data.XT = XT;
Data.Pos = POS;
clear XT LT POS
%this assert is very mem consuming
%assert(~any(XT(:)~='A' & XT(:)~='G' & XT(:)~='C' & XT(:)~='T')) ;
% XT = to_upper_case(char(XT)) ;
if ~isempty(SC)
  Data.SC = char(SC) ;
  assert(~any(SC(:)~='A' & SC(:)~='G' & SC(:)~='C' & SC(:)~='T')) ;
  
  % SC = to_upper_case(char(SC)) ;
end


assert(size(Data.XT,2)==length(Data.LT))
assert(length(Data.LT)==length(Data.Pos.pos))
assert(length(Data.LT)==length(Data.Pos.region_id))
