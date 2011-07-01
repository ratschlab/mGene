function [Data,par] = load_content_data(Content,fn_filter_settings,fn_examples, which, filter_mode,subsample_neg,subsample_pos)
% [Data,par] = load_content_data(Content,fn_filter_settings,fn_examples, which, filter_mode,subsample_neg,subsample_pos)


method = Content.method;

XT = []; LT = [] ;  POS.region_id = []; POS.pos=[]; POS_tmp = []; 

if nargin<4
  subsample_neg=1;
end
if nargin<5
  subsample_pos=1;
end

use_cons = method.use_cons;
% if nargin<6
%   use_cons = 0;
% end
for w = which
  % create filename of input file
  filename = gen_example_filename(Content,fn_filter_settings,fn_examples,filter_mode,w);
  
  fprintf('loading XT, LT and POS from file %s \n',filename)
  L = load(filename,'LT','XT','POS');
  %  L.XT = load_matrix(filename, 'XT');
  % [L.XT, L.LT, P.POS, L.SC_] = subsample(L.XT, L.LT, P.POS, L.SC_, subsample_neg, subsample_pos) ;

  if isfield(L,'LT')
    if size(L.LT,1)>1,
      L.LT = L.LT';
    end
    LT = [LT L.LT] ;
  end
  POS.region_id =[POS.region_id L.POS.region_id];
  POS.pos = [POS.pos L.POS.pos];
  
  % generate feature vectors
  %----------------------------------------------------------------
  dim = 0;
  for o=1:length(method.kernels)
    dim = dim + 4^method.kernels{o}.wordlen;
  end
  xt = zeros(dim,length(L.XT));
  for j=1:length(L.XT) 
    dim = 0;
    for o=1:length(method.kernels)
      [mask hist] = compute_mers(L.XT{j}, method.kernels{o}.wordlen, method.kernels{o}.stepping, method.kernels{o}.offset);
      xt(dim+[1:4^method.kernels{o}.wordlen],j) = hist';
      dim = dim+ 4^method.kernels{o}.wordlen;
    end
  end
  clear L
  XT=[XT xt];
  
end ;

if ~isempty(LT)
  fprintf('%i positive, %i negative examples (%.2f %% positives), sequences of length %i \n',sum(LT==1),sum(LT==-1),sum(LT==1)/length(LT)*100, size(XT,1));
end


Data.LT = LT;
% Data.XT = char(XT) ;
Data.XT = XT;
Data.Pos = POS;
clear XT LT POS

assert(size(Data.XT,2)==length(Data.LT))
% assert(length(Data.LT)==length(Data.Pos.pos))
% assert(length(Data.LT)==length(Data.Pos.region_id))






