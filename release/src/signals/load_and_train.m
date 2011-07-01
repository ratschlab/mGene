function [x y] = load_and_train(P, Train_Data)

% load_and_train(P)
%
% INPUT: P with fields:   
%         - fn_set,             file name of sets to be loaded 
%         - SETs.which_svm,     specifies the sets on which training is
%                               performed (e.g [1 2 3])
%         - SETs.train_subsample_neg  number between 0-1: for subsampling of negative examples, 
%                                0: setting label of examples labeled with 0 to -1   
%                               -1: removing examples labeled with 0 
%                               default(1) no subsampling 
%         - SETs.train_subsample_pos  number between 0-1: for subsampling of positive examples,                             
%                               default(1) no subsampling 
%         - par,                window around consensus to be used  
%         - use_cons            specifies if conservation is to be used
%         - fn_SVMs             file nmame to save results  
%  
%   
% OUTPUT: Train_Data with fields:
%  (XT,LT,SC,pos,reg,str,P)
%  saved to file P.fn_SVMs
  
paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load data of training set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SETs = P.SETs;
Signal = P.Signal;
fn_filter_settings = P.fn_filter_settings;
fn_examples = P.fn_examples;

which_train = SETs.which_train;
train_fct = Signal.method.train_fct;

if ~isfield(SETs,'train_subsample_pos') && isfield(SETs,'subsample_pos')
  subsample_pos = SETs.subsample_pos ;
  error('warning: please rename field from subsample_pos to train_subsample_pos') ;
else
  subsample_pos = SETs.train_subsample_pos ;
end

if ~isfield(SETs,'train_subsample_neg') && isfield(SETs,'subsample_neg')
  subsample_neg = subsample_neg ;
  error('warning: please rename field from subsample_neg to train_subsample_neg') ;
else
  subsample_neg = SETs.train_subsample_neg ;
end

if ~exist('Train_Data', 'var')
  Train_Data = feval(Signal.load_fct,Signal,fn_filter_settings,fn_examples, which_train, 'train',subsample_neg, subsample_pos);
end
Train_Data.P = P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRAIN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = cputime;
feval(train_fct,Train_Data);
fprintf('\ntraining took %.2f min\n',(cputime-t)/60);
% if PAR.method.name=='SVM'
%   tssTrain(Train_Data)
% elseif PARn.method.name=='PWM'
%   hmm_train(Train_Data)
% else 
%   error('unknown method');
% end

% dummy return argument for rproc
x=0;
y=0;
