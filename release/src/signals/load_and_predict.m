function  [dummy1, dummy2]=load_and_predict(P, Test_Data)
% load_and_predict(P, Test_Data)

% used for both evaluation on the leftout set of the training region and
% prediction on test regions
% 
% * IF PAR.eval_on is set to 1 prediction is performed on the evaluation
%   set specified by PAR.SETs.which_val. The model saved in PAR.fn_SVMs is used for
%   prediction. Additionally RSV and RFV are determined.
%  
% * IF PAR.eval_on is set to 0  only prediction is performed on the test set
%   specified by PAR.SETs.which_pred.
% 
% 

paths

  dummy1 = [] ;
  dummy2 = [] ;

fn_SVMs = P.fn_SVMs ;
SETs = P.SETs;
pred_from_seq = P.pred_from_seq;
Signal = P.Signal;
fn_examples = P.fn_examples;
fn_filter_settings = P.fn_filter_settings;

if isfield(P,'avg_all_svms') 
  avg_all_svms = P.avg_all_svms;
else
  avg_all_svms = 0;
end

if isfield(P,'check_overlap')
  check_overlap = P.check_overlap;
else
  check_overlap = 1;
end

if isfield(P,'eval_on') && P.eval_on
  which_set = SETs.which_val;
  if isfield(P,'fn_eval')
    fn_pred = P.fn_eval;
  else
    fn_pred = [];
  end
elseif isfield(P,'plif_on') && P.plif_on
  which_set = SETs.which_plif;
  if isfield(P,'fn_plifs')
    fn_pred = P.fn_plifs;
  else
    fn_pred = [];
  end
else
  which_set = SETs.which_pred;
  if isfield(P,'fn_pred')
    fn_pred = P.fn_pred;
  else
    fn_pred = [];
  end
end

if pred_from_seq
  fn_pos = P.fn_pos;
  if isequal(fn_pos(1),'~')
    fn_pos = sprintf('%s/%s', getenv('HOME'), fn_pos(2:end));
  end
  region = P.region;
  fn_pred = P.fn_pred;
  chunk_idx = 1;% P.chunk_idx;
  RPROC = struct;
end

load(fn_SVMs, 'Train') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(SETs,'test_subsample_pos')
  subsample_pos = 1;
else
  subsample_pos = SETs.test_subsample_pos;
end
if ~isfield(SETs,'test_subsample_neg')
  subsample_neg = 1;
else
  subsample_neg= SETs.test_subsample_neg;
end


if ~pred_from_seq
  pred_from_seq=0;
  if ~exist('Test_Data', 'var')
    Test_Data = feval(Signal.load_fct,Signal,fn_filter_settings,fn_examples, which_set, 'test',subsample_neg, subsample_pos);
  end
  
  %  if isequal(PAR.method.name,'SVM')
    % [RFV RSV]=tssTest(Train, Test_Data, fn_SVMs, fn_pred);
    % varargout{1} = RFV;
    % varargout{2} = RSV;
 %  elseif isequal(PAR.method.name,'PWM')
    % hmm_test(Train, Test_Data, fn_SVMs,fn_pred);
 %  end
else %genomweite Vorhersage
  
  region = load_sequence(region);
  region = retrieve_signals(region,fn_pos,Signal.name,{'svm'});
  Test_Data.seq = region.seq;
  Test_Data.pos = region.Signals.(Signal.name).pos';
  Test_Data.svm = region.Signals.(Signal.name).svm';
  idx = find((Test_Data.pos<Signal.lwin_big+1)|(Test_Data.pos>length(Test_Data.seq)-Signal.rwin_big-1));
  if length(idx)==length(Test_Data.pos), 
    % error('no examples in valid region found'); 
  end
  Test_Data.pos(idx)=[];
  Test_Data.svm(idx)=[];
  % load(PAR.fn_pos,'pos','seq','region_id','svm_num') ;
  if ~all(ismember(unique(Test_Data.seq),['ACGT']))
    warning('sequence alphabet different from ACGT; setting all Ns to A')
    Test_Data.seq(~ismember(Test_Data.seq,['ACGT']))='A';
  end
  if ~avg_all_svms
    idx = find(Test_Data.svm==which_set);
  else
    idx = 1:length(Test_Data.svm);
  end
  if isequal(Signal.method.predict_fct,'svm_content_predict')
    idx = [idx length(Test_Data.svm)];
  end
  if length(idx)==0, warning('no examples for specified svm found'); end
  Test_Data.pos = Test_Data.pos(idx);
  Test_Data.svm = Test_Data.svm(idx);
  if isfield(RPROC, 'divide4pred') && RPROC.divide4pred
    junk_length = ceil(length(Test_Data.pos)/RPROC.numjunks) ;
    idx_sub = (chunk_idx-1)*junk_length+1:min(chunk_idx*junk_length,length(Test_Data.pos));
    Test_Data.pos = Test_Data.pos(idx_sub);
    Test_Data.svm = Test_Data.svm(idx_sub);
    fn_pred = sprintf('%s%i.mat',fn_pred, chunk_idx);
  end
  %   Test_Data.region_id = Test_Data.region_id(idx);
end

Test_Data.PAR = P;

feval(Signal.method.predict_fct, Train, Test_Data,fn_SVMs,fn_pred,pred_from_seq, check_overlap);
