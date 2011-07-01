function [RFV, RSV, train_test_overlap] = svm_content_predict(Train, Data, svmFileName,outFileName,pred_from_seq, check_overlap);
% [RFV, RSV] = svm_predict(Train, Data, svmFileName,outFileName);
%
%------------------------------------------------------------------------------%
% deal with parameters 
%------------------------------------------------------------------------------%
% for usage with rproc

if nargin==1
  Data  = Train.Data;
  svmFileName = Train.svmFileName;
  outFileName = Train.outFileName;
  Train = Train.Train;
end

method = Train.PAR.Signal.method ;

if ~exist('check_overlap', 'var'),
  check_overlap=1 ;
end ;

sg_ver = sg('get_version');

assert(size(Data.pos,1)==1|size(Data.pos,2)==2|isempty(Data.pos))
if size(Data.pos,1)~=1
  Data.pos = Data.pos';
end

if length(Train.w)~=sum(4.^[3:6]),
  fprintf('Assumed degrees 3:6 are not used in trained content predictor (dim %i)\n', length(Train.w)) ;
  Train.w(sum(4.^[3:6]))=0 ;
end ;

assert(isequal(sort(Data.pos), Data.pos));
sg('init_dyn_prog', 1);
sg('precompute_content_svms', Data.seq, int32(Data.pos-1), Train.w);
output = sg('get_lin_feat');
% output = output - output(end)/Data.pos(end).*Data.pos;

%------------------------------------------------------------------------------%
% save output
%------------------------------------------------------------------------------%

PAR = Data.PAR;
PAR.Signal = Train.PAR.Signal;

if isfield(Data, 'LT')
  RSV = calcrocscore(output, Data.LT); 
  RFV= calcrfcscore(output, Data.LT);
  LT = Data.LT;
  %   save( outFileName, 'output','LT','PAR', 'svmFileName','RSV','RFV' );
  POS = Data.Pos;
  fprintf( 'save output to %s\n',outFileName );
  inventory = {'output', 'POS'};
  %save( outFileName, 'inventory', 'output','LT','POS', 'Train', 'svmFileName','RSV','RFV','sg_ver','PAR' );
  save( outFileName, 'inventory', 'output','LT','POS', 'svmFileName','RSV','RFV','sg_ver','PAR' );
else
  pos = Data.pos;
  fprintf( 'save output to %s\n',outFileName );
  inventory = {'output', 'pos'};
  %save( outFileName, 'inventory', 'output', 'pos', 'Train', 'svmFileName','sg_ver','PAR' );
  save( outFileName, 'inventory', 'output', 'pos', 'svmFileName','sg_ver','PAR' );
end

return
%%%%TEST



%L = load(['galaxy/testit/content_cds_exon_classifier_files/examples_settings=2_split=1.mat']);
%LL = load('galaxy/testit/content_cds_exon_classifier_files/output_SVMCOMBINED/MS/partition=1.mat')
%
%for i=1:length(L.LT)
%  sg('init_dyn_prog', 1);
%  sg('precompute_content_svms', L.XT{i},int32(0:length(L.XT{i})-1),LL.Train.w);
%  output{i} = sg('get_lin_feat');
%  output1(i) = output{i}(end);
%end
%
%%%%
%
%
%XT = zeros(5440,length(L.XT));
%method = LL.Train.PAR.Signal.method;
%for j=1:length(L.XT) 
%  dim = 0;
%  for o=1:length(method.kernels)
%    [mask hist] = compute_mers(L.XT{j}, method.kernels{o}.wordlen, method.kernels{o}.stepping, method.kernels{o}.offset);
%    XT(dim+[1:4^method.kernels{o}.wordlen],j) = hist';
%    dim = dim+ 4^method.kernels{o}.wordlen;
%  end
%end
%init_svm_and_kernels(method,0);
%sg('clean_features', 'TRAIN');
%D.XT = LL.Train.XT;
%D.LT = LL.Train.LT;
%setFeatures(method, D,'TRAIN',0,0);
%clear D
%
%sg('clean_features', 'TEST' );
%D.XT = XT;
%setFeatures(method, D,'TEST',0,0);
%
%sg('set_svm',LL.Train.b,LL.Train.alphas);
%output2 = sg( 'svm_classify' );
