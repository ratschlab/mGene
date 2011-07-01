function [ alphas, b, betas, PAR1 ] = svm_train(Data);

P = Data.P;

method  = P.Signal.method;
fn_SVMs = P.fn_SVMs;
pred_from_seq = P.pred_from_seq;
use_preproc = P.use_preproc;


%------------------------------------------------------------------------------%
% SVM already trained?
%------------------------------------------------------------------------------%

%if( ~isfield(PAR,'conservWeighting') )
method.conservWeighting = 0;
%end;
if( ~exist('testDataPar','var') )
  testDataPar = [];
end;


%------------------------------------------------------------------------------%
% train SVM
%------------------------------------------------------------------------------%


if( method.conservWeighting == 0 )
  Data.W = [];
else
  Data.W = res.W;
  assert( ~isempty(W) );
end;

% keyboard
init_svm_and_kernels(method,pred_from_seq);


if isequal(method.feature_type,'double') && method.normalize_features
  idx_rm = find(std(Data.XT')<1e-8);
  Data.XT(idx_rm,:)  = [];
  mean_tr = mean(Data.XT');
  Data.XT = Data.XT - repmat(mean_tr,size(Data.XT,2),1)';
  std_tr = std(Data.XT');
  Data.XT = Data.XT./repmat(std_tr,size(Data.XT,2),1)';
end

setFeatures(method, Data,'TRAIN',pred_from_seq,use_preproc);

m = length( Data.LT );
% PAR1 = P;
% PAR1.POS = Data.Pos;
% PAR1.lenTrain = size( Data.XT, 1 );

sg_ver=sg('get_version');

%sg('loglevel', 'DEBUG')

% === train
nofKernels = length(method.kernels);
t0 = cputime;
sg('set_kernel_optimization_type', 'FASTBUTMEMHUNGRY' );
%sg( 'set_kernel_optimization_type', 'SLOWBUTMEMEFFICIENT' );
if( method.use_mkl )
  sg('use_mkl', true );
  sg('mkl_parameters', 1e-3, 0 );
  betas = sg( 'get_subkernel_weights' );
  assert( size(betas) == [nofKernels,1] );
end;
sg('svm_train' );
fprintf( '\n' );
t = cputime - t0;

%------------------------------------------------------------------------------%
% get alphas and betas
%------------------------------------------------------------------------------%

% === alphas
[ b, alphas ] = sg( 'get_svm' );


% === betas (kernel weights)
if(method.use_mkl)
  betas = sg( 'get_subkernel_weights' );
  assert( size(betas) == [nofKernels,1] );
  betas
else
  betas = ones( nofKernels, 1 );
end;


%------------------------------------------------------------------------------%
% save SVM
%------------------------------------------------------------------------------%
Train.alphas = alphas;
Train.b = b;
Train.PAR = P;
Train.POS = Data.Pos;
sup_vec_idx = alphas(:,2)+1;
Train.XT = Data.XT(:,sup_vec_idx);
Train.LT = Data.LT(sup_vec_idx);

w = zeros(size(Train.XT,1),1);
for a=1:length(Train.LT)
  w = w+Train.alphas(a,1).*Train.XT(:,a);
end
Train.w = w;

Train.POS.pos = Data.Pos.pos(sup_vec_idx);
Train.POS.region_id = Data.Pos.region_id(sup_vec_idx);

if isequal(method.feature_type,'double') && method.normalize_features
  Train.idx_rm = idx_rm
  Train.mean = mean_tr;
  Train.std = std_tr;
end

clear P Data
%sg( 'clean_features', 'TRAIN' );
%sg( 'clean_features', 'TEST' );

Train.time=t;
Train.sg_ver = sg_ver ; 
fprintf('save: %s',fn_SVMs)
save(fn_SVMs, '-V7', 'Train')





