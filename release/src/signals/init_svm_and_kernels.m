function init_svm_and_kernels(method,pred_from_seq)
% init_svm_and_kernels(method,pred_from_seq)
  
  
  conservWdsWeighting =0;
  %------------------------------------------------------------------------------%
  % init svm
  %------------------------------------------------------------------------------%
  
  % === cleaning
%  sg( 'clean_features', 'TRAIN' );
%  sg( 'clean_preproc' );
%  sg( 'clean_kernel');
  %clear sg;
  
  %sg( 'loglevel', 'ALL');   % DEBUG
  if ~isfield(method, 'svm_name')
    warning( 'no svm specified, set to default: LIGHT')
    method.svm_name='LIGHT'
  end
  
  fprintf('LD_LIBRARY_PATH was set to : %s\n', getenv('LD_LIBRARY_PATH'))

  % === SVM
  sg( 'new_svm', method.svm_name);
  %sg( 'svm_epsilon', 1e-4 );
  %%%sg( 'svm_epsilon', PAR.svmEps);
  sg('svm_epsilon', 1e-4);
  %%%sg( 'send_command', sprintf( 'use_linadd %d', PAR.useLinadd ) );
  sg('use_linadd', true);
  %%%sg( 'svm_qpsize', PAR.qpSize );
  if ~isfield(method, 'qpsize')
    warning( 'qpsize not specified, set to default: 50')
    method.qpsize=50;
  end
  sg('svm_qpsize', method.qpsize) ; 
  if( conservWdsWeighting > 0 )
    sg('use_batch_computation', false);
  else
    sg('use_batch_computation', true);
  end
  
  % === parallel training?
  if isfield(method,'threads')
    sg('threads', method.threads);
  else
    fprintf('number of threads not specified, set to default: 16')
    sg('threads', 16);
  end;
  
  
  %------------------------------------------------------------------------------%
  % init kernels
  %------------------------------------------------------------------------------%
  
  nofKernels = length(method.kernels);
  if nofKernels>1 && ~isfield(method.kernels{2},'name') && ~isfield(method.kernels{2}, 'kernel_name')
    % for content sensors; several orders are trained simulatuously.
    % but this is done by concatination of the feature vector directly.
    nofKernels = 1;
  end
  
  assert( nofKernels >= 1 );
  if( nofKernels == 1 )
    addKernel = {'set_kernel'};
  else
    addKernel = {'add_kernel', 1};
    sg( 'set_kernel', 'COMBINED', 100);
  end;
  
  for k = 1:nofKernels
    kernel = method.kernels{k};
    % for compatibility with older versions
    if ~isfield(kernel,'kernel_name') && isfield(kernel, 'name')
      kernel.kernel_name = kernel.name;
    end
    fprintf(1, '%s: %s\n', sprintf('%s %1.1f', addKernel{:}), kernel.kernel_name)
    
    switch kernel.kernel_name
     case 'WD'
      if ~isfield(kernel, 'mismatch')
        kernel.mismatch=0;
      end
      assert(length(kernel.order)==1)
      assert(length(kernel.mismatch)==1)
      sg( addKernel{:}, 'WEIGHTEDDEGREE', 'CHAR', 10, kernel.order, kernel.mismatch, true, 1) ;
      
     case 'WDS'
      if ~isfield(kernel, 'mismatch')
        kernel.mismatch=0;
      end
      if ~isfield(kernel, 'shift')
        kernel.shift=0;
      end
      if ~isfield(kernel, 'shift_const')
        kernel.shift_const=0;
      end
      assert(length(kernel.order)==1)
      assert(length(kernel.mismatch)==1)
      len = kernel.rwin-kernel.lwin+1;
      shifts = ceil(kernel.shift * abs( (kernel.lwin:kernel.rwin)+1 ));
      shifts = shifts + kernel.shift_const;
      %shifts = sprintf( '%i ', shifts );
      sg( addKernel{:}, 'WEIGHTEDDEGREEPOS3', 'CHAR', 10, kernel.order, kernel.mismatch, len, 1, int32(shifts)) ;
      clear shifts len
      
     case 'SP'
      use_sign = false;
      % use_sign and prediction from position lists does not work
      %assert(~use_sign|~PAR.tasks.pred_from_seq)
      
      %addKernel_ = {'add_kernel', 0};
      
      sg( addKernel{:}, 'COMMSTRING', 'WORD', 10, use_sign);
      if pred_from_seq
        %sg('loglevel','ALL')
        sg( 'use_diagonal_speedup',true)
      end
     case 'RBF'
      if ~isfield(kernel, 'cache')
        kernel.cache = 40;
      end
      assert(length(kernel.sigma)==1)
      sg( addKernel{:}, 'GAUSSIAN', 'REAL', kernel.cache, kernel.sigma) ;
      
     case 'POLY'
      if ~isfield(kernel, 'cache')
        kernel.cache = 40;
      end
      if ~isfield(kernel, 'type')
        kernel.type = 'REAL';
      end 
      if ~isfield(kernel, 'norm')
        kernel.type = 'R';
      end
      assert(length(kernel.degree)==1)
      sg( addKernel{:}, 'POLY', Kernel.type, kernel.cache, kernel.degree, kernel.norm) ;
      
     case 'LINEAR'
      if ~isfield(kernel, 'cache')
        kernel.cache = 40;
      end
      if ~isfield(kernel, 'type')
        kernel.type = 'REAL';
      end
      sg( addKernel{:}, 'LINEAR', kernel.type, kernel.cache) ;
      
     otherwise
      fprintf(1, 'not yet implemented: kernel_name: %s\n', kernel.kernel_name)
    end;
    
  end;%loop over kernels
  
  if isfield(method, 'subkernel_weights') 
    assert(length(method.subkernel_weights)==nofKernels)
    sg('set_subkernel_weights', method.subkernel_weights);
  end
  
  % for compatibility with older versions
  if ~isfield(method.par_ms,'C') && isfield(method.kernels{1}, 'C')
    method.par_ms.C = method.kernels{1}.C;
  end
  sg('c', method.par_ms.C) ;
  
  disp('done with kernel setup');
