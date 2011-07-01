function setFeatures(method, Data, svmMode, pred_from_seq, use_preproc);
% setFeatures(method,Data, svmMode,pred_from_seq,use_preproc);

  
%------------------------------------------------------------------------------%
% check
%------------------------------------------------------------------------------%

assert( strcmp(svmMode,'TRAIN') || strcmp(svmMode,'TEST') );

sg('clean_features', svmMode)

%assert( isscalar(Data.pos0) );
%assert( size(XT,2) > 0 );
%assert( ischar( XT ) );

if nargin<5
  use_preproc=1;
end
if ~isfield(method,'conservWdsWeighting')
  conservWdsWeighting=0;
else
  conservWdsWeighting= method.conservWdsWeighting;
end
%------------------------------------------------------------------------------%
% pass data
%------------------------------------------------------------------------------%

%sg('loglevel', 'ALL');
if pred_from_seq
  seq = Data.seq;
  idx = find(seq~='A'&seq~='C'&seq~='G'&seq~='T');
  if ~isempty(idx)
    warning('non AGCT-Symbols found; changed to As')
    seq(idx)='A';
  end
  pos = Data.pos;
else
  XT = Data.XT;
  %pos0=Data.pos0;
  if isfield(Data, 'LT')
    LT = Data.LT;
  end
  if isfield(Data, 'W')
    W = Data.W;
  end
end

%==============================================================================%
% set labels
%==============================================================================%
if strcmp(svmMode,'TRAIN')
  sg( 'set_labels', 'TRAIN', LT ); 
end
%==============================================================================%
% set features
%==============================================================================%



nofKernels = length(method.kernels);
if nofKernels>1 && ~isfield(method.kernels{2},'name') && ~isfield(method.kernels{2}, 'kernel_name')
  % for content sensors; several orders are trained simulatuously.
  % but this is done by concatination of the feature vector directly.
  nofKernels = 1;
end
assert( nofKernels >= 1 );
if( nofKernels == 1 )
  addFeats = 'set_features';
else
  addFeats = 'add_features';
end;


for k = 1:nofKernels
  kernel = method.kernels{k};
  % for compatibility with older versions
  if ~isfield(kernel, 'kernel_name') && isfield(kernel,'name')
    kernel.kernel_name = kernel.name;
  end
  if ~pred_from_seq   
    if ~isfield(method,'center') || isempty(method.center)
      select= [1:size(XT,1)];
    else
      if ~isfield(kernel,'lwin') 
        warning('no lwin specified for kernel')
        kernel.lwin = -method.center;
      end
      if ~isfield(kernel,'rwin') 
        warning('no rwin specified for kernel')
        kernel.rwin = -method.center+size(XT,1)-1;
      end
      select = (kernel.lwin:kernel.rwin)+method.center+1;
    end
    if ~isequal(select,[1:size(XT,1)])
      XT_temp = XT(select,:);
    else
      XT_temp = XT;
    end
    %keyboard;
  else
    %positiones should be a 0-based index
    pos_temp=pos+kernel.lwin-1;
    winsize_temp=kernel.rwin-kernel.lwin+1;
    assert(length(pos_temp)>0)
    assert(length(seq)>max(pos_temp)+winsize_temp)
    assert(winsize_temp>0)

    % work-around for from_position_list bug
    if length(pos_temp)==1,
      pos_temp = [pos_temp pos_temp] ;
    end ;
  end
  %------------------------------------------------------------------------------%
  % kernel: WD or WDS
  %------------------------------------------------------------------------------%

  if strcmp(kernel.kernel_name,'WD') || strcmp(kernel.kernel_name,'WDS')
    if ~pred_from_seq
      fprintf( '--> XT[%d:%d], num_examples=%i,  order=%d\n', min(select), max(select),size(XT,2), kernel.order );
      sg( addFeats, svmMode, XT_temp, 'DNA' );
    else
      %sg( addFeats, svmMode,{seq}, 'DNA' );
      %sg('from_position_list','TEST', winsize_temp, int32(pos_temp));
      sg( addFeats, svmMode, {seq}, 'DNA', 'from_position_list', winsize_temp, int32(pos_temp));
    end
  end;
  

  %------------------------------------------------------------------------------%
  % kernel: Spectrum
  %------------------------------------------------------------------------------%

  if strcmp(kernel.kernel_name,'SP' )
    if~pred_from_seq 
      fprintf( '--> XT[%d:%d], wordlen=%d\n', min(select), max(select), kernel.wordlen );
      sg( addFeats, svmMode, XT_temp, 'DNA');
    
      if( ~isfield(kernel,'gaps') )
        sg( 'convert', svmMode, 'STRING', 'CHAR', 'STRING', 'WORD', kernel.wordlen, kernel.wordlen-1 ) ;
      else
        fprintf('--> convert\n')
        sg( 'convert', svmMode, 'STRING', 'CHAR', 'STRING', 'WORD', kernel.wordlen, kernel.wordlen-1, kernel.gaps );
      end;
      
      if( strcmp(svmMode,'TRAIN') && use_preproc )
        fprintf('--> clean_preproc\n')
        sg( 'clean_preproc' );
        fprintf('--> add_preproc SORTWORDSTRING\n')
        sg( 'add_preproc', 'SORTWORDSTRING' );
        fprintf('--> attach_preproc %s\n', svmMode)
        sg( 'attach_preproc', svmMode );
      end;
    else
      % from_position_list doesnot work for training for spectrum kernel  
      assert(strcmp(svmMode,'TEST')) ;
      %sg( addFeats, svmMode, {seq}, 'DNA' );
      sg( addFeats, svmMode, {seq}, 'DNA', 'from_position_list', winsize_temp-kernel.wordlen+1, int32(pos_temp+kernel.wordlen-1) );
      % window is cut out from sequence if convertion is made after calling 'from_position_list
      % if convertion is done at this place the window size must be reduced by -order-1;
      % (not yet tested) this would save a lot of time and space. 
      % 
      if( ~isfield(kernel,'gaps') )
        sg( 'convert', svmMode, 'STRING', 'CHAR', 'STRING', 'WORD', kernel.wordlen, 0);
      else
        fprintf('--> convert\n')
        sg( 'convert', svmMode, 'STRING', 'CHAR', 'STRING', 'WORD', kernel.wordlen, 0, kernel.gaps );
      end;
      %sg('from_position_list','TEST', winsize_temp- kernel.wordlen+1, int32(pos_temp+kernel.wordlen-1));
    end
  end;
 
  %------------------------------------------------------------------------------%
  % kernel:  RBF 
  %------------------------------------------------------------------------------%
 
  if strcmp(kernel.kernel_name,'RBF')
    fprintf( '--> XT[%d:%d], num_examples=%i,  sigma=%d\n', min(select), max(select),size(XT,2), kernel.sigma );
    sg( addFeats, svmMode, XT_temp, 'REAL');    
  end;
  
  %------------------------------------------------------------------------------%
  % kernel:  POLYNOMIAL
  %------------------------------------------------------------------------------%
 
  if strcmp(kernel.kernel_name,'POLY')
    fprintf( '--> XT[%d:%d], num_examples=%i,  degree=%d\n', min(select), max(select),size(XT,2), kernel.degree );
    sg( addFeats, svmMode, XT_temp, 'REAL');    
  end;
  
  %------------------------------------------------------------------------------%
  % kernel:  LINEAR
  %------------------------------------------------------------------------------%
  
  if strcmp(kernel.kernel_name,'LINEAR')
    fprintf( '--> XT[%d:%d], num_examples=%i,\n', min(select), max(select),size(XT,2) );
    sg( addFeats, svmMode, XT_temp, 'REAL');    
  end;
end%loop over kernels
%------------------------------------------------------------------------------%
% done
%------------------------------------------------------------------------------%

if( use_preproc & strcmp(svmMode,'TEST')  & ~pred_from_seq )
   fprintf('--> attach_preproc %s', svmMode)
  sg( 'attach_preproc', svmMode );
end;

clear XT;
%sg( 'init_kernel', svmMode );

% --- set feature weights
if( conservWdsWeighting > 0 )
  W1 = W( select1, : );
  if( conservWdsWeightNorming )
    fprintf( '--> conserv(WDS): using normalized feature weights!\n' );
    meanW = mean( W1, 1 );
    W1 = W1 .* repmat( 1./meanW, size(W1,1), 1 );
  else
    fprintf( '--> conserv(WDS): using raw feature weights!\n' );
    W1 = 3/size(W1,1) * W1;
  end;
  sg( 'set_subkernel_weights_combined', W1, 1 );
  clear W1;
end;


