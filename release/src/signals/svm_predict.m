function [RFV, RSV, train_test_overlap] = svm_predict(Train, Data, svmFileName, outFileName, pred_from_seq, check_overlap);
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
  
  fprintf('working for outputs to be written to %s\n', outFileName);
  
  method = Train.PAR.Signal.method ;
  
  if ~exist('check_overlap', 'var'),
    check_overlap=1 ;
  end ;
  
  if( ~isfield(method,'conservWeighting') )
    method.conservWeighting = 0;
  end;
  if( method.conservWeighting == 0 )
    W1 = [];
    meanCons = [];
  else
    W1 = res.W;
    assert( ~isempty(W1) );
    meanCons = mean( W1, 1 );
  end;
  
  sg('loglevel', 'DEBUG')

  %------------------------------------------------------------------------------%
  % init svm and set kernel parameters
  %------------------------------------------------------------------------------%
  init_svm_and_kernels(method,pred_from_seq);
  sg_ver = sg('get_version');
  %------------------------------------------------------------------------------%
  % init training features with labels and examples that only contain those of 
  % support vectors
  %------------------------------------------------------------------------------%
  %sg('clean_features', 'TRAIN');
  % D.pos0 = method.center;
  D.XT = Train.XT;
  D.LT = Train.LT;
  
  setFeatures(method, D,'TRAIN',0,0);

  % this helped to avoid the error "Kernel not initialized with training examples"
  % I guess this is due to the init_kernel call that is done in get_kernel_matrix
  %%km=sg('get_kernel_matrix', 'TRAIN');


  clear D km 
  
  if pred_from_seq
	max_size = 5e7;
	if size(Data.seq, 2)>max_size
		num_splits = floor(size(Data.seq, 2)/max_size);
		bins = 0:max_size:size(Data.seq, 2);
		bins(end) = size(Data.seq, 2);
	else
		num_splits = 1;
		bins = [0 size(Data.seq, 2)];
	end

	win_size = 2000;% make shure there is enough sequence to predict on all positions in split
	output = [];
	for split = 1:num_splits
	  
	  if num_splits==1
      	D.seq = Data.seq(bins(split)+1:bins(split+1));
	  elseif split==1
      	D.seq = Data.seq(bins(split)+1:bins(split+1)+win_size);
	  elseif split==num_splits
      	D.seq = Data.seq(bins(split)+1-win_size:bins(split+1));
	  else
        D.seq = Data.seq(bins(split)+1-win_size:bins(split+1)+win_size);
      end

      assert(size(Data.pos,1)==1|size(Data.pos,2)==2|isempty(Data.pos))
      if size(Data.pos,1)~=1
        D.pos = Data.pos';
      else
        D.pos = Data.pos;
      end

      % select positions to predict
	  split_idx = find(D.pos>=bins(split)+1&D.pos<=bins(split+1));
      D.pos = D.pos(split_idx);
      if split>1
        D.pos = D.pos-bins(split)+win_size;
      end

      if ~isempty(D.pos)
        setFeatures(method, D, 'TEST', pred_from_seq, 0);
        D_pos = D.pos ;
        clear D
        assert(size(Train.alphas, 1)==size(Train.XT,2))
        Train.alphas(:,2)= (0:size(Train.alphas, 1)-1)';
        sg('set_svm', Train.b, Train.alphas);
        
        if( method.use_mkl )
          sg( 'set_subkernel_weights', betas );
        end;
        output_split = sg( 'svm_classify' );
        
        % work-around for from_position_list bug
        if length(D_pos) == 1 && length(output_split)==2,
          output_split=output_split(1) ;
        end ;
        clear D_pos
      else
        output_split = [];
      end
      output = [output output_split];
      sg('clean_features', 'TEST' );
	end
  else
    %------------------------------------------------------------------------------%
    % init test features with examples if pred_from_seq==0 or with seq and positions
    % in the case of genome wide predictions
    %------------------------------------------------------------------------------%
    %sg('clean_features', 'TEST' );
    D.XT = Data.XT;
    if isfield(method, 'feature_type') && isequal(method.feature_type,'double') && method.normalize_features
      D.XT(Train.idx_rm,:)  = [];
      D.XT = D.XT - repmat(Train.mean,size(D.XT,2),1)';
      D.XT = D.XT./repmat(Train.std,size(D.XT,2),1)';
    end
    % D.pos0 = method.center;
    D.pos = Data.Pos;
    
    if ~isempty(D.pos)
      setFeatures(method, D, 'TEST', pred_from_seq, 0);
      D_pos = D.pos ;
      clear D
      %------------------------------------------------------------------------------%
      % setup svm with trainingdata
      %------------------------------------------------------------------------------%
      %Train.XT should contain only support vector examples
      assert(size(Train.alphas, 1)==size(Train.XT,2))
      % alphas are in sparse format col_1=value; col_2=0-based index
      % change index to 0:size(alphas,2)-1 consistent with XTs  
      Train.alphas(:,2)= (0:size(Train.alphas, 1)-1)';
      sg('set_svm', Train.b, Train.alphas);
      
      %------------------------------------------------------------------------------%
      % classify
      %------------------------------------------------------------------------------%
      if( method.use_mkl )
        sg( 'set_subkernel_weights', betas );
      end;
      output = sg( 'svm_classify' );
      
      % work-around for from_position_list bug
      if length(D_pos) == 1 && length(output)==2,
        output=output(1) ;
      end ;
      clear D_pos
    else
      output = [];
    end
  end
  % check for overlaps between training and testing data
  
  if check_overlap && ~isequal(determine_engine, 'octave'),
    
    train_subset = randperm(size(Train.XT, 2)) ;
    train_subset = train_subset(1:min(length(train_subset), 100)) ;
    
    if ~pred_from_seq
      [tmp, idx, idx_] = intersect(Train.XT(:, train_subset)', Data.XT', 'rows') ;
    else
      P=[] ; 
      for i=1:length(train_subset), 
        P = [P strfind(Data.seq, Train.XT(:,train_subset(i))') + method.center] ; 
      end ;
      [tmp,idx,idx_]=intersect(P, Data.pos);
    end ;
    if length(idx)>0,
      warning('Overlap between Train and Test sets discovered (%i/%i sequences)\n', length(idx), length(train_subset)) ;
      if length(idx)/length(train_subset)>0.1,
        warning('Large overlap between Train and Test sets discovered (%i/%i sequences)\n', length(idx), length(train_subset)) ;
      end ;
    end ;
  end 
      
      
      
  %------------------------------------------------------------------------------%
  % save output
  %------------------------------------------------------------------------------%
  
  PAR = Data.PAR;
  PAR.Signal = Train.PAR.Signal;
  
  if isfield(Data, 'LT')
    RSV = calcrocscore(output, Data.LT); 
    if isnan(RSV), 
      RSV=0.5 ;
    end ;
    RFV= calcrfcscore(output, Data.LT);
    if isnan(RFV),
      RFV=0 ;
    end;
    LT = Data.LT;
    %   save( outFileName, 'output','LT','PAR', 'svmFileName','RSV','RFV' );
    POS = Data.Pos;
    fprintf( 'save output to %s\n',outFileName );
    inventory = {'output', 'POS'};
    %save( outFileName, '-V7', 'inventory', 'output','LT','POS', 'Train', 'svmFileName','RSV','RFV','sg_ver','PAR' );
    save( outFileName, '-V7', 'inventory', 'output','LT','POS', 'svmFileName','RSV','RFV','sg_ver','PAR' );
  else
    pos = Data.pos;
    fprintf( 'save output to %s\n',outFileName );
    %save( outFileName, '-V7', 'output', 'pos', 'Train', 'svmFileName','sg_ver','PAR' );
    %if isfield(PAR.Signal, 'export_settings') && isfield(PAR.Signal.export_settings, 'resolution') && PAR.Signal.export_settings.resolution>1
    %  pos_short = zeros(1,floor(length(pos)/2));
    %  output_short = zeros(1,floor(length(pos)/2));
    %  % reduce number of predictions (take the better one out of two) 
    %  for j=2:2:length(pos)
    %    if output(j-1)>output(j)
    %      output_short(j/2) = output(j-1);
    %      pos_short(j/2) = pos(j-1);
    %    else
    %      output_short(j/2) = output(j);
    %      pos_short(j/2) = pos(j);
    %    end
    %  end
    %  output = output_short;
    %  pos = pos_short;
    %end
    
    inventory = {'output', 'pos'};
    assert(isequal(size(pos),size(output)) || (isempty(pos) && isempty(output)));
    save( outFileName, '-V7', 'inventory', 'output', 'pos', 'svmFileName','sg_ver','PAR' );
  end
      
      
