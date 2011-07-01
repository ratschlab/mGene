function x = save_content_examples(blocks,split,P)
%save_examples(blocks,split,P)
%
% INPUT 	regions with field Signals with field `Signal.name` with fields:
%			-pos
%			-anno_cover
%			-est_cover
%			-cdna_cover
%	PAR with fields
%			-organism
%			-Content.filter_label
%			-Content.name
%			-Content.consensus
%			-FN.input_sig.fn_examples
%			-FN.input.fn_split
%			-FN.input.fn_regions

%dummy return argument
x=0;

%only for RPROC
if nargin<3
  %addpath ..
  % use the rproc.addpath way to add paths in client processes!!
  %paths
  P = blocks.P;
  split = blocks.split;
  blocks = blocks.blocks;
end

Content = P.Content ; 
num_splits = P.num_splits;

fn_training_blocks = P.fn_training_blocks;
fn_candsites = P.fn_candsites;
fn_filter_settings = P.fn_filter_settings;
fn_examples = P.fn_examples;
fn_example_statistics = P.fn_example_statistics;

Content_name = Content.name;

fields = fieldnames(Content.filter_label);
fprintf( 1, 'save_examples... (%d fields, %d splits)\n', length(fields), num_splits );

min_seg_length = 7;

for ff = 1:length(fields)
   for ss=1:num_splits
    %---------------------------------------------
    %CREATE OUTPUT FILENAME
    %----------------------------------------------  
    filename = gen_example_filename(Content,fn_filter_settings,fn_examples,ff,ss);
    if ~fexist(filename)
      num_examples = 0;
      POS.fn_training_blocks = fn_training_blocks;	 
      split_reg_ids = getfield(split,sprintf('split_%i',ss));
      blocks_idx = find(ismember([blocks.id],split_reg_ids));
      split_regs = blocks(blocks_idx);
      clear blocks_idx split_reg_ids
     

      %-------------------------------------------------- 
      % determine the number of positive and negative 
      % examples in this split
      %-------------------------------------------------- 
      t = cputime;
      num_pos_ex = 0;
      num_neg_ex = 0;
      for r=1:length(split_regs)
        temp_region = retrieve_signals(split_regs(r),fn_candsites,Content_name,{'label'});
	num_pos_ex = num_pos_ex+sum(temp_region.Signals.(Content_name).label == 1);
	num_neg_ex = num_neg_ex+sum(temp_region.Signals.(Content_name).label == -1);
      end
      fprintf(1,'it took %.2fs to determine the number of examples\n', cputime-t);
      vec=subsample_policy([num_pos_ex num_neg_ex], Content_name);
      max_num_pos = vec(1);
      max_num_neg = vec(2);
      reduction_factor_pos = max_num_pos/num_pos_ex;
      reduction_factor_neg = max_num_neg/num_neg_ex;
      num_pos = max_num_pos+max_num_neg;

 
      %-------------------------------------------------- 
      % initialize arrays 
      %-------------------------------------------------- 
      XT = cell(1,num_pos);
      LT = zeros(1,num_pos);
      POS.pos=zeros(1,num_pos);
      POS.region_id=zeros(1,num_pos);
      
      if isfield(Content.filter_label.(fields{ff}), 'USE_ALL')&&...
            Content.filter_label.(fields{ff}).USE_ALL
        Content.Conf_names = {} ;
      end ;
      
      t = cputime;
      fprintf(1,'starting with split %i\n',ss)
      for r=1:length(split_regs)
        if mod(r,100)==0,
          fprintf('region %i/%i\r', r, length(split_regs)) ;
        end ;

        
        if fn_candsites(1)=='~'
          fn_candsites = sprintf('%s/%s', getenv('HOME'), fn_candsites(2:end));
        end
        if ~isfield(split_regs, 'Contents')
          temp_region = retrieve_signals(split_regs(r),fn_candsites,Content_name,{'pos2','label',Content.Conf_names{:}});
        end
        
        
        sig_field = getfield(temp_region.Signals,Content_name);
        if isfield(Content.filter_label.(fields{ff}),'use_label')
          sig_field.label = sig_field.(Content.filter_label.(fields{ff}).use_label);
        end
        
        if ~isfield(split_regs,'seq') 		
          temp_region = load_sequence(temp_region);
        end
        seq = temp_region.seq;

        %---------------------
        %FILTER EXAMPLES
        %---------------------
        if ~isfield(Content.filter_label.(fields{ff}), 'USE_ALL')||...
              ~Content.filter_label.(fields{ff}).USE_ALL
          if ~isfield(Content.filter_label.(fields{ff}),'conf_pos') ...
                && isfield(Content.filter_label.(fields{ff}),'conf')
            error('fix this: filtering will not work')
          end
          idx = feval(Content.filter_fct, temp_region, Content, fields{ff},1);
        else
          idx = 1:length(sig_field.pos);
        end
        idx1 = find(sig_field.pos2<=length(seq));

        idx = intersect(idx,idx1);

	pos_idx = find(sig_field.label==1);
        if reduction_factor_pos<1
	  keep = randperm(length(pos_idx));
	  pos_idx = pos_idx(keep(1:floor(length(pos_idx)*reduction_factor_pos)));
	end
	neg_idx = find(sig_field.label==-1);
	if reduction_factor_neg<1
	  keep = randperm(length(neg_idx));
	  neg_idx = neg_idx(keep(1:floor(length(neg_idx)*reduction_factor_neg)));
	end
	if reduction_factor_pos<1||reduction_factor_neg<1
          if ~isempty(pos_idx) && ~isempty(neg_idx),
            un=[pos_idx', neg_idx'] ;
          elseif isempty(pos_idx),
            un = neg_idx' ;
          else
            un = pos_idx' ;
          end ;
	  num_tmp = length(idx);
	  idx = intersect(idx, un);
	  %fprintf('reduce the number of examples from %i to %i (pos:%i->%i, neg:%i->%i)\n', num_tmp, length(idx), sum(sig_field.label==1), length(pos_idx), sum(sig_field.label==-1), length(neg_idx)); 
          clear num_tmp
	end

        if isempty(idx); continue; end
        
        %--------------------
        %CREATE XTs
        %--------------------
        XT_temp = cell(1,length(idx));
        rid_temp = zeros(1,length(idx))+temp_region.id;
        for num=1:length(idx)
          if any(~ismember(seq(max(1,sig_field.pos(idx(num))):sig_field.pos2(idx(num))),'ACGT'))
            %warning(sprintf(' examples containing non ACGT symbols excluded'))
            continue
          end
          if length(seq(max(1,sig_field.pos(idx(num))):sig_field.pos2(idx(num))))<min_seg_length
            %warning(sprintf('very short segments are excluded (<5bp)'))
            continue
          end
          if isempty(seq(max(1,sig_field.pos(idx(num))):sig_field.pos2(idx(num))))
            warning(sprintf('empty range'))
          end
          XT_temp{num} = seq(max(1,sig_field.pos(idx(num))):sig_field.pos2(idx(num)));
        end%loop over rows
        LT_temp= sig_field.label(idx)';
        pos_temp = sig_field.pos(idx)';
        pos2_temp = sig_field.pos2(idx)';
        assert(size(LT,2)==num_examples||LT(num_examples+1)==0)
        assert(num_examples==0||LT(num_examples)~=0)
        
        %[size(XT,2) num_examples+length(idx) ceil(size(XT,2)*r/length(split_regs))]
        if size(XT)<num_examples+length(idx), 
          warning('extending XT') ;
          LT(ceil(size(XT,2)*1.2)) = 0 ;
          
          XT(1,ceil(size(XT,2)*1.2)) = 0 ;
        end ;
        XT(num_examples+1:num_examples+length(idx)) =  XT_temp(:);	
        LT(num_examples+1:num_examples+length(idx)) =  LT_temp;
        POS.pos(num_examples+1:num_examples+length(idx))= pos_temp;
        POS.region_id(num_examples+1:num_examples+length(idx))=rid_temp;
        num_examples = num_examples+length(idx);
      end%loop over blocks
      
      
     
      fprintf(1, 'time required for one split: %.2f\n',cputime-t)
      fprintf(1,'save file %s\n',filename)
      %=== remove zeros
      assert(all(LT(num_examples+1:end)==0))
      LT = LT(1:num_examples);
      XT = XT(1:num_examples);
      
      
      % balance length distribution of negative examples
      %---------------------------------------------------------------
      fprintf(1,'start fit_length_hist\n');
      bins = 100;
      reduce_pos=0;% throw away positive examples if necessary to achieve similar histograms
      
      neg_examples = XT(LT==-1);
      pos_examples = XT(LT==1);
      % if ~use_frame
      neg_examples = fit_length_hist(pos_examples,neg_examples,bins, reduce_pos);
      % end
      LT = [-ones(size(neg_examples)) ones(size(pos_examples))];
      XT = {neg_examples{:} pos_examples{:}};
      idx = randperm(length(LT));
      LT = LT(idx);
      XT = XT(idx);
      %----------------------------------------------- 
      % remove examples that appear more than once
      %----------------------------------------------- 
      if 0
        [tmp1 idx tmp2] = unique(XT','rows');
        idx = idx(randperm(length(idx)));
        XT = XT(:,idx);
        LT = LT(idx);
        POS.pos = POS.pos(idx);
        POS.region_id = POS.region_id(idx);
      end
     
      %----------------------------------------------- 
      % save matrix in pieces into one file
      %----------------------------------------------- 
      % assert(sum(LT==1)>0);
      % fprintf( 'trying to save %s\n', filename );
      save(filename,'XT');
      % save_matrix(filename, XT, 'XT');
      % save(filename,'-append','LT','POS','PAR');
      save_append(filename, 1,'LT', LT, 'POS', POS, 'P', P);
      % save_append(filename, 1,'LT', LT, 'POS', POS, 'PAR', PAR);
    
    end%if ~fexist(filename)
    
   end%loop over splits
   
   stat = save_example_statistics(Content, fn_filter_settings, ...
                                  fn_examples, fn_example_statistics, ...
                                  num_splits);
end%loop over filter_label fields

