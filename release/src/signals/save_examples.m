function save_examples(blocks,split,P)
%save_examples(blocks,split,P)
%
% INPUT 	regions with field Signals with field `Signal.name` with fields:
%			-pos
%			-anno_cover
%			-est_cover
%			-cdna_cover
%	PAR with fields
%			-organism
%			-Signal.filter_label
%			-Signal.name
%			-Signal.consensus
%			-FN.input_sig.fn_examples
%			-FN.input.fn_split
%			-FN.input.fn_regions


%only for RPROC
if nargin<3
	%addpath ..
	% use the rproc.addpath way to add paths in client processes!!
	%paths
	P = blocks.P;
	split = blocks.split;
	blocks = blocks.blocks;
end

Signal = P.Signal ; 
num_splits = P.num_splits;

fn_training_blocks = P.fn_training_blocks;
fn_candsites = P.fn_candsites;
fn_filter_settings = P.fn_filter_settings;
fn_examples = P.fn_examples;
fn_example_statistics = P.fn_example_statistics;

Signal_name = Signal.name;
lwin = Signal.lwin_big;
rwin = Signal.rwin_big;

print=1; % print output only once

fields = fieldnames(Signal.filter_label);
fprintf( 1, 'save_examples... (%d fields, %d splits)\n', length(fields), num_splits );
for ff = 1:length(fields)
	 for ss=1:num_splits
		%---------------------------------------------
		%CREATE OUTPUT FILENAME
		%----------------------------------------------	
		filename = gen_example_filename(Signal,fn_filter_settings,fn_examples,ff,ss);
		if ~fexist(filename)
			num_examples = 0;
			POS.fn_training_blocks = fn_training_blocks;	 
			split_reg_ids = getfield(split,sprintf('split_%i',ss));
			blocks_idx = find(ismember([blocks.id],split_reg_ids));
			split_regs = blocks(blocks_idx);
			clear blocks_idx split_reg_ids
			
			%if ~isfield(split_regs,'seq') 
			%	fprintf(1,'loading sequence...\n')
			%	split_regs = load_sequence(split_regs);
			%	fprintf(1,'\n')
			%end
			
			%---------------------------------------------
			% retreive label information
			%---------------------------------------------
			
			% if ~isfield(split_regs, 'Signals')
			%	 temp_region = retrieve_signals(split_regs(1),fn_candsites,Signal_name,{'label',Sig.Conf_names{:}});
			% end
			
			
			%---------------------------------------------
			%----heuristic--->estimate size of chararray; 
			% -some positions are filtered out by the filerfunction: Signal.filter_fct
			% -if the size is to small then this function is very time-consuming
			% -otherwise it is very space-consuming
			%---------------------------------------------
			%total_num_pos = retrieve_num_pos(split_regs, fn_candsites) ;
			%num_pos = ceil(total_num_pos) ;
		 
			%--------------------------------------------------
			% determine the maximal number of training examples
			% for the given sensor
			%--------------------------------------------------
			t = cputime;
			num_pos_ex = 0;
			num_neg_ex = 0;
			for r=1:length(split_regs)
				if ~isfield(split_regs, 'Signals')
					temp_region = retrieve_signals(split_regs(r),fn_candsites,Signal_name,{'label'});
				end
				num_pos_ex = num_pos_ex+sum(temp_region.Signals.(Signal_name).label == 1);
				num_neg_ex = num_neg_ex+sum(temp_region.Signals.(Signal_name).label == -1);
			end
			fprintf(1,'it took %fs to determine the true number of positive and negative examples\n', cputime-t);
			vec=subsample_policy([num_pos_ex num_neg_ex], Signal_name);
			max_num_pos = vec(1);
			max_num_neg = vec(2);
			reduction_factor_pos = max_num_pos/num_pos_ex;
			reduction_factor_neg = max_num_neg/num_neg_ex;
			num_pos = max_num_pos+max_num_neg;

			if print==1
				fprintf(1,'reducing number of positive examples by factor %f\n',reduction_factor_pos);
				fprintf(1,'reducing number of negative examples by factor %f\n',reduction_factor_neg);
				print=0;
			end

			XT=zeros([lwin+rwin,num_pos],'int8');
			LT=zeros(1,num_pos);
			POS.pos=zeros(1,num_pos);
			POS.region_id=zeros(1,num_pos);

			if isfield(Signal.filter_label.(fields{ff}), 'USE_ALL')&&...
						Signal.filter_label.(fields{ff}).USE_ALL
				Signal.Conf_names = {} ;
			end ;
			
			t = cputime;
			fprintf(1,'starting with split %i\n',ss)
			for r=1:length(split_regs)
				if mod(r,100)==0,
					fprintf('region %i/%i\r', r, length(split_regs)) ;
				end ;
				
				if ~isfield(split_regs, 'Signals')
					temp_region = retrieve_signals(split_regs(r),fn_candsites,Signal_name,{'label',Signal.Conf_names{:}});
				end
				
				sig_field = getfield(temp_region.Signals,Signal_name);
				if isfield(Signal.filter_label.(fields{ff}),'use_label')
					sig_field.label = sig_field.(Signal.filter_label.(fields{ff}).use_label);
				end
				
				if ~isfield(split_regs,'seq') 		
					if isfield(Signal.filter_label.(fields{ff}), 'USE_ALL')&&...
								Signal.filter_label.(fields{ff}).USE_ALL
						if temp_region.strand == '+',
							loffset=0 ; % remember the offset for later
							tmp = max(1,temp_region.start-lwin) ;
							loffset = temp_region.start - tmp ;
							temp_region.start = tmp ;
							temp_region.stop = temp_region.stop + rwin ;
						else
							tmp = temp_region.start-rwin ;
							temp_region.start = tmp ;
							temp_region.stop = temp_region.stop + lwin ;
							loffset = lwin ;
						end ;
					end ;
					temp_region = load_sequence(temp_region);
				end
				seq = temp_region.seq;
				
				%---------------------
				%FILTER EXAMPLES
				%---------------------
				if ~isfield(Signal.filter_label.(fields{ff}), 'USE_ALL')||...
							~Signal.filter_label.(fields{ff}).USE_ALL
					if ~isfield(Signal.filter_label.(fields{ff}),'conf_pos') ...
								&& isfield(Signal.filter_label.(fields{ff}),'conf')
						error('fix this: filtering will not work')
					end
					idx = feval(Signal.filter_fct, temp_region, Signal, fields{ff},1);
				else
					sig_field.pos = sig_field.pos + loffset ;
					idx = find(sig_field.pos>lwin&sig_field.pos<length(seq)-rwin+1);
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
						idx = intersect(idx, un);
					end
				end
				if isempty(idx); continue; end
				
				%--------------------
				%CREATE XTs
				%--------------------
				XT_temp = zeros(lwin+rwin,length(idx),'int8');
				rid_temp= zeros(1,length(idx))+temp_region.id;
				for row=1:lwin+rwin
					XT_temp(row,:)= seq(sig_field.pos(idx)-lwin+row-1);
				end%loop over rows
				LT_temp= sig_field.label(idx)';
				pos_temp = sig_field.pos(idx)';
				assert(size(LT,2)==num_examples||LT(num_examples+1)==0)
				assert(num_examples==0||LT(num_examples)~=0)
				
				%[size(XT,2) num_examples+length(idx) ceil(size(XT,2)*r/length(split_regs))]
				if size(XT,2)<num_examples+length(idx), 
					warning('extending XT') ;
					POS.pos(ceil(size(XT,2)*1.2)) = 0 ;
					POS.region_id(ceil(size(XT,2)*1.2)) = 0 ;
					LT(ceil(size(XT,2)*1.2)) = 0 ;
					XT(1,ceil(size(XT,2)*1.2)) = 0 ;
				end ;
				
				XT(:, num_examples+1:num_examples+length(idx))= XT_temp;	
				LT(num_examples+1:num_examples+length(idx))= LT_temp;
				POS.pos(num_examples+1:num_examples+length(idx))= pos_temp;
				POS.region_id(num_examples+1:num_examples+length(idx))=rid_temp;
				num_examples = num_examples+length(idx);
			end%loop over blocks
			fprintf(1, 'time required for one split: %.2f\n',cputime-t)
			fprintf(1,'save file %s\n',filename)
			%=== remove zeros
			assert(all(LT(num_examples+1:end)==0))
			LT = LT(1:num_examples);
			XT = XT(:,1:num_examples);
			POS.pos = POS.pos(1:num_examples);
			POS.region_id = POS.region_id(1:num_examples);
			
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
			% exclude examples containing non ACGT-Symbols 
			%-----------------------------------------------
			%xtpos=XT~='A'&XT~='C'&XT~='G'&XT~='T';
			%save space
			xtpos=XT~='A';
			for s='CGT'
				xtpos=xtpos&XT~=s;
			end
			cols = find(sum(xtpos));
			clear xtpos
			if length(cols)>0
				warning(sprintf('%i examples containing non ACGT symbols excluded',length(cols)))
				XT(:,cols)=[];
				LT(cols)=[];
				POS.pos(cols)=[];
				POS.region_id(cols)=[];
			end
			%----------------------------------------------- 
			% save matrix in pieces into one file
			%----------------------------------------------- 
			assert(sum(LT==1)>0);
			fprintf( 'trying to save %s\n', filename );
			save_matrix(filename, XT, 'XT');
			%save(filename,'-append','LT','POS','PAR');
			save_append(filename, 1,'LT', LT, 'POS', POS, 'P', P);
			% save_append(filename, 1,'LT', LT, 'POS', POS, 'PAR', PAR);
		end%if ~fexist(filename)
	 end%loop over splits
	 
	 stat = save_example_statistics(Signal, fn_filter_settings, ...
																	fn_examples, fn_example_statistics, ...
																	num_splits);
end%loop over filter_label fields

