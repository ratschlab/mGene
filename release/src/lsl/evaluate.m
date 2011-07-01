function res = evaluate(fn_blocks, fn_pred_iter, model)



%iteration1 = find_iteration(fn_pred);
%l = load(sprintf('%siteration_%i/blocks.mat', PAR1.FN.output_lsl.fn_pred,iteration1)) ;
%blocks1 = load_struct(PAR1.FN.input_lsl.fn_val_blocks, 'blocks');
if ischar(fn_blocks) && ischar(fn_pred_iter),
  blocks1 = load_struct(fn_blocks, 'blocks');

  for j=1:length(blocks1)
    fname = sprintf('%s_viterbi_block%i.mat',fn_pred_iter,blocks1(j).id);
    load(fname,'pred');
    blocks1(j).pred = pred;
    clear pred
  end
elseif isstruct(fn_blocks) && ischar(fn_pred_iter),
	blocks1 = fn_blocks; 
	clear fn_blocks
	for j=1:length(blocks1)
		fname = sprintf('%s_viterbi_block%i.mat',fn_pred_iter,blocks1(j).id);
		load(fname,'pred');
		blocks1(j).pred = pred;
		clear pred
	end
else
  blocks1 = fn_blocks ;
  assert(isempty(fn_pred_iter)) ;
end ;
%blocks1 = l.blocks;
%fn_pred_iter = [fn_pred 'iteration_' num2str(iteration1) '/'];

%blocks1([blocks1.strand]=='-')=[];%used for a check
%blocks1([blocks1.strand]=='+')=[];%used for a check
%perm = randperm(length(blocks1));
%blocks1 = blocks1(perm(1:floor(length(blocks1)*7/10)));

  res_all.cds_nucleotides.num_obsv =0; 
  res_all.cds_nucleotides.num_pred =0; 
  res_all.cds_nucleotides.num_corr =0; 
  res_all.cds_exons.num_obsv       =0; 
  res_all.cds_exons.num_pred       =0; 
  res_all.cds_exons.num_corr       =0; 
  res_all.cds_transcripts.num_obsv =0; 
  res_all.cds_transcripts.num_pred =0; 
  res_all.cds_transcripts.num_corr =0; 

%transcripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(blocks1)
  segments = blocks1(j).truth(1).segments;
  res_all.cds_transcripts.num_obsv =res_all.cds_transcripts.num_obsv+max(segments(:,4)); 
  if ~isempty(blocks1(j).pred.genes{1})
    res_all.cds_transcripts.num_pred =res_all.cds_transcripts.num_pred+length(blocks1(j).pred.genes); 
  end
  for k=1:max(segments(:,4))
    for l=1:length(blocks1(j).pred.genes)
      if isempty(blocks1(j).pred.genes{1})
        continue;
      end
      trans1 = blocks1(j).pred.genes{l}(blocks1(j).pred.genes{l}(:,3)==model.segments.cds_exon,:);
      trans2 = segments(segments(:,3)==model.segments.cds_exon&segments(:,4)==k,1:3);
      if isequal(trans1,trans2)
         res_all.cds_transcripts.num_corr = res_all.cds_transcripts.num_corr +1;
      end
    end
  end
end

%exons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(blocks1)
  all_exons1 = [];
  all_exons2 = blocks1(j).truth(1).segments(blocks1(j).truth(1).segments(:,3)==model.segments.cds_exon, 1:3);
  for k=1:length(blocks1(j).pred.genes)
    if isempty(blocks1(j).pred.genes{1})
      continue;
    end
    trans1 = blocks1(j).pred.genes{k}(blocks1(j).pred.genes{k}(:,3)==model.segments.cds_exon,:);
    all_exons1 = [all_exons1;trans1];
  end

  common_exons = intersect(all_exons1, all_exons2, 'rows');
  res_all.cds_exons.num_pred=res_all.cds_exons.num_pred+size(all_exons2,1);
  res_all.cds_exons.num_obsv=res_all.cds_exons.num_obsv+size(all_exons1,1);
  res_all.cds_exons.num_corr=res_all.cds_exons.num_corr+size(common_exons,1);
  if size(common_exons,1)<size(all_exons2,1)
    %keyboard
  end
end


% Compute sensitivity and specificity
%res_all.cds_nucleotides.SN = res_all.cds_nucleotides.num_corr  /res_all.cds_nucleotides.num_obsv;
%res_all.cds_nucleotides.SP = res_all.cds_nucleotides.num_corr  /res_all.cds_nucleotides.num_pred;

res_all.cds_exons.SN       = res_all.cds_exons.num_corr        /res_all.cds_exons.num_obsv;
res_all.cds_exons.SP       = res_all.cds_exons.num_corr        /res_all.cds_exons.num_pred;


res_all.cds_transcripts.SN = res_all.cds_transcripts.num_corr  /res_all.cds_transcripts.num_obsv;
res_all.cds_transcripts.SP = res_all.cds_transcripts.num_corr  /res_all.cds_transcripts.num_pred;
res_all.cds_transcripts.F  = 2*res_all.cds_transcripts.SN*res_all.cds_transcripts.SP/(res_all.cds_transcripts.SP+res_all.cds_transcripts.SN)

res = res_all;

%save(sprintf('%sresults',PAR.fn_pred),'res');

return
