function [blocks,res]=predict(PAR, blocks)
% [blocks,res]=predict(PAR, blocks)

%addpath(PAR.FN.source_directory)
paths

%if nargin<3 || strcmp(set,'test')
%  testset=1;
%else
%  testset=0;
%end
%
%%PAR.submit=1;
%if ~isfield(PAR, 'iteration')||isempty(PAR.iteration) 
%  PAR.iteration = find_iteration(PAR);
%end
%if testset&&(~isfield(PAR, 'fn_pred')||isempty(PAR.fn_pred))
%  PAR.fn_pred = [PAR.FN.output_lsl.fn_pred 'iteration_' num2str(PAR.iteration) '/'];
%  PAR.fn_blocks = PAR.FN.input_lsl.fn_val_blocks;
%elseif ~testset&&(~isfield(PAR, 'fn_pred')||isempty(PAR.fn_pred))
%  PAR.fn_pred = [PAR.FN.output_lsl.fn_pred 'trainset_iteration_' num2str(PAR.iteration) '/'];
%  PAR.fn_blocks = PAR.FN.input_lsl.fn_training_blocks;
%end

PAR.model = fix_model(PAR.model) ;

if nargin<2
  blocks = load_struct(PAR.fn_blocks, 'blocks');
  blocks = fix_blocks(blocks, PAR.model) ;
end

unix(sprintf('mkdir -p %s',PAR.fn_pred));
%blocks = pred_path(PAR, blocks) ;
pred_path(PAR, blocks) ;

res = [] ;

if isfield(PAR, 'evaluate') && ~PAR.evaluate
  return
end

%blocks = [];

%res = evaluate(PAR.fn_blocks, PAR.fn_pred, PAR.model)
res = evaluate(blocks, PAR.fn_pred, PAR.model)
fprintf('save results to %sresults.mat\n',PAR.fn_pred) 
save(sprintf('%sresults.mat',PAR.fn_pred), 'res', '-v7');

%printf('Nucleotide level accuracy:\n') ;
%print_struct(1, res.cds_nucleotides) ;

fprintf('Exon level accuracy:\n')
print_struct(1, res.cds_exons) ;

fprintf('Transcript level accuracy:\n')
print_struct(1, res.cds_transcripts) ;

fprintf('\n\n\n') ;

[engine]=determine_engine() 
if isequal(engine, 'octave'),
  struct_levels_to_print(2) ;
end ;

res


if isequal(engine, 'octave'),
  struct_levels_to_print(0) ;
end ;

return

genome_info = init_genome(PAR.FN.input.fn_genome_config);

for j=1:length(blocks)
  load(sprintf('%s_viterbi_block%i',PAR.fn_pred,blocks(j).id),'pred');
  blocks(j).pred = pred;
  clear pred
end

if 1
  load(PAR.FN.output.fn_genes_merged,'genes');
else
  load(PAR.FN.output.fn_genes_anno,'genes');
  for j=1:length(genes)
    genes(j).anno_id = genes(j).id;
    if strcmp(genes(j).transcripts{1},'dummy')
      genes(j).mask_on = 1;
    end
%    genes(j).mask_on = max(genes(j).transcript_status)<1;
  end
end
  res_all.cds_nucleotides.num_obsv =0; 
  res_all.cds_nucleotides.num_pred =0; 
  res_all.cds_nucleotides.num_corr =0; 
  res_all.cds_exons.num_obsv       =0; 
  res_all.cds_exons.num_pred       =0; 
  res_all.cds_exons.num_corr       =0; 
  res_all.cds_transcripts.num_obsv =0; 
  res_all.cds_transcripts.num_pred =0; 
  res_all.cds_transcripts.num_corr =0; 

for j=1:length(blocks)
  if ~isempty(blocks(j).pred.genes{1})
    genes_pred = blocks2genes(blocks(j),PAR.model,'pred',genome_info);
  else
    genes_pred = init_genes(1);
  end
  for k=1:length(genes_pred)
    genes_pred(k).id = k;
    %if isempty(genes_pred(j).exons{1})
    %  warning('exons field empty');
    %end
    %if isempty(genes_pred(j).cds_exons{1})
    %  warning('exons field empty');
    %end
  end
  [res,map] = evaluate_predictions(blocks(j),genes,genes_pred,genome_info,0);
  res_all = addresults(res_all, res);
end

% Compute sensitivity and specificity
%res_all.nucleotides.SN     = res_all.nucleotides.num_corr      /res_all.nucleotides.num_obsv;
%res_all.nucleotides.SP     = res_all.nucleotides.num_corr      /res_all.nucleotides.num_pred;

res_all.cds_nucleotides.SN = res_all.cds_nucleotides.num_corr  /res_all.cds_nucleotides.num_obsv;
res_all.cds_nucleotides.SP = res_all.cds_nucleotides.num_corr  /res_all.cds_nucleotides.num_pred;

res_all.cds_exons.SN       = res_all.cds_exons.num_corr        /res_all.cds_exons.num_obsv;
res_all.cds_exons.SP       = res_all.cds_exons.num_corr        /res_all.cds_exons.num_pred;

%res_all.exons.SN           = res_all.exons.num_corr            /res_all.exons.num_obsv;
%res_all.exons.SP           = res_all.exons.num_corr            /res_all.exons.num_pred;
%
%res_all.transcripts.SN     = res_all.transcripts.num_corr      /res_all.transcripts.num_obsv;
%res_all.transcripts.SP     = res_all.transcripts.num_corr      /res_all.transcripts.num_pred;
%
res_all.cds_transcripts.SN = res_all.cds_transcripts.num_corr  /res_all.cds_transcripts.num_obsv;
res_all.cds_transcripts.SP = res_all.cds_transcripts.num_corr  /res_all.cds_transcripts.num_pred;

res = res_all;

ret = 0%res.cds_transcripts.mean_SNSP;

save(sprintf('%sresults.mat',PAR.fn_pred), 'res', 'ret', '-v7');

return

function res = addresults(r1, r2)

  res.cds_nucleotides.num_obsv 	= r1.cds_nucleotides.num_obsv	+r2.cds_nucleotides.num_obsv;
  res.cds_nucleotides.num_pred 	= r1.cds_nucleotides.num_pred	+r2.cds_nucleotides.num_pred;
  res.cds_nucleotides.num_corr 	= r1.cds_nucleotides.num_corr	+r2.cds_nucleotides.num_corr;
  res.cds_exons.num_obsv 	= r1.cds_exons.num_obsv		+r2.cds_exons.num_obsv;
  res.cds_exons.num_pred 	= r1.cds_exons.num_pred		+r2.cds_exons.num_pred;
  res.cds_exons.num_corr 	= r1.cds_exons.num_corr		+r2.cds_exons.num_corr;
  res.cds_transcripts.num_obsv 	= r1.cds_transcripts.num_obsv	+r2.cds_transcripts.num_obsv;
  res.cds_transcripts.num_pred 	= r1.cds_transcripts.num_pred	+r2.cds_transcripts.num_pred;
  res.cds_transcripts.num_corr 	= r1.cds_transcripts.num_corr	+r2.cds_transcripts.num_corr;
return


function iteration = find_iteration(PAR)

iteration = 1 ;
while (1)
  fname = sprintf('%s_iteration%i.mat', PAR.FN.output_lsl.fn_train,iteration) ;
  if ~fexist(fname), iteration=iteration-1 ; break ; end ;
  iteration=iteration+1 ;
end ;
return
