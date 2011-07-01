function [blocks,PAR,P] = pred_path(PAR, blocks)
% [blocks,PAR,P] = pred_path(PAR [, blocks])
%
%
check_split = 0 ;

if ~isfield(PAR, 'fn_log')
  fid=1;
PAR.evaluate = 0;
else  
  fid = fopen( sprintf('%spred.txt',PAR.fn_log),'a+') ;
  if fid<0
	fid = 1;
  end
end
fprintf(fid,'\n\npred_path \n');
fprintf(fid,'date: %s\n',date);
%fprintf(fid,'\n start \n');

if ~isfield(PAR, 'gen_paths_func')
  PAR.gen_paths_func = 'gen_paths' ; 
else
  PAR.gen_paths_func 
end ;

%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(PAR, 'fn_gene_predictor')
  if ~isfield(PAR, 'iteration')
    iteration = 1 ;
    while (1)
      fname = sprintf('%s_iteration%i.mat', PAR.FN.output_lsl.fn_train,iteration) ;
      if ~fexist(fname), iteration=iteration-1; break ; end 
      iteration=iteration+1 ;
    end ;
  else
    iteration = PAR.iteration;
  end ;
  fname = sprintf('%s_iteration%i.mat', PAR.FN.output_lsl.fn_train,iteration) ;
  fprintf('Loading iteration %i\n', iteration) ;
  if check_split,
    load(fname, 'QP', 'TRAIN_PAR', 'num_examples');
    TRAIN_PAR.model = fix_model(TRAIN_PAR.model) ;
    blocks_train = load_struct(sprintf('%s_iteration%i_blocks.mat', PAR.FN.output_lsl.fn_train,iteration)) ;
    blocks_train = fix_blocks(blocks_train, TRAIN_PAR.model) ;
    idx = [blocks_train.id] ;
    clear blocks_train
  else
    load(fname, 'QP', 'TRAIN_PAR', 'num_examples');
    TRAIN_PAR.model = fix_model(TRAIN_PAR.model) ;
  end ;
else 
  load(PAR.fn_gene_predictor, 'QP', 'TRAIN_PAR', 'num_examples');
  TRAIN_PAR.model = fix_model(TRAIN_PAR.model) ;
end

%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
  load(PAR.fn_blocks, 'blocks') ;
  blocks = fix_blocks(blocks, TRAIN_PAR.model) ;
end
if isfield(PAR,'fn_splits')
  load(PAR.fn_splits, 'blocks_split') ; 
  block_idx = [blocks_split{PAR.pred_split_idx}] ;
  if check_split
    if ~isempty(intersect(idx, block_idx))
      disp('test and training are overlapping') ;
      keyboard
      disp('correcting indices') ;
      block_idx = setdiff([blocks_split{:}], idx) ;
    end ;
  end ;
  blocks = blocks(block_idx) ;
end
fprintf('Predicting on %i blocks\n', length(blocks)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(TRAIN_PAR.model, 'use_conservation') && TRAIN_PAR.model.use_conservation 
  disp('Converting conservation levels') ;
  for id = 1:length(blocks),
    blocks(id).genestr_conserv = conservation2DNA_alphabet(blocks(id).conservation) ;
  end ;  
end ;

if isfield(PAR.LSL.method, 'add_lin_feat') && PAR.LSL.method.add_lin_feat
  assert(iscell(PAR.LSL.method.add_feat_fun))
  assert(iscell(PAR.FN.input_lsl.fn_lin_feat_data))
  assert(length(PAR.FN.input_lsl.fn_lin_feat_data)==length(PAR.LSL.method.add_feat_fun))
  for feat_type=1:length(PAR.LSL.method.add_feat_fun)
    fprintf('load linear features with function: %s from file: %s\n', PAR.LSL.method.add_feat_fun{feat_type}, PAR.FN.input_lsl.fn_lin_feat_data{feat_type}) ;
    blocks = feval(PAR.LSL.method.add_feat_fun{feat_type}, blocks, PAR.FN.input_lsl.fn_lin_feat_data{feat_type}, 80);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if isfield(PAR, 'use_gc') && ~PAR.use_gc
%  disp('Removing GC splice sites') ;
%  for id = 1:length(blocks),
%    idx=find(blocks(id).genestr(blocks(id).details.pos.don+1)=='T') ;
%    blocks(id).details.pos.don=blocks(id).details.pos.don(idx) ;
%    blocks(id).details.score.don=blocks(id).details.score.don(idx) ;
%    blocks(id).details.output.don=blocks(id).details.output.don(idx) ;
%  end ;  
%end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% avoid multiple loading of the big constant parameter file
global GLOBAL_cPAR ;
GLOBAL_cPAR=[] ;


cPAR=[] ;
cPAR.model = TRAIN_PAR.model ;

if isfield(PAR.model, 'intron_list'),
  cPAR.model.intron_list = PAR.model.intron_list ;
end ;

%cPAR.blocks_fname = ['~/' tempname '.mat'] ;
%cPAR.blocks = blocks ;
[xis, cPAR.parameters] = res2param(QP.res, num_examples, TRAIN_PAR.model, ...
                                   QP, 0, TRAIN_PAR.weight_names);
cPAR.parameters_fields = fieldnames(cPAR.parameters) ;

%% ADD offset to tss, tis, cdsStop or cleave transitions to enhance sensitivity
  %cPAR.parameters.transitions(find(cPAR.model.a_trans(2,:)==cPAR.model.state_ids.tss-1)) = cPAR.parameters.transitions(find(cPAR.model.a_trans(2,:)==cPAR.model.state_ids.tss-1))+PAR.LSL.method.tss_offset ;
if isfield(PAR, 'cleave_offset')
	tis_idx = find(cPAR.model.a_trans(2,:)==cPAR.model.state_ids.tis(1)-1);
	cPAR.parameters.transitions(tis_idx) = cPAR.parameters.transitions(tis_idx)+PAR.cleave_offset ;
end
P=[] ;
P.cPAR_fname = sprintf('~/tmp/%i-%1.10f.mat', round(rand(1)*1e8), now) ;
assert(~fexist(P.cPAR_fname)) ;
save(P.cPAR_fname, 'cPAR', '-v7') ;

%save_struct(cPAR.blocks_fname, blocks) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.viterbi_nbest = 1;

run_locally = rproc_policy('pred_path:gen_paths', [], length(blocks)) ;

[engine, environment] = determine_engine() ;

num = 0 ; model = TRAIN_PAR.model ;
jobinfo = rproc_empty(1) ;
if run_locally,
  fprintf(fid,'\nComputing %i dynamic programs ...',length(blocks));
else
  fprintf(fid,'\nSubmitting %i dynamic programs ...',length(blocks));
end ;
for block_id = 1:length(blocks),

  %idx = find(blocks(block_id).truth(1).segments(:,3)~=model.segments.intergenic) ;
  %max_len = max(500, max(blocks(block_id).truth(1).segments(idx,2)- ...
  %                       blocks(block_id).truth(1).segments(idx,1))) ;
  %max_len = max_len*1.2+1000 
  %P.max_len = max_len ;
  fprintf('%i ... ', block_id) ;

  P.block_id = block_id ;
  P.block = blocks(block_id) ;
  P.fn_save = sprintf('%s_viterbi_block%i.mat', PAR.fn_pred, blocks(block_id).id);
  system(sprintf('rm -f %s', P.fn_save)) ;

  %if fexist(P.fn_save),
    %continue; 
  %end 
  P.prediction = 1 ;
  P.FN = PAR.FN;
  options = PAR.RPROC.options;
  if isequal(engine, 'matlab') && isequal(environment, 'galaxy'),
    options.force_octave = 1 ;
  end ;
  options.identifier = ['G'  num2str(block_id)] ;
  options.start_dir = get_base_dir();
  options.express = 0;
  options.verbosity = 2 ;
  options.addpaths = {shogun_settings(), fileparts(which('gen_paths'))};

  if ~run_locally,
    [mem_req, time_req, opts]=rproc_memtime_policy('gen_paths', length(P.block.seq), options) ;
    num = num+1 ;
    if ~isfield(options, 'immediately_bg') || options.immediately_bg==0,
      if PAR.RPROC.exm_per_batch==1,
        jobinfo(num) = rproc(PAR.gen_paths_func, P, mem_req, opts, time_req) ;
      else
        jobinfo(num) = rproc_create(PAR.gen_paths_func, P, mem_req, opts, time_req) ;
      end ;
    else
      jobinfo(num) = rproc(PAR.gen_paths_func, P, mem_req, opts, time_req) ;
    end ;
  else
    pred = feval(PAR.gen_paths_func, P) ;
  end ;
end
fprintf('Done\n') ;

if ~run_locally,
  if ~isfield(options, 'immediately_bg') || options.immediately_bg==0,
    if PAR.RPROC.exm_per_batch==1,
        [jobinfo, num_crashed] = rproc_wait(jobinfo, 60, 1, -1, 1);
    else
      [jobinfo, meta_jobinfo] = rproc_submit_batch(jobinfo, PAR.RPROC.exm_per_batch) ;
      [meta_jobinfo,num_crashed] = rproc_wait(meta_jobinfo, 60, 1, -1, 1);
    end ;
  else
    [jobinfo,num_crashed] = rproc_wait(jobinfo, 60, 1, -1, 1);
  end ;
end ;

if nargout>0
	blocks(1).pred = struct;
	for block_id = 1:length(blocks),
	  block = blocks(block_id) ;
	  while ~fexist(sprintf('%s_viterbi_block%i.mat', PAR.fn_pred, block.id ))
	    fprintf('.') ;
	    pause(10)
	  end
	  load(sprintf('%s_viterbi_block%i.mat', PAR.fn_pred, block.id ),'pred');   
	  blocks(block_id).prediction = pred; 
	end
end
unix(['rm -f ' P.cPAR_fname]) ;

%save([PAR.fn_pred 'blocks'], 'blocks', 'PAR', 'P') ;
if fid>2
  fclose(fid);
end
