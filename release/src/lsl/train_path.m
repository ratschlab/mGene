function [fn_predictor, iteration] = train_path(PAR)
% [fn_predictor, iteration] = train_path(PAR)

fid = 1% fopen( sprintf('%strain.txt',PAR.FN.output_lsl.fn_log),'a+') ;
fprintf(fid,'\n\n train_path \n');
fprintf(fid,'date: %s\n',date);
fprintf(fid,'\n initializing \n');


PAR.model = fix_model(PAR.model) ;
model = PAR.model ;

PAR.weight_names = init_weight_names(model) ;

method = PAR.LSL.method;

lpenv = init_solver ;


%% determint iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(PAR,'iteration'),
  iteration = PAR.iteration 
  PAR = rmfield(PAR,'iteration') ;
else
  iteration = 1 ;
  while (1)
    fname = sprintf('%s_iteration%i.mat', PAR.FN.output_lsl.fn_train,iteration) ;
    if ~fexist(fname), iteration=iteration-1 ; break ; end ;
    iteration=iteration+1 ;
  end ;  
end ;


%% load blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = sprintf('%s_iteration%i.mat', PAR.FN.output_lsl.fn_train,iteration) 
fn_parameters = sprintf('%s_parameters.mat', PAR.FN.output_lsl.fn_train)
if fexist(fname),
  fprintf('restarting at iteration %i\n', iteration) ;
  load(fname, 'TRAIN_PAR', 'QP','num_examples') ;
  %blocks_train = load_struct(sprintf('%s_iteration%i_blocks.mat', PAR.FN.output_lsl.fn_train,iteration)) ;
  blocks_train = load_struct(PAR.FN.input_lsl.fn_training_blocks, 'blocks') ;
  blocks_train = fix_blocks(blocks_train, model) ;
  %init_path(PAR, blocks_train, lpenv, PAR.weight_names);
  assert(fexist(fn_parameters)==1)
else
  assert(iteration==0)
  blocks_train = load_struct(PAR.FN.input_lsl.fn_training_blocks, 'blocks') ;
  blocks_train = fix_blocks(blocks_train, model) ;
  if isfield(PAR.LSL, 'train_blocks_num') && PAR.LSL.train_blocks_num<length(blocks_train)
	fprintf('reduce number of blocks to %i\n', PAR.LSL.train_blocks_num) ;
	blocks_train = blocks_train(1:PAR.LSL.train_blocks_num);
  end
  %if isfield(PAR.LSL.method, 'add_lin_feat')&&PAR.LSL.method.add_lin_feat
  %  assert(iscell(PAR.LSL.method.add_feat_fun))
  %  assert(iscell(PAR.FN.input_lsl.fn_lin_feat_data))
  %  assert(length(PAR.FN.input_lsl.fn_lin_feat_data)==length(PAR.LSL.method.add_feat_fun))
  %  for feat_type=1:length(PAR.LSL.method.add_feat_fun)
  %    fprintf('load linear features with function: %s from file: %s\n', PAR.LSL.method.add_feat_fun{feat_type}, PAR.FN.input_lsl.fn_lin_feat_data{feat_type}) ;
  %    blocks_train = feval(PAR.LSL.method.add_feat_fun{feat_type}, blocks_train, PAR.FN.input_lsl.fn_lin_feat_data{feat_type}, 80);
  %  end
  %end
  %blocks_train = initialize_blocks(blocks_train, PAR.model, PAR.LSL.method.max_neg, PAR.FN);
  blocks_train = fix_blocks(blocks_train, model) ;
  [QP,parameters] = init_path(PAR, blocks_train, lpenv, PAR.weight_names);
  num_examples = length(blocks_train) ;
  solve_cnt = 0;
  solve_time = 0;
  decode_time = 0;
  param_history={};
  save('-v7', fn_parameters, 'decode_time', 'solve_time', 'solve_cnt', 'param_history')
  fprintf(fid,'\nStarting from iteration %i  \n\n',iteration);
end ;

QP.lpenv = lpenv;
global last_p_lp last_num_con
last_p_lp=[] ;
last_num_con=0 ;

%% sort blocks such that bigger jobs are submitted first
%[tmp sort_idx] = sort([blocks_train.stop]-[blocks_train.start], 'descend');
len = zeros(1, length(blocks_train));
for j = 1:length(blocks_train)
	len(j) = length(blocks_train(j).seq);%start and stop for these artificial blocks is arbitrary
end
[tmp sort_idx] = sort(len, 'descend');
blocks_train = blocks_train(sort_idx);

% make sure parameters fit to QP.res
fprintf(fid,'call res2param\n');
[xis,parameters] = res2param(QP.res, num_examples, model, QP, 0, PAR.weight_names);

% rproc options
[shogun_extra_path, shogun_envstr] = shogun_settings();
options = PAR.RPROC.options;
options.express = 0;
options.start_dir = get_base_dir();
options.verbosity = 0;
options.hard_time_limit= 10;
options.addpaths = {shogun_settings(), fileparts(which('gen_paths'))};
options.envstr = shogun_envstr ;

%[engine, environment] = determine_engine() ;
if 1%isequal(engine, 'matlab') && isequal(environment, 'galaxy'),
  options.force_octave = 1 ;
  % options.immediately_bg = 1 ;
end ;

% save blocks
blocks_name = ['~' tempname]
nice_mkdir(fileparts(blocks_name));
for i=1:length(blocks_train),
  block = blocks_train(i);
  save(sprintf('%s_block_%i',blocks_name, block.id), 'block');
end ;
clear block

% remove fields that are not used for training
run_locally=0;
run_locally= run_locally |  rproc_policy('train_path:gen_paths:approx', [], length(blocks_train)/10)
run_locally= run_locally |  rproc_policy('train_path:gen_paths:noapprox', [], length(blocks_train))
if ~run_locally
	blocks_train = remove_unused_fields(blocks_train);
	never_submit = 1;
end

% initialize approximation flags
for i=1:length(blocks_train),
  blocks_train(i).pred_use_approximation = 1 ;
  %blocks_train(i).pred_use_approximation = 0 ;
end ;
clear block

jobinfo_predict = rproc_empty(0) ;

%% this is the folder where the predictions are saved to.
%% 
if (PAR.FN.output_lsl.fn_train_block_preds(end) == '/')
	PAR.FN.output_lsl.fn_train_block_preds = PAR.FN.output_lsl.fn_train_block_preds(1:end-1);
end
fn_train_block_preds_orig = PAR.FN.output_lsl.fn_train_block_preds;

% as soon as the number of new constraints dropped below some 
% threshold do not use approximation any more
never_approx = 0;


%% start with iteration over training examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
old_OBJ = -100 ;
idx_blocks = 1:num_examples;
viterbi_nbest = ones(1,num_examples) ;
while iteration < method.max_num_iterations % | abs(old_OBJ-OBJ)<1e-3,
  iteration = iteration + 1 ;
  num_new_constraints_iteration=0;
  OBJ(iteration) = sum(QP.res.*QP.f)+0.5*QP.res'*QP.Q*QP.res ;
  fprintf(fid,'\nStarting iteration: %i\n-----------------------\n',iteration);

  % save predictions for each iteration to enable reusing them for another 
  % gene finder instance 
  PAR.FN.output_lsl.fn_train_block_preds = sprintf('%s_iter%i/', fn_train_block_preds_orig, iteration);
  nice_mkdir(PAR.FN.output_lsl.fn_train_block_preds);

  num_correct = 0;
  num_crashed = 0;
  jobinfo = rproc_empty(0) ;

  for block_id_start = 1%:method.exm_per_solve:num_examples,
    %block_ids = block_id_start:min(block_id_start+method.exm_per_solve-1, num_examples) ;
    block_ids = 1:num_examples ;
    
    %%%%%%%%%%%%%%%%%%%%%
    global GLOBAL_cPAR ;
    GLOBAL_cPAR=[] ;
    
    cPAR=[] ;
    cPAR.model = model ;
    %cPAR.blocks_fname = ['~/' tempname '.mat'] ;
    [xis, cPAR.parameters] = res2param(QP.res, num_examples, model, QP, 0, PAR.weight_names);
    cPAR.parameters_fields = fieldnames(cPAR.parameters) ;
    
    P=[] ;
    P.cPAR_fname = sprintf('~/tmp/%i-%1.10f.mat', round(rand(1)*1e8), now) ;
    assert(~fexist(P.cPAR_fname)) ;
    save(P.cPAR_fname, 'cPAR', '-v7') ;

    %save_struct(cPAR.blocks_fname, blocks_train) ;

    %%%%%%%%%%%%%%%%%%%%%
    
    % only submit jobs if sufficiently many blocks need the full computation
    frac_full = mean([blocks_train.pred_use_approximation]==0);
    compute_full_gen_paths = frac_full > 0.5 ;

    fprintf('Dynamic programming approximation appropriate for %i/%i jobs\n', sum([blocks_train.pred_use_approximation]==1), length(blocks_train)) ;
    if compute_full_gen_paths,
      fprintf('=> Compute exact dynamic programs for all blocks\n') ;
    else
      fprintf('=> Compute approximate dynamic programs for all blocks\n') ;
    end ;

    if compute_full_gen_paths,
      [run_locally]= rproc_policy('train_path:gen_paths:noapprox', [], length(block_ids)) ;
    else
      [run_locally]= rproc_policy('train_path:gen_paths:approx', [], length(block_ids)/10) ;
    end ;
    
    if run_locally,
      fprintf('\nRunning %i dynamic programs: ', length(block_ids)) ;
    else
      fprintf('\nSubmitting %i dynamic programs: ', length(block_ids)) ;
    end ;
    jobinfo=rproc_empty(0) ; ii=0 ;

	%% start loop over blocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for block_id = block_ids,
      idx = find(blocks_train(block_id).truth(1).segments(:,3)~=model.segments.intergenic) ;

      P.block_id = block_id ;
      P.FN = PAR.FN;
      P.block_name = sprintf('%s_block_%i',blocks_name, blocks_train(block_id).id); 
      P.viterbi_nbest = viterbi_nbest(block_id) ;
      P.fn_save = sprintf('%spred_block%i.mat', PAR.FN.output_lsl.fn_train_block_preds, blocks_train(block_id).id); 
      P.prediction = 0 ;

      % get a quick approximation of the best path 
      if ~compute_full_gen_paths 
        P.approximate = 1 ;
        
        if rand(1)<0.05, % randomly disable some blocks ... to guarantee convergence
          blocks_train(block_id).pred_use_approximation = 0 ;
        end ;
      else
        unix(['rm -f ' P.fn_save]) ;
        P.approximate = 0 ;
        blocks_train(block_id).pred_use_approximation = 1 ;
        %blocks_train(block_id).pred_use_approximation = 0 ;
      end ;
      
      options.identifier = ['G'  num2str(block_id)] ;
      ii=ii+1 ;
      if mod(ii, 10)==0,
        fprintf('%i ... ', ii) ;
      end ;

      if ~run_locally,
        
        [mem_req, time_req, opts]= rproc_memtime_policy('gen_paths', length(blocks_train(block_id).seq), options, 0) ;
        
        if ~isfield(options, 'immediately_bg') || options.immediately_bg==0,
          if 1%PAR.RPROC.exm_per_batch==1,
            jobinfo(end+1) = rproc('gen_paths', P, mem_req, opts, time_req) ;
          else
            jobinfo(end+1) = rproc_create('gen_paths', P, mem_req, opts, time_req) ;
          end ;
        else
          jobinfo(end+1) = rproc('gen_paths', P, mem_req, opts, time_req) ;
        end ;
      else
        P.block = blocks_train(block_id) ;
        gen_paths(P) ;
      end
    end ;
    fprintf('Done.\n') ;
    
    if ~isempty(jobinfo),
      if ~isfield(options, 'immediately_bg') || options.immediately_bg==0,
        if 1%PAR.RPROC.exm_per_batch==1,
          %[jobinfo, num_crashed] = rproc_wait(jobinfo, 15, 1, -1, 1);
          [jobinfo, num_crashed] = rproc_wait(jobinfo, 30, 1, -1, 1);
        else
          [jobinfo, meta_jobinfo] = rproc_submit_batch(jobinfo, PAR.RPROC.exm_per_batch) ;
          [meta_jobinfo,num_crashed] = rproc_wait(jobinfo, 15, 1, -1, 1);
        end ;
        %[meta_jobinfo,num_crashed] = rproc_wait(jobinfo, 20, 1, 0);
        %unix(['rm -f ' P.cPAR_fname ' ' cPAR.blocks_fname]) ;
      else
        [jobinfo] = rproc_wait(jobinfo, 30, 1, -1, 1);
      end ;
    end ;

    blocks = blocks_train(block_ids) ;
    idx_temp = idx_blocks(block_ids) ;
    
    %% add new constraints 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    [blocks, QP, viterbi_nbest(block_ids), num_new_constraints, time] = add_constraints(blocks, idx_temp, num_examples, PAR, QP, viterbi_nbest(block_ids));
    blocks_train(block_ids) = blocks ;

    if ~isempty(jobinfo),
      rproc_cleanup(jobinfo) ;
      if ~isfield(options, 'immediately_bg') || options.immediately_bg==0,
        if 1%PAR.RPROC.exm_per_batch==1,
          rproc_cleanup(jobinfo) ;
        else
          rproc_cleanup(meta_jobinfo) ;
        end ;
      end  ;
    end ;

    %fprintf(fid,'\n number of new constraints: %i\n',num_new_constraints );
    %if num_new_constraints==0, continue; end
    num_new_constraints_iteration = num_new_constraints_iteration+num_new_constraints;

    %% solve the new QP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    load(fn_parameters, 'decode_time', 'solve_time', 'solve_cnt', 'param_history')
    solve_cnt = solve_cnt+1;
    fprintf(fid, '\n solving qp: %i\n', solve_cnt);
    t = cputime;
    QP = solve_qp_iter(QP) ;
    %QP = solve_qp(QP, PAR) ;
    solve_time(solve_cnt) = cputime - t;
    decode_time(solve_cnt)=sum(time);
    [xis,parameters,hvar,H,obj] = res2param(QP.res, num_examples, model, QP, 0, PAR.weight_names);
    %obj
    param_history{end+1} = parameters;
    parameters_fields = fieldnames(parameters) ;
    save('-v7', fn_parameters, 'decode_time', 'solve_time', 'solve_cnt', 'param_history', 'parameters_fields')

    if exist('obj_solve')&&~isempty(fieldnames(obj_solve))%not an empty struct
      obj_solve(solve_cnt) = obj;
    else
      tmp11(solve_cnt) = obj;
      obj_solve = tmp11;
      clear tmp11;
    end
    print_struct(fid, obj);
    save('-v7', [PAR.FN.output_lsl.fn_train,'.mat'], 'PAR', 'QP','num_examples') ;
    %draw_pred(blocks_train(block_id).truth, blocks_train(block_id).pred,model) 
  end

  TRAIN_PAR = PAR;
  if ~exist('obj_solve'),
    obj_solve=struct ;
  end ;
  save('-v7', sprintf('%s_iteration%i.mat', PAR.FN.output_lsl.fn_train,iteration), 'TRAIN_PAR', 'QP','num_examples','obj_solve') ;
  if isfield(PAR, 'debug_mode') && PAR.debug_mode == 1
    save_struct(sprintf('%s_iteration%i_blocks.mat', PAR.FN.output_lsl.fn_train,iteration), blocks_train) ;
  end

  %% evaluation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if mod(iteration, 10)==0 && ~never_submit
    jobinfo_predict = eval_iter(PAR, iteration); 
  end ;

  %% termination
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if length(obj_solve)>20 && std([obj_solve(end-4:end).OBJ])<0.05
    fprintf('\n\n objective change very small %i\n\n', iteration);
	if ~never_submit
    	jobinfo_predict = eval_iter(PAR, iteration); 
	end
    if compute_full_gen_paths
      break;
    else
      fprintf('\ndo another iteration without approximation\n');
      for bb = 1:length(blocks_train)
        blocks_train(bb).pred_use_approximation = 0;
      end
    end
  end
  if num_new_constraints_iteration==0,
    fprintf('\n\nno new constraints during iteration %i', iteration);
    if compute_full_gen_paths
      fprintf('\nfinishing optimization.\n', iteration);
	  if ~never_submit
      	jobinfo_predict = eval_iter(PAR, iteration); 
	  end
      break;
    else
      fprintf('\ndo another iteration without approximation\n');
	  never_approx = 1;
    end
  elseif num_new_constraints_iteration<50
	never_approx = 1;
  end
  if never_approx
    for bb = 1:length(blocks_train)
      blocks_train(bb).pred_use_approximation = 0;
    end
  end
end
if fid>2
  fclose(fid) ;
end

fprintf('\nWaiting for predictions on validation/training set\n') ;

jobinfo_predict = rproc_wait(jobinfo_predict, 60, 1, 0, 2) ;

for i=1:iteration,
  fn_pred = [PAR.FN.output_lsl.fn_pred 'iteration_' num2str(i) '/results.mat'];
  if fexist(fn_pred),
    L=load(fn_pred) ;
    if isfield(L, 'res'),
      fprintf('\nResults on validation/training set from iteration %i:\n', i) ;
      fprintf('Exon level accuracy:\n')
      print_struct(1, L.res.cds_exons) ;
      fprintf('Transcript level accuracy:\n')
      print_struct(1, L.res.cds_transcripts) ;
      %if isequal(engine, 'octave'), struct_levels_to_print(2) ; end ;
      %L.pred
      %if isequal(engine, 'octave'), struct_levels_to_print(0) ; end ;
    end ;
  end ;
end ;
rproc_cleanup(jobinfo_predict) ;

fn_predictor = sprintf('%s_iteration%i.mat', PAR.FN.output_lsl.fn_train, iteration);

if isfield(PAR, 'save_fname') && ~isempty(PAR.save_fname),
  unix(sprintf('cp %s %s', fn_predictor, PAR.save_fname)) ;
  save_append(PAR.save_fname, 1, 'iter', iteration, 'fn_predictor', fn_predictor) ;
end ;


fprintf('\nClean up temporary files\n') ;
unix(sprintf('rm %s*', blocks_name));

fprintf('Trained gene predictor in %s\n\n', fn_predictor) ;

return

function jobinfo_predict = eval_iter(PAR, iteration)
  fprintf('\nevaluate iteration %i\n', iteration);
  myoptions = PAR.RPROC.options;
  myoptions.force_octave = 0 ;
  myoptions.express = 0;
  myoptions.start_dir = get_base_dir();
  myoptions.addpaths = {shogun_settings(), fileparts(which('predict'))};
  myoptions.verbosity = 2 ;
  PAR.fn_blocks = PAR.FN.input_lsl.fn_val_blocks; 
  PAR.iteration = iteration;
  PAR.fn_pred = [PAR.FN.output_lsl.fn_pred 'iteration_' num2str(iteration) '/'];
  PAR.evaluate = 1;
  [mem_req, time_req, opts] = rproc_memtime_policy('predict', 0, myoptions) ;
  myoptions.verbosity=2 ;
  %jobinfo_predict = rproc('predict', PAR, mem_req, opts, time_req);
  jobinfo_predict = rproc('predict', PAR, 20000, opts, time_req);

  for i=1:iteration,
    fn_pred = [PAR.FN.output_lsl.fn_pred 'iteration_' num2str(i) '/results.mat'];
    if fexist(fn_pred),
      L=load(fn_pred) ;
      if isfield(L, 'pred'),
        fprintf('\nResults on validation/training set from iteration %i:\n', i) ;
        fprintf('Exon level accuracy:\n')
        print_struct(1, L.res.cds_exons) ;
        fprintf('Transcript level accuracy:\n')
        print_struct(1, L.res.cds_transcripts) ;
      end ;
    end ;
  end ;

return

function blocks_train = remove_unused_fields(blocks_train)
	if isfield(blocks_train, 'tracks')
		blocks_train = rmfield(blocks_train, 'tracks');
	end
	if isfield(blocks_train, 'segment_lists')
		blocks_train = rmfield(blocks_train, 'segment_lists');
	end
	if isfield(blocks_train, 'segment_scores')
		blocks_train = rmfield(blocks_train, 'segment_scores');
	end
	if isfield(blocks_train, 'Signals')
		blocks_train = rmfield(blocks_train, 'Signals');
	end
	if isfield(blocks_train, 'features')
		blocks_train = rmfield(blocks_train, 'features');
	end
return
