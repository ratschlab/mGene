function jobinfo = ms_gfwd_starter(PAR)

% ms_gfwd_starter(PAR)
%
% used for modelselection as well as simple training and evaluation 
% INPUT : struct PAR as generated by create_PAR
%         fields used: 
%         - Signal_name: string specifying the signal that is actually computed
%           (e.g. 'acc','don') 
% 
%         - Signals.(Signal_name) : struct containing all necessary information to
%           compute the signal,
% 
%         - SETs : struct with information about split and partitions, 
%   
%         - FN : struct containing all used filenames including: 
%           FN.input_sig.(Signal_name).fn_filter_settings; 
%           FN.input_sig.(Signal_name).fn_examples; 
%           FN.output_sig.(Signal_name).fn_models;
%           FN.output_sig.(Signal_name).fn_SVMs;
%           FN.output_sig.(Signal_name).fn_motifs;
%   
%         - tasks.signals struct specifying what will be done
%           fields include run_locally (default=1), run_all,
%  
%         - RPROC (only necessary when tasks.signals.run_locally=0)
  
  
fid = 1;%fopen(PAR_1.FN.fn_log ,'a+');
FN = PAR.FN;
if ~isempty(PAR.Signal_name)
  Signal_name = PAR.Signal_name;
  Signal = PAR.Signals.(Signal_name);
  input_fn = FN.input_sig.(Signal_name);
  output_fn = FN.output_sig.(Signal_name);
  tasks = PAR.tasks.signals;
else
  Signal_name = PAR.Content_name;
  Signal = PAR.Contents.(Signal_name);
  input_fn = FN.input_cont.(Signal_name);
  output_fn = FN.output_cont.(Signal_name);
  tasks = PAR.tasks.contents;
end
fn_filter_settings = input_fn.fn_filter_settings;
fn_examples = input_fn.fn_examples;


method = Signal.method;
SETs = PAR.SETs;

if ~isfield(tasks,'run_locally')
  tasks.run_locally = 1;
end
if ~isfield(tasks,'run_all')
  tasks.run_all = 0;
end 

if ~tasks.run_locally.train
  [sg_path, envstr] = shogun_settings();
  RPROC = PAR.RPROC;
  organism = PAR.organism.name;
  RPROC.options.identifier = sprintf('%s%s',Signal_name,upper(organism)) ;
  RPROC.options.start_dir = get_base_dir();
  RPROC.options.addpaths = {sg_path, fileparts(which('load_and_train'))};
  RPROC.options.envstr = envstr;
  RPROC.options.ncpus = 1;
end  

fprintf(fid,'\nstarting ms_gfwd_starter \n');
fprintf(fid,'date: %s\n',date);
fprintf(fid,'\nSITE: %s\n',Signal_name);


%------------------------------------------------------------------------%
% write numeric kernelpars into field for model selection
%------------------------------------------------------------------------%

if isfield(method, 'kernels') && ~isempty(method.kernels{1})
  method.par_ms = collect_ms_parameters(method.kernels, method.par_ms);
end

%------------------------------------------------------------------------%
% initiate parameters
%------------------------------------------------------------------------%
Model_param = getfield(method,'par_ms');
Model_param = orderfields(Model_param) ;

[MODELS,param_names] = generate_all_models(orderfields(method.par_ms));
if ~tasks.modelsel
  fprintf('no model selection. taking first model;\n')
  MODELS = MODELS(:,1);
end

fprintf(fid,'number of models to be tested %s: %i\n',Signal_name,size(MODELS,2));
num_models_tested = size(MODELS,2);
RSV_eval = zeros(SETs.num_partitions,size(MODELS,2));
RFV_eval = zeros(SETs.num_partitions,size(MODELS,2));

%------------------------------------------------------------------------%
% check for saved model files and save models if non are found
%------------------------------------------------------------------------%
if ~fexist(output_fn.fn_models) && ~fexist([output_fn.fn_models '.mat']) 
  save('-V7', output_fn.fn_models,'MODELS','Model_param')
elseif (fexist(output_fn.fn_models) || fexist([output_fn.fn_models '.mat']) ) 
  LL = load(output_fn.fn_models);
  if isfield(LL, 'MODELS') && isfield(LL, 'Model_param')
    assert(isequal(LL.Model_param,Model_param))
    if ~isequal(LL.MODELS,MODELS)
      MODELS=append_new_models(MODELS, LL.MODELS);
      old=LL;
      save('-V7',output_fn.fn_models,'MODELS','old')
      warning('different other MODELS have already been trained, new models are appended')
    end
    clear old
  else
    save('-V7',output_fn.fn_models,'MODELS','Model_param')
  end
end
num_not_eval = SETs.num_partitions*sum(num_models_tested) ;

fprintf(fid,'total number of jobs for training: %i\n',sum(SETs.partition_train_on)*sum(num_models_tested));
fprintf(fid,'total number of jobs for evaluation: %i\n',num_not_eval);


P.SETs = SETs;
P.Signal = Signal;
P.fn_filter_settings = fn_filter_settings;
P.fn_examples = fn_examples;
P.train_error_on = tasks.get_train_error;
P.pred_from_seq = 0;

P.use_preproc = 1;

num_crashed = 1;
num_crashed_1 = 0;
rep = 1;

P.check_overlap = 1;
while num_crashed + num_crashed_1 + num_not_eval>0 & rep<5
  q = 0 ; nq = 0 ;
  jobinfo = rproc_empty(0) ;
  jobinfo_1 = rproc_empty(0) ;
  num_not_eval = SETs.num_partitions*sum(num_models_tested) ;
  for p = 1:size(SETs.partitions,1) %% loop over partitions  
    fprintf(fid,'\n\n-----------\n');
    fprintf(fid,'PARTITION %i\n',p);
    fprintf(fid,'-----------\n');
    
    if tasks.run_all
      Train_data=[];
      Test_data=[];
    end
    
    P.partition = p;
    P.train_on = SETs.partition_train_on(p);
    P.SETs.which_train = find(ismember(SETs.partitions(p,:),SETs.info.train));
    P.SETs.which_val = find(ismember(SETs.partitions(p,:),SETs.info.eval));
    P.SETs.which_pred = find(ismember(SETs.partitions(p,:),SETs.info.pred));
    
    fn_save = [sprintf('%sms_set=[',output_fn.fn_SVMs),...
               sprintf('%i',P.SETs.which_train)];
    fn_save_pred  = [fn_save '_' sprintf('%i',P.SETs.which_val)];
    
    for num_models = 1:size(MODELS,2)
      %--------------------------------------------------------------%
      % assign parameters of model number num_models to the current parameter
      % structure
      %--------------------------------------------------------------%
      
      P.Signal.method = assign_model(method, MODELS, param_names,num_models); 
      
      %--------------------------------------------------------------%
      %% generate filename for SVMs and predictions
      %--------------------------------------------------------------%
      
      P.fn_SVMs = sprintf('%s]_model=%i.mat',fn_save,num_models) ;
      P.fn_eval = sprintf('%s]_model=%i.mat',fn_save_pred,num_models) ;
      
      %--------------------------------------------------------------%
      %%%%% TRAIN
      %--------------------------------------------------------------%
      
      %if ~fexist(P.fn_SVMs) & P.train_on
      %  fprintf(fid,'waiting for file %s...\n', P.fn_SVMs);
      %pause(1)
      %end      

      if ~fexist(P.fn_SVMs) & P.train_on
        q = q+1 ;
        fprintf(fid,'\n start training on set ');
        fprintf(fid,'%i',P.SETs.which_train);
        fprintf(fid,'\n')
        
        train_fct = 'load_and_train';
        
        if tasks.run_all
          Data = load_data(Signal, fn_filter_settings, fn_examples, P.SETs.which_train, 'train', SETs.train_subsample_neg, SETs.train_subsample_pos);
          Train_data.(sprintf('partition_%i',p))= Data ; clear Data 
          feval(train_fct,P,Train_data.(sprintf('partition_%i',p)));
        elseif tasks.run_locally.train==0,
          fprintf('creating job (training) %i/%i  \r', q,sum(SETs.partition_train_on)*sum(num_models_tested)) ;
          RPROC.options.identifier = sprintf('Tr%i_%s%s',num_models, Signal_name, upper(organism)) ;

          num_data = count_data(Signal, fn_filter_settings, fn_examples, P.SETs.which_train, 'train', SETs.train_subsample_neg, SETs.train_subsample_pos);

          [mem_req, time_req, opts] = rproc_memtime_policy(train_fct,  num_data, RPROC.options) ;

          jobinfo(q) = rproc_create(train_fct, P, mem_req, opts, time_req) ;
        else
          feval(train_fct,P);
        end 
      	
      %--------------------------------------------------------------%
      %%%%% EVAL 
      %--------------------------------------------------------------%
      elseif fexist(P.fn_SVMs) 
        fprintf(fid,'\nSVM already trained on set ');
        fprintf(fid,'%i',P.SETs.which_train);
	fprintf(fid,'\n')
       
	if ~fexist(P.fn_eval)
          fprintf(fid,'waiting for file %s...\n', P.fn_eval);
	  pause(1)
        end
 
        if ~fexist(P.fn_eval) 
          P.eval_on = 1;
          fprintf(fid,'\n\nstart evaluation on set ');
          fprintf(fid,'%i',P.SETs.which_val);
          fprintf(fid,'\n')
          nq = nq+1;
          
          if isfield(method,'pred_fct')&& strcmp(method.pred_fct,'load_content_and_predict')
            pred_fct = 'load_content_and_predict';
          else
            pred_fct = 'load_and_predict';
          end
  
          if tasks.run_all
            Data = load_data(Signal,fn_filter_settings,fn_examples, P.SETs.which_val, 'test',SETs.test_subsample_neg,SETs.test_subsample_pos);
            Test_data.(sprintf('partition_%i',p)) = Data; clear Data;
            feval(pred_fct,P,Test_data.(sprintf('partition_%i',p)));
          elseif tasks.run_locally.eval==0,
            fprintf('creating job (evaluation) %i/%i  \r', nq,num_not_eval) ;
            RPROC.options.identifier = sprintf('Ev%i_%s%s',num_models,Signal_name,upper(organism)) ;

            num_data = count_data(Signal,fn_filter_settings,fn_examples, P.SETs.which_val, 'test', SETs.test_subsample_neg, SETs.test_subsample_pos);

            [mem_req, time_req, opts] = rproc_memtime_policy(pred_fct,  num_data, RPROC.options) ;

            jobinfo_1(nq) = rproc_create(pred_fct, P, mem_req, opts, time_req) ;
          else
            feval(pred_fct,P);
          end
        
        %--------------------------------------------------------------%
        %%%%% Collext Results 
        %--------------------------------------------------------------%
        elseif fexist(P.fn_eval) 
          fprintf(fid,'\n model already evaluated on set ');
          fprintf(fid,'%i',P.SETs.which_val);
          
          L = load(P.fn_eval, 'RSV', 'RFV');%,'fn_SVMs') ;
          num_not_eval = num_not_eval-1;
              
          assert(isnumeric(L.RFV) && prod(size(L.RFV))==1) ;
          assert(isnumeric(L.RSV) && prod(size(L.RSV))==1) ;
          RFV_eval(p,num_models)= L.RFV ;
          RSV_eval(p,num_models)= L.RSV ;

        else
          fprintf(fid,'\n only training, no evaluation ');
          num_not_eval = num_not_eval-1;
        end
      end   %% if trained eval
    end   %% loop over num_models      
    
    fprintf(fid,'partition %i: DONE\n', P.partition);
  end  %% loop over partitions
  fprintf(fid,'experiment %s: DONE\n',Signal_name );
      
  
  %--------------------------------------------------------------%
  %%%%% RPROC STUFF
  %--------------------------------------------------------------%

  %if (~isempty(jobinfo) || ~isempty(jobinfo_1)) && isfield(RPROC, 'collect_jobs') && RPROC.collect_jobs,
  %  jobinfo = [jobinfo jobinfo_1] ;
  %  return ;
  %end ;

  if ~isempty(jobinfo),
    for ii=1:length(jobinfo), jobinfo(ii) = rproc_resubmit(jobinfo(ii)); end ;
  else
    num_crashed = 0 ;
  end
  if ~isempty(jobinfo_1),
    %[jobinfo_1, meta_jobinfo_1] = rproc_submit_batch(jobinfo_1, RPROC.exm_per_batch) ;
    for ii=1:length(jobinfo_1), jobinfo_1(ii) = rproc_resubmit(jobinfo_1(ii)); end ;
  else
    num_crashed_1 = 0;
  end ;
  if ~isempty(jobinfo),
    [jobinfo,num_crashed] = rproc_wait(jobinfo, 60, 1, -1);
    if num_crashed==length(jobinfo)
      fprintf('all jobs seem to be crashed. Waiting 2min.\n') ;
      pause(120) ;
      [jobinfo,num_crashed] = rproc_wait(jobinfo, 60, 1, -1);
      if num_crashed==length(jobinfo),
        error('all jobs crashed')
      end ;
    end
    if num_crashed == 0,
      rproc_cleanup(jobinfo) ;
    end ;
  end
  if ~isempty(jobinfo_1),
    [jobinfo_1,num_crashed_1] = rproc_wait(jobinfo_1, 60, 1, -1);
    if num_crashed_1==length(jobinfo_1)
      fprintf('all jobs seem to be crashed. Waiting 2min.\n') ;
      pause(120) ;
      [jobinfo,num_crashed] = rproc_wait(jobinfo, 60, 1,-1);
      if num_crashed==length(jobinfo),
        error('all jobs crashed') ;
      end ;
    end
    if num_crashed_1 == 0,
      rproc_cleanup(jobinfo_1) ;
    end ;
  end ;
  
  %--------------------------------------------------------------%
  %%%%% FINAL STUFF
  %--------------------------------------------------------------%
  
  try
    L=load(output_fn.fn_models, 'RSV_eval', 'RFV_eval') ;
  catch
    L=struct ;
  end ;
  if ~isfield(L, 'RSV_eval') || ~isfield(L, 'RFV_eval') || ~isequal(L.RSV_eval, RSV_eval) || ~isequal(L.RFV_eval, RFV_eval),
    %save(output_fn.fn_models,'-append','RSV_eval','RFV_eval')
    save_append(output_fn.fn_models, 1,'RSV_eval',RSV_eval, 'RFV_eval',RFV_eval)
    % save_append(output_fn.fn_models, 1,'RSV_eval',RSV_eval, 'RFV_eval',RFV_eval)
  end 
  if num_crashed+num_crashed_1>0
    rep = rep+1;
  end
end  %% loop over crashed jobs

%fclose(fid)



