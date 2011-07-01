function ret = predict_testplif_geterror(PAR)

  
SETs = PAR.SETs;
RPROC = PAR.RPROC;
FN = PAR.FN;

if ~isempty(PAR.Signal_name)
  Signal_name = PAR.Signal_name;
  Signal = PAR.Signals.(Signal_name);
  output_fn = FN.output_sig.(Signal_name);
  tasks = PAR.tasks.signals;
else
  Signal_name = PAR.Content_name;
  Signal = PAR.Contents.(Signal_name);
  output_fn = FN.output_cont.(Signal_name);
  tasks = PAR.tasks.contents;
end


fid = 1;%fopen(PAR.FN.fn_log ,'a+');
fprintf(fid,'\n\n starting test on left out set \n');
fprintf(fid,'date: %s\n',date);


%PAR = setfield(PAR,PAR.method,Model_param);


load(output_fn.fn_models,'RSV_eval','RFV_eval','MODELS','Model_param')
param_names = fieldnames(Model_param); 

fprintf(fid,'\n\n-----------\n');

RSV_test = zeros(SETs.num_partitions,1);
RFV_test = zeros(SETs.num_partitions,1);
best_models = zeros(SETs.num_partitions,size(MODELS,1));

[sg_path, envstr] = shogun_settings();
options = RPROC.options;
organism = PAR.organism.name;
options.start_dir = get_base_dir();
options.identifier = [ Signal_name '_' Signal.method.name];
options.addpaths = {sg_path, fileparts(which('load_and_predict'))};
options.envstr = envstr;


jobinfo = rproc_empty(0) ;
jobinfo_1 = rproc_empty(0) ;  


for p = 1:size(SETs.partitions,1) %% loop over partitions  
  fprintf(fid,'\n\n-----------\n');
  fprintf(fid,'PARTITION %i\n',p);
  fprintf(fid,'-----------\n');
  fprintf(fid,'\nSITE: %s\n',Signal_name);
  
  which_train = find(ismember(SETs.partitions(p,:),SETs.info.train));
  which_val = find(ismember(SETs.partitions(p,:),SETs.info.eval));
  which_pred = find(ismember(SETs.partitions(p,:),SETs.info.pred));
  which_plif = find(ismember(SETs.partitions(p,:),SETs.info.plif));

  fn_best_SVM = sprintf('%s_best_partition=%i.mat', output_fn.fn_SVMs,p) ;
  fn_SVM_train = [sprintf('%sms_set=[',output_fn.fn_SVMs),...
                  sprintf('%i',which_train)];
  fn_SVM_eval  = [fn_SVM_train '_' sprintf('%i',which_val)];
  fn_pred = sprintf('%spartition=%i.mat',output_fn.fn_SVMs,p);
  fn_plifs = sprintf('%spartition=%i_plifs.mat',output_fn.fn_SVMs,p);
  
  %--------------------------------------------------------------%
  %%% select for each partition best model:
  %--------------------------------------------------------------%
  if ~all(RFV_eval(p,:)~=0) 
    warning('some models not evaluated')
  end
  [tmp,idx] = max(RFV_eval(p,:));
  fn_SVM_eval = sprintf('%s]_model=%i.mat',fn_SVM_eval,idx);
  fn_SVM_train = sprintf('%s]_model=%i.mat',fn_SVM_train,idx)
  load(fn_SVM_train,'Train');
  assert(isequal(Train.PAR.SETs.which_train,which_train))
  best_models(p,:) = MODELS(:,idx) ;
  for ii = 1:length(param_names)
    assert(isequal(getfield(Train.PAR.Signal.method.par_ms,param_names{ii}), MODELS(ii,idx)))
  end
  
  fprintf(fid,'\nbest model (%s,partition=%i): ',Signal_name,p);
  fprintf('\nbest model (%s,partition=%i): \n',Signal_name,p);
  for num_param = 1:size(MODELS,1)
    fprintf(1,'%s = %d; \t',param_names{num_param},getfield(Train.PAR.Signal.method.par_ms,param_names{num_param}));
    fprintf(fid,'%s = %d; \t',param_names{num_param},getfield(Train.PAR.Signal.method.par_ms,param_names{num_param}));
  end
  L = load(fn_SVM_eval) ;
  if ~nan(L.RFV)
    assert(L.RFV == RFV_eval(p,idx)) ;
    assert(L.RSV == RSV_eval(p,idx)) ;
  end
  fprintf(fid,'\nvalidation error: RFV:%.4f \t RSV:%.4f\n',L.RFV,L.RSV);
  %fprintf('\nvalidation error: RFV:%.4f \t RSV:%.4f\n',L.RFV,L.RSV);
  
  %--------------------------------------------------------------%
  %%%%% Predict on TEST SET
  %--------------------------------------------------------------%  
  
  Train.PAR.partition = p;
  Train.PAR.SETs.which_val = which_val;
  Train.PAR.SETs.which_pred = which_pred ;
  Train.PAR.SETs.which_plif = which_plif ;
  Train.PAR.fn_plifs = fn_plifs;
  Train.PAR.fn_SVM_eval = fn_SVM_eval;
  Train.PAR.fn_pred = fn_pred;
  
  PAR_new = Train.PAR;
  PAR_new.eval_on = 0;
  PAR_new.Signal.filter_label.test = Signal.filter_label.test;

  if ~fexist(PAR_new.fn_pred),
    if ~isfield(tasks, 'run_locally') || tasks.run_locally.predict==0,

      [mem_req, time_req, opts] = rproc_memtime_policy('load_and_predict',  0, options) ;%RPROC.MEMREQ ;

      jobinfo(end+1) = rproc('load_and_predict',PAR_new, mem_req, opts, time_req);
      % load_and_predict(PAR_new)
    else
      %D=load_data(PAR_new,which_pred,'test',PAR.SETs.test_subsample_neg,PAR.SETs.test_subsample_pos); 
      load_and_predict(PAR_new);%, D);
      %clear D
    end
  end ;
  %--------------------------------------------------------------%
  %%%%% Predict on PLIF SET
  %--------------------------------------------------------------%  
  
  PAR_new.eval_on = 0;
  PAR_new.plif_on = 1;
  
  
  
  if ~fexist(PAR_new.fn_plifs)
    if ~isfield(tasks, 'run_locally') || tasks.run_locally.predict==0,
      [mem_req, time_req, opts] = rproc_memtime_policy('load_and_predict',  0, options) ;

      jobinfo_1(end+1) = rproc('load_and_predict', PAR_new, mem_req, opts, time_req);
      % load_and_predict(PAR_new)
    else
      %D=load_data(PAR,which_plif,'test',PAR.SETs.test_subsample_neg,PAR.SETs.test_subsample_pos); 
      load_and_predict(PAR_new);%, D);
      %clear D
    end
  end ;
end  %% loop over partitions
if ~isempty(jobinfo),
  [jobinfo,num_crashed] = rproc_wait(jobinfo, 20, 1, -1) ;
end ;
if ~isempty(jobinfo_1),
  [jobinfo_1,num_crashed] = rproc_wait(jobinfo_1, 20, 1, -1) ;
end ;

if ~isempty(jobinfo)||~isempty(jobinfo_1),
  %unix('sync')
  pause(30)
  %unix('sync')
end ;

for p = 1:size(SETs.partitions,1) %% loop over partitions 
  fprintf(fid,'\n\n-----------\n');
  fprintf(fid,'PARTITION %i\n',p);
  fprintf(fid,'-----------\n');
  fprintf(fid,'\nSITE: %s\n',Signal_name);
    
  fn_best_SVM = sprintf('%sbest_partition=%i.mat', output_fn.fn_SVMs,p) ;
  fn_pred = sprintf('%spartition=%i.mat',output_fn.fn_SVMs,p);
  L = load(fn_pred,'RSV','RFV','svmFileName','PAR') ;
  load(L.svmFileName,'Train')
  Train.PAR = L.PAR;

  if ~fexist(fn_best_SVM)
    save_append(fn_best_SVM, 0, L)
    save_append(fn_best_SVM, 2, 'Train', Train)
  else
    fprintf('file %s already exists\n',fn_best_SVM),
  end
  fprintf(fid,'\ntest error: RFV:%.4f \t RSV:%.4f\n',L.RFV,L.RSV);
  fprintf('\ntest error: RFV:%.4f \t RSV:%.4f\n',L.RFV,L.RSV);
  
  RFV_test(p)= L.RFV ;
  RSV_test(p)= L.RSV ;
  
  clear PAR_new
end  %% loop over partitions

L=load(output_fn.fn_models, 'best_models','Model_param','RSV_test','RFV_test') ;
if ~isfield(L, 'best_models') || ~isfield(L, 'Model_param') || ~isfield(L, 'RSV_test') || ~isfield(L, 'RFV_test') || ...
      ~isequal(L.best_models, best_models) || ~isequal(L.Model_param, Model_param) || ~isequal(L.RSV_test, RSV_test) || ~isequal(L.RFV_test, RFV_test)
  %save(output_fn.fn_models,'-append','best_models','Model_param','RSV_test','RFV_test')
  save_append(output_fn.fn_models, 1, 'best_models', best_models, ...
              'Model_param', Model_param, 'RSV_test', RSV_test, ...
              'RFV_test', RFV_test)
end ;

fprintf(fid,'\n%s\n',Signal_name);
fprintf(fid,'RSV (mean and std over %i partitions): %.4f +/- %.4f\n',length(RSV_test),mean(RSV_test),std(RSV_test));
fprintf(fid,'RFV (mean and std over %i partitions): %.4f +/- %.4f\n',length(RSV_test),mean(RFV_test),std(RFV_test));
ret = [];
%fclose(fid)




