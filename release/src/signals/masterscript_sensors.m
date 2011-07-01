function [PAR, jobinfo ,all_done]=masterscript_sensors(PAR,check_files,keep_all_files,verb)
% [PAR, jobinfo ,all_done] = masterscript_sensors(PAR,check_files,verb)

%cd ~/svn/projects/genefinding/sensors
%paths

tasks = PAR.tasks.signals;  
%run_locally = tasks.run_locally;

FN = PAR.FN;  
fn_genome_config = PAR.FN.input.fn_genome_config;
 
if nargin<2 
  check_files =0;
end
if nargin<3 
  keep_all_files = 0;
end
if nargin<4 
  verb = 0;
end

fid =1;
fprintf(fid,'\n\n-----------\n');
fprintf(fid,'SENSOR :  %s\n',PAR.Signal_name);
fprintf(fid,'-----------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Train SENSORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAR.method.submit_batch=1;
all_done=0;
%train_ok = 0;

jobinfo=struct ;

if tasks.prepare_labels
  %fprintf('start training for sensor %s on organism %s with method %s: %s\n\n',sensor,organism,method)
  gen_training_examples4signal(PAR) ;
end
if tasks.modelsel | tasks.train
  jobinfo = ms_gfwd_starter(PAR);
  if ~isempty(jobinfo), return ; end ;
end
if tasks.get_train_error
  jobinfo = predict_testplif_geterror(PAR)
  if ~isempty(jobinfo), return ; end ;
end
if tasks.learn_plifs
  calc_plifs(PAR);
end


if ~keep_all_files
  unix(sprintf('rm -rf %s', FN.input_sig.(PAR.Signal_name).fn_candsites));
  unix(sprintf('rm -rf %s*', FN.input_sig.(PAR.Signal_name).fn_examples));
  unix(sprintf('rm -rf %s*', FN.input_sig.(PAR.Signal_name).fn_example_statistics(1:end-4)));  
  unix(sprintf('rm -rf %s', FN.input_sig.(PAR.Signal_name).fn_pos));
  unix(sprintf('rm -rf %s', FN.input_sig.(PAR.Signal_name).fn_filter_settings));

  unix(sprintf('rm -rf %s', FN.output_sig.(PAR.Signal_name).fn_pred));
  unix(sprintf('rm -rf %s/ms_set*', FN.output_sig.(PAR.Signal_name).fn_SVMs));
  unix(sprintf('rm -rf %s/partition=*', FN.output_sig.(PAR.Signal_name).fn_SVMs));  
end

if ~tasks.pred_on_genome
  all_done = 1 ;%train_ok ;
  return
end
