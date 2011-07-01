function [ret1, ret2]=signal_wrapper(P)
% [ret1, ret2] = signal_wrapper(P)

paths
%%%%%%%%%%%%%%%
try
  pause('on')
  warning('off', 'MATLAB:typeaheadBufferOverflow');  
catch
  % works only for matlab
end

ret1 = [] ;
ret2 = [] ;

signal_name = P.signal_name;  
config_dir = P.config_dir; 

if P.train_on 
  annotation_fname = P.annotation_fname;; 
  annotation_dir = P.annotation_dir; 

  spf_label_fname = P.spf_label_fname; 
  spf_label_dir = P.spf_label_dir; 
  covered_genic_positions_only = P.covered_genic_positions_only;
  
  tsp_train_info = P.tsp_train_info;
end

tsp_train_dir = P.tsp_train_dir;

spf_pred_fname = P.spf_pred_fname;
spf_pred_dir = P.spf_pred_dir;  

eval_fname = P.eval_fname;
use_old_files = P.use_old_files;

if isfield(P,'rproc_options')
  options = P.rproc_options
else
  options = [];
end

save_as.spf_ascii = 0 ;
save_as.spf_binary = 1 ;
save_as.wiggle = 0 ;

run_locally = P.run_locally;


%%%%%%%%%%%%%%%

if P.train_on
  % label generation
  if ~use_old_files
    unix(sprintf('rm -rf %s', spf_label_fname)) ;
  end
  anno2signallabel(annotation_fname, annotation_dir,config_dir,spf_label_fname, spf_label_dir, signal_name, covered_genic_positions_only ) ;

  % training 
  if ~use_old_files
    unix(sprintf('rm -rf %s', tsp_train_dir)) ;      
  end
  signal_train(spf_label_fname,spf_label_dir,  signal_name, config_dir,tsp_train_info, tsp_train_dir,run_locally, options) ;
end
   

% prediction
if ~use_old_files
  unix(sprintf('rm -rf %s',spf_pred_dir )) ;
end
signal_predict(tsp_train_dir, config_dir, spf_pred_fname, spf_pred_dir,run_locally, save_as,options) 

if 0% P.train_on
  %  evaluation
  unix(sprintf('rm -rf %s', eval_fname)) ;
  [ROC, PRC] = signal_eval(spf_label_fname,spf_label_dir,spf_pred_fname, spf_pred_dir, eval_fname) ;
end
disp('done') ;

% fd = fopen(eval_fname, 'r') ;
% text = fread(fd,'uint8=>char') ;
% fclose(fd) ;
% disp(text)
