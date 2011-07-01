function [ret1,ret2] = content_wrapper(P)
% [ret1,ret2]=content_wrapper(P)  

paths
try
  pause('on')
  warning('off', 'MATLAB:typeaheadBufferOverflow');
catch
  % works only for matlab
end

ret1 = [] ;
ret2 = [] ;

%%%%%%%%%%%%%%%

content_name = P.content_name;  
config_dir = P.config_dir; 

if P.train_on
  annotation_fname = P.annotation_fname;; 
  annotation_dir = P.annotation_dir; 

  cpf_label_fname = P.cpf_label_fname; 
  cpf_label_dir = P.cpf_label_dir; 
  covered_genic_positions_only = P.covered_genic_positions_only;
  vs_content_names = P. vs_content_names;
  intergenic_dist_to_genes = P.intergenic_dist_to_genes;

  tcp_train_info = P.tcp_train_info;
end
tcp_train_dir = P.tcp_train_dir;

cpf_pred_fname = P.cpf_pred_fname;
cpf_pred_dir = P.cpf_pred_dir;  

eval_fname = P.eval_fname;

use_old_files = P.use_old_files;
run_locally = P.run_locally;;
if isfield(P,'rproc_options')
  options = P.rproc_options;
else
  options = [];
end

save_as.spf_ascii = 0 ;
save_as.spf_binary = 1 ;
save_as.wiggle = 0 ;

run_locally = P.run_locally;


%%%%%%%%%%%%%%%

% [engine,env]=determine_engine;

if P.train_on
  % label generation
  if ~use_old_files
    unix(sprintf('rm -rf %s', cpf_label_fname)) ;
  end
  anno2contentlabel(annotation_fname, annotation_dir, config_dir, cpf_label_fname, cpf_label_dir, content_name, vs_content_names, intergenic_dist_to_genes) ;
  
  % training
  if ~use_old_files
    unix(sprintf('rm -rf %s', tcp_train_dir)) ;      
  end
  content_train(cpf_label_fname, cpf_label_dir,  content_name, config_dir, annotation_dir,tcp_train_info, tcp_train_dir, run_locally,options) ;
  
end

% keyboard

% prediction
if ~use_old_files
  unix(sprintf('rm -rf %s',cpf_pred_dir )) ;
end
content_predict(tcp_train_dir, config_dir, cpf_pred_fname, cpf_pred_dir, run_locally, save_as, options);

% evaluation
if P.train_on
unix(sprintf('rm -rf %s', eval_fname)) ;
content_eval(cpf_label_fname, cpf_label_dir, cpf_pred_fname, cpf_pred_dir, eval_fname) ;

end

disp('done') ;

% fd = fopen(eval_fname, 'r') ;
% text = fread(fd,'uint8=>char') ;
% fclose(fd) ;
% disp(text)
