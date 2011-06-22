function [run_locally] = rproc_policy(task, genome_info, param)
% run_locally = rproc_policy(task, genome_info)
% determines based on genome properties, whether the task to be solved 
%   should use subjobs or should run sequentially

if ~isempty(genome_info),
  num_contigs = length(genome_info.contig_names) ;
  total_size = 0 ;
else
  num_contigs = 0 ;
  total_size = 0 ;
end ;

if isequal(task, 'gff2anno') 
  % param = number of genes
  if num_contigs<1000 || num_contigs>25 || num_contigs<=2,
    run_locally = 1 ;
  else
    run_locally = 1 ;
  end ;
elseif isequal(task, 'signal_wrapper') || isequal(task, 'content_wrapper'),
  run_locally = 1 ;
elseif isequal(task, 'signal_train:pos') || isequal(task, 'content_train:pos'),
  run_locally = 1 ;
elseif isequal(task, 'signal_train:train') || isequal(task, 'content_train:train'),
  run_locally = 1 ;
elseif isequal(task, 'signal_train:eval') || isequal(task, 'content_train:eval'),
  run_locally = 1 ;
elseif isequal(task, 'signal_train:conf') || isequal(task, 'content_train:conf'),
  run_locally = 1 ;
elseif isequal(task, 'signal_train:predict') || isequal(task, 'content_train:predict'),
  run_locally = 1 ;
elseif isequal(task, 'signal_train:save') || isequal(task, 'content_train:save'),
  run_locally = 1 ;
elseif isequal(task, 'signal_predict:pos') || isequal(task, 'content_predict:pos')
  run_locally = 1 ;
elseif isequal(task, 'signal_predict:predict') || isequal(task, 'content_predict:predict')
  if num_contigs>=3 && num_contigs<=30,
    run_locally = 1 ;
  else
    run_locally = 1 ;
  end ;
elseif isequal(task, 'signal_predict:conf') || isequal(task, 'content_predict:conf')
  run_locally = 1 ;
elseif isequal(task, 'signal_predict:concat') || isequal(task, 'content_predict:concat')
  run_locally = 1 ;
elseif isequal(task, 'signal_predict:save') || isequal(task, 'content_predict:save')
  run_locally = 1 ;
elseif isequal(task, 'gene_train') 
  run_locally=1 ;  
elseif isequal(task, 'gene_predict') 
  if num_contigs>=3 && num_contigs<=300,
    run_locally = 1 ;
  else
    run_locally = 1 ;
  end ;
elseif isequal(task, 'mgene_train') 
  run_locally = 1 ;  
elseif isequal(task, 'mgene_predict') 
  run_locally = 1 ;    
elseif isequal(task, 'train_path:gen_paths:noapprox') 
  if param>3000,
    run_locally = 1 ;    
  else
    run_locally = 1 ;    
  end ;
elseif isequal(task, 'train_path:gen_paths:approx')
  run_locally= 1; 
elseif isequal(task, 'pred_path:gen_paths') 
  if param == 1 % this is the case for genome wide predictions
    run_locally = 1;
  else	
    run_locally = 1;
  end
else
  fprintf('did not find rproc policy for task %s\n', task) ;  
  run_locally = 1 ;
end ;
