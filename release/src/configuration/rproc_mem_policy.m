function [mem_mb, resubmit] = rproc_mem_policy(function_name, data_size, param) ;
% [mem_mb, resubmit] = rproc_mem_policy(function_name, data_size, repetition) ;

mem_mb = 7000 ;
resubmit = 0 ;

if isequal(function_name, 'load_and_train'),

elseif isequal(function_name, 'load_and_predict'),

elseif isequal(function_name, 'concat_genomewide_predictons_helper'),

elseif isequal(function_name, 'Out2Confhelper'),

elseif isequal(function_name, 'save_label_files_starter'),

elseif isequal(function_name, 'save_examples'),

end ;

fprintf('rproc_mem_policy for %s with data_size=%i and repetition=%i: %iMb\n', ...
	function_name, data_size, param, mem_mb) ;
