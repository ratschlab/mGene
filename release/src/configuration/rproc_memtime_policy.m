function [mem_mb, time_min, opts] = rproc_memtime_policy(function_name, data_size, options, verbosity) ;

if nargin<4,
  verbosity=1 ;
end ;
if nargin<3,
  options=struct ;
end ;
  
mem_mb = 4000 ;
time_min = 2000 ;

resubmit = 2 ;
opts = options ;
opts.resubmit = resubmit ;
opts.mem_req_resubmit = [8000 12000] ;
opts.time_req_resubmit = [4000 8000] ;
opts.data_size = data_size ;
opts.hard_time_limit = 720; % in minutes
opts.waitonfull = 1; % wait if to many jobs are in the queue
opts.maxjobs = 2990;


if isequal(function_name, 'gen_all_positions_rproc'),
    mem_mb = 10000 ;
    time_min = 4000 ;
    opts.resubmit = 1 ;
    opts.mem_req_resubmit = [20000] ;
    opts.time_req_resubmit = [8000] ;

elseif isequal(function_name, 'load_and_train'),
    mem_mb = 10000 ;
    time_min = 4000 ;
    opts.resubmit = 1 ;
    opts.mem_req_resubmit = [30000] ;
    opts.time_req_resubmit = [4000] ;
	opts.hard_time_limit = 1440; % in minutes

elseif length(function_name)>=length('load_and_predict') && isequal(function_name(1:length('load_and_predict')), 'load_and_predict')
    mem_mb = 20000 ;
    time_min = 4000 ;
    opts.resubmit = 1 ;
    opts.mem_req_resubmit = [40000] ;
    opts.time_req_resubmit = [8000] ;

	if isequal(function_name, 'load_and_predict:cleave')||isequal(function_name, 'load_and_predict:tss'),
	    %mem_mb = 90000 ;%maxvmem=84.656G
	    mem_mb = 30000 ;
    	time_min = 4000 ;
    	opts.resubmit = 1 ;
    	opts.mem_req_resubmit = [50000] ;
    	opts.time_req_resubmit = [8000] ;
	end

elseif isequal(function_name, 'concat_genomewide_predictons_helper'),
	opts.hard_time_limit = 60; % in minutes
    mem_mb = 8000 ;
    time_min = 4000 ;
    opts.resubmit = 1 ;
    opts.mem_req_resubmit = [30000] ;
    opts.time_req_resubmit = [8000] ;

elseif isequal(function_name, 'Out2Confhelper'),
	opts.hard_time_limit = 60; % in minutes

elseif isequal(function_name, 'save_label_files_starter'),
	opts.hard_time_limit = 60; % in minutes

elseif isequal(function_name, 'save_examples'),
	opts.hard_time_limit = 120; % in minutes

elseif length(function_name)>=length('signal_wrapper') && isequal(function_name(1:length('signal_wrapper')), 'signal_wrapper')
	opts.hard_time_limit = 2880; % in minutes
    mem_mb = 10000 ;
    time_min = 4000 ;
    opts.resubmit = 1 ;
    opts.mem_req_resubmit = [15000] ;
    opts.time_req_resubmit = [8000] ;

  if isequal(function_name, 'signal_wrapper:tss') || isequal(function_name, 'signal_wrapper:cleave')
    mem_mb = 15000 ;
    time_min = 4000 ;
    opts.resubmit = 1 ;
    opts.mem_req_resubmit = [20000] ;
    opts.time_req_resubmit = [8000] ;
  end ;

elseif isequal(function_name, 'gen_paths'),

  mem_mb = 1500 ;
  time_min = 100 ;
  
  if abs(data_size)>1e7
    mem_mb = mem_mb*10 ;
  elseif abs(data_size)>1e6,
    mem_mb = mem_mb*5 ;
  elseif abs(data_size)>5e5,
    mem_mb = mem_mb*3 ;
  end ;
  
  resubmit = 2 ;
  opts = options ;
  opts.resubmit = resubmit ;
  opts.mem_req_resubmit = [mem_mb*2 mem_mb*4] ;
  opts.time_req_resubmit = [500 2000] ;
  opts.data_size = data_size ;
  opts.hard_time_limit = 20; % in minutes
  %opts.hard_time_limit = 40; % in minutes
  
elseif isequal(function_name, 'predict'),
  
  mem_mb = 8000 ;
  time_min = 1000 ;
  
  resubmit = 0 ;
  opts = options ;
  opts.resubmit = resubmit ;
  opts.mem_req_resubmit = [] ;
  opts.time_req_resubmit = [] ;
  opts.data_size = data_size ;

end ;


if verbosity>=1
  if ~isempty(data_size),
    fprintf('rproc_memtime_policy for %s with data_size=%i: %i Mb and %i minutes\n', ...
            function_name, data_size(1), mem_mb, time_min) ;
  else
    fprintf('rproc_memtime_policy for %s: %i Mb and %i minutes\n', ...
            function_name, mem_mb, time_min) ;
  end ;
end ;

