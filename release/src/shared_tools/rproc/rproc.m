function [jobinfo]=rproc(ProcName, P1, Mem, options, time)
% [jobinfo]=rproc(ProcName, P1, Mem, options, time)
%
% time in minutes
% mem in mb

% currently only has adversal effects, change this when a new
% gridengine version has been installed
time = 10000 ; 

[engine, environment] = determine_engine() ;

if ~isfield(options, 'force_octave'),
  options.force_octave = 0 ;
end ;
if ~isfield(options, 'force_matlab'),
  options.force_matlab = 0 ;
end ;

if (isequal(engine, 'matlab') && options.force_octave==0) || options.force_matlab==1,
  use_engine='matlab' ;
else
  use_engine='octave' ;
end ;


if nargin<3, Mem=300 ; end ;


if isequal(use_engine, 'matlab'),
  if isequal(environment, 'internal'),
    Matlab_licenses = 1 ; 
  else
    Matlab_licenses = 0 ; 
  end ;

  if Mem<200, warning('You specified to allocate less than 200Mb for your job. Matlab might not be able to start.') ; end ;
else
  Matlab_licenses = 0 ; 
  if Mem<300, warning('You specified to allocate less than 300Mb for your job. Octave might not be able to start.') ; end ;
end ;

if isequal(environment, 'internal'),
    home_str = '/fml/ag-raetsch/home/';
elseif isequal(environment, 'galaxy'),
    home_str = '/home/' ;
end

use_reservation = 0 ;
if isfield(options,'ncpus') && ~isempty(options.ncpus) && ~equal(options.ncpus, 1)
  myassert(options.ncpus>1) ;
  Mem=Mem/options.ncpus ;
  Matlab_licenses = Matlab_licenses/options.ncpus ;
  use_reservation = 1 ;
end ;
if Mem>20000, use_reservation = 1 ; end ;

if nargin<4, options=[] ; end ;
if nargin<5, time=100*24*60 ; end ;
if Mem<100, Mem=100 ; end ;
if ~isfield(options,'verbosity')
  options.verbosity = 1 ;
end ;
if ~isfield(options,'maxjobs')
  options.maxjobs=3000; %hardcoded limit by sebastian
end
if ~isfield(options,'waitonfull')
  options.waitonfull=1 ;
end ;
if ~isfield(options,'immediately')
  options.immediately=0 ;
end ;
if ~isfield(options,'immediately_bg')
  options.immediately_bg=0 ;
end ;
if ~isfield(options,'submit_now')
  options.submit_now=1 ;
end ;
if ~isfield(options,'nicetohave')
  options.nicetohave=0 ;
end ;
if ~isfield(options, 'start_dir')
  dirctry=cwd ;
else
  dirctry=options.start_dir ;
end ;
if ~isfield(options, 'resubmit'),
  options.resubmit = 0 ;
  options.time_req_resubmit = [] ;
  options.mem_req_resubmit = [] ;
end ;
if ~isfield(options, 'data_size'),
  options.data_size = [] ;
end ;
if ~isfield(options, 'hard_time_limit')
  options.hard_time_limit=inf;
end
jobinfo=rproc_empty ;

jobinfo.ProcName = ProcName ;
jobinfo.P1 = P1 ;
jobinfo.Mem = Mem ;
jobinfo.options = options ;
jobinfo.time = time ;
jobinfo.created = 1 ;
jobinfo.resubmit = options.resubmit ;
jobinfo.mem_req_resubmit = options.mem_req_resubmit ;
jobinfo.time_req_resubmit = options.time_req_resubmit ;
jobinfo.data_size = options.data_size ;
jobinfo.hard_time_limit = options.hard_time_limit ;

% wait for the previous matlab 
% the started matlab will remove the call information file
if ~exist('~/tmp/sge')
  username = whoami ;
  base_dir = sprintf('%s/sge/tmp/%s',home_str, username) ;
  if exist(base_dir)~=7,
    [succ,tmp,tmp]=mkdir(base_dir) ; assert(succ) ;
  end ;
  tmp_dir = sprintf('%s/sge/tmp/%s/tmp',home_str, username) ;
  [succ,tmp,tmp]=mkdir(tmp_dir) ; assert(succ) ;
  sge_tmp_dir = sprintf('%s/sge/tmp/%s/tmp/sge',home_str, username) ;
  [succ,tmp,tmp]=mkdir(sge_tmp_dir) ; assert(succ) ;

  [succ,tmp,tmp]=mkdir('~/tmp') ; assert(succ) ;
  unix(sprintf('ln -s %s ~/tmp/sge',sge_tmp_dir)) ;
end ;
myassert(exist('~/tmp/sge')==7) ;

if ~exist([dirctry '/sge'])
  username = whoami ;
  sge_base_dir = strrep(dirctry, sprintf('/fml/ag-raetsch/home/%s', username), sprintf('%s/sge/tmp/%s',home_str, username)) ;
  if exist(sge_base_dir)~=7,
    succ=unix(sprintf('mkdir -p %s',sge_base_dir)) ; 
    assert(succ==0) ;
  end ;
  sge_dir = sprintf('%s/sge', sge_base_dir) ;
  [succ,tmp,tmp]=mkdir(sge_dir) ; assert(succ) ;

  unix(sprintf('ln -s %s %s/sge', sge_dir, dirctry)) ;
end ;

%option_str=sprintf(' -l h_vmem=%iM -soft -l h_cpu=%1.0f -hard ', Mem, max(60,time*60)) ;

if use_reservation,
  if Matlab_licenses>0,
    option_str=sprintf(' -R y -l h_vmem=%iM -l matlab=%1.2f -soft -l h_cpu=%1.0f -hard ', Mem, Matlab_licenses, max(60,time*60)) ;
  else
    option_str=sprintf(' -R y -l h_vmem=%iM -soft -l h_cpu=%1.0f -hard ', Mem, max(60,time*60)) ;
  end ;
else
  if Matlab_licenses>0,
    option_str=sprintf(' -l h_vmem=%iM -l matlab=%1.2f -soft -l h_cpu=%1.0f -hard ', Mem, Matlab_licenses, max(60,time*60)) ;
  else
    option_str=sprintf(' -l h_vmem=%iM -soft -l h_cpu=%1.0f -hard ', Mem, max(60,time*60)) ;
  end ;
end ;

if isequal(environment, 'galaxy'),
  option_str=[option_str ' -l parent=0.0 '] ;
end ;
  
%option_str=sprintf(' -l h_vmem=%iM -soft -l h_cpu=%1.0f -hard ', Mem, max(60,time*60)) ;

%option_str=sprintf(' -l matlab=0 -l vf=%iM', Mem) ;
%option_str=sprintf(' -l matlab=1 -l h_cpu=%i ', time*60) ;
%if Mem>1000, option_str = [option_str ' -l BIG '] ; end ;
if isfield(options,'hold')
  if options.hold==1,
	option_str=[option_str ' -h u '] ;
  end
end ;
if isfield(options,'queue') 
  option_str=[option_str sprintf(' -q "%s" ', options.queue)] ;
end ;
if isfield(options,'nicetohave') && isequal(options.nicetohave,1)
  option_str=[option_str sprintf(' -l nicetohave=1 ')] ;
end ;  

if isfield(options,'priority') 
  option_str=[option_str sprintf(' -p %i ', options.priority)] ;
end ;
if isfield(options,'express') && isequal(options.express,1)
  option_str=[option_str ' -l express '] ;
end ;

if isfield(options,'huanghonet')
  option_str=[option_str ' -l huanghonet ' num2str(options.huanghonet) 'M'] ;
end ;

if isequal(use_engine, 'matlab')
  if isequal(environment, 'internal'),
    if isfield(options,'matlab_v7_0') && isequal(options.matlab_v7_0, 1)
      bin_str='/fml/ag-raetsch/share/software/matlab-7.0/bin/matlab_rproc -nojvm -nodisplay' ;
    else
      bin_str = '/fml/ag-raetsch/share/software/matlab-7.6/bin/matlab_rproc -nojvm -nodisplay' ;
    end ;
  else
    bin_str = '/home/galaxy/matlab-7.6/bin/matlab_rproc -nojvm -nodisplay' ;
  end ;
elseif isequal(use_engine, 'octave'),
  if isequal(environment, 'internal'),
    %bin_str = '/fml/ag-raetsch/share/software/octave/x86_64/bin/octave' ;
    %bin_str = '/fml/ag-raetsch/share/software/octave-3.0.3/bin/octave';
    %bin_str = '/usr/bin/octave-3.0.0';
    bin_str = '/usr/bin/octave';
  elseif isequal(environment, 'galaxy'),
    bin_str = '/home/galaxy/octave/bin/octave' ;
  else
    error('unknown environment') ;
  end ;
 else
   error('unknown engine') ;
end ;
% request a specific system architecture
if ~isfield(options,'arch')
  option_str=[option_str ' -soft -l arch=lx26-amd64 -hard '] ;
end ;

% request cplex license
if isfield(options,'cplex') && isequal(options.cplex,1)
  option_str=[option_str ' -l cplex=1 '] ;
end ;

% request several cpus
if isfield(options,'ncpus') && ~isempty(options.ncpus) && ~equal(options.ncpus, 1)
  option_str=[option_str sprintf(' -pe "*" %i ', options.ncpus)] ;
end ;

if isfield(options, 'identifier'),
  identifier = options.identifier ;
else
  identifier = 'GR' ;
end ;

cc=round(rand(1)*100000) ;
prefix=sprintf('%s%i-%1.10f', identifier,cc(1),now) ;
mat_fname=sprintf('~/tmp/sge/%s.mat',prefix) ;
result_fname=sprintf('~/tmp/sge/%s_result.mat',prefix) ;
m_fname=sprintf('~/tmp/sge/%s.m',prefix) ;
while fexist(mat_fname) | fexist(result_fname) | fexist(m_fname),
  cc=round(rand(1)*100000) ;
  prefix=sprintf('%s%i-%1.10f', identifier,cc(1), now) ;
  mat_fname=sprintf('~/tmp/sge/%s.mat',prefix) ;
  result_fname=sprintf('~/tmp/sge/%s_result.mat',prefix) ;
  m_fname=sprintf('~/tmp/sge/%s.m',prefix) ;
end ;


clo = clock ;
log_fname = sprintf('%s/sge/%s_%s_%i_%i.rproc', dirctry, prefix, date,	clo(4), clo(5)) ;
qsublog_fname = sprintf('%s.qsubout',log_fname) ;

jobinfo.prefix = prefix ;
jobinfo.mat_fname = mat_fname ;
jobinfo.result_fname = result_fname ;
jobinfo.m_fname = m_fname ;
jobinfo.log_fname = log_fname ;
jobinfo.qsublog_fname = qsublog_fname ;

% save the call information
save(mat_fname, 'ProcName', 'P1', 'dirctry', 'options', '-v7') ;

if isequal(use_engine, 'matlab'),
  if isequal(environment, 'internal'),
    evalstring_=['addpath(''/fml/ag-raetsch/share/software/matlab_tools/rproc'') ;'] ;
    evalstring_=[evalstring_ 'addpath(''/fml/ag-raetsch/share/software/matlab_tools/utils'') ;'] ;
  elseif isequal(environment, 'galaxy'),
    evalstring_=['addpath(''/home/galaxy/svn/tools/rproc'') ;'] ;
    evalstring_=[evalstring_ 'addpath(''/home/galaxy/svn/tools/utils'') ;'] ;
  else
    error('unknown environment') ;
  end ;
else
  if isequal(environment, 'internal'),
    evalstring_=['addpath(''/fml/ag-raetsch/share/software/octave_tools/rproc'') ;'] ;
    evalstring_=[evalstring_ 'addpath(''/fml/ag-raetsch/share/software/octave_tools/utils'') ;'] ;
  elseif isequal(environment, 'galaxy'),
    evalstring_=['addpath(''/home/galaxy/svn/tools/rproc'') ;'] ;
    evalstring_=[evalstring_ 'addpath(''/home/galaxy/svn/tools/utils'') ;'] ;
  else
    error('unknown environment') ;
  end ;
end 
evalstring=['start_proc(''' mat_fname ''')'] ;
evalstring=[evalstring_ 'cd ' dirctry '; ' evalstring '; exit'] ;
fd=fopen(m_fname,'w+') ;
fprintf(fd,'%s\n',evalstring) ;
fclose(fd) ;

cl=clock ;
clstr=sprintf('%i:%i:%i--%i.%i.%i', round(cl([6 5 4 3 2 1]))) ;
mem=sprintf('%iMb',Mem) ;

if isfield(options, 'envstr'),
  envstr=options.envstr ;
  if ~isempty(envstr),
    envstr(end+1)=';' ;
  end ;
else
  envstr = '' ; 
end ;
if options.immediately,
  str=[ envstr 'cd ~/matlab; cat ' m_fname '| ' bin_str ' >>' log_fname ] ;
elseif options.immediately_bg,
  str=[ envstr 'cd ~/matlab; cat ' m_fname '| ' bin_str  ' >>' log_fname '&' ] ;
else
  % start matlab
  str=['echo "' envstr 'hostname; cd ~/matlab; cat ' m_fname '| ' ...
       bin_str ' >>' log_fname '" | qsub -o "' qsublog_fname '" -j y -r y ' ...
       option_str ' -N ' prefix ' >>& ' log_fname] ;
end ;

if options.submit_now,
  if options.verbosity>1
    disp(str)
  end ;
end ;

%unix('sync') ;

% wait until we are allowed to submit again, i.e. #jobs < maxjobs
if ~options.immediately && ~options.immediately_bg && isequal(options.waitonfull,1)
  while 1,
    [s,num_queued]=unix('qstat -u `whoami`  2> /dev/null | grep `whoami` | wc -l | tr -d " "');
    if (s~=0)
      warning('could not determine how many jobs are scheduled');
      break;
    end
    num_queued=str2double(num_queued);
    
    %keep 50 spare jobs if multiple rprocs are scheduling...
    if (num_queued < options.maxjobs)
      break;
    else
      if options.verbosity>=1
        fprintf('queue full, sleeping 60 seconds (%i/%i)\n', ...
                num_queued, options.maxjobs);
      end ;
      system('sleep 60');
    end
  end
end ;


if options.submit_now,
  if options.immediately && options.verbosity>0,
    disp('immediatedly starting job on local machine') ;
  end ;
  if options.immediately_bg && options.verbosity>0,
    disp('immediatedly starting job on local machine in background') ;
  end ;
  if options.immediately_bg,
    while(1)
      [tmp,str_]=unix('uptime');
      idx=strfind(str_, 'average:') ;
      assert(~isempty(idx)) ;
      b=separate(str_(idx+8:end), ',') ; 
      cpu_load = str2num(b{1}) ;
      if cpu_load>10,
        if options.verbosity>=1
          fprintf('load too high: %1.2f\n', cpu_load) ;
        end ;
        pause(10) ;
      else
        break ;
      end ;
    end
    pause(2)
  end ;
  
  unix(['echo ''' str '''|tcsh']) ;
  jobinfo.submission_time = now ;

  if ~options.immediately && ~options.immediately_bg,
    fd = fopen(log_fname,'r') ;
    jobinfo.jobid = -1 ;
    if fd>=0,
      s=fgetl(fd) ;
      items=separate(s,' ');
      %p = textscan(s,'Your job %d%s',1);
      if ~(isequal(items{1}, 'Your') && isequal(items{2}, 'job'))
	error('error: submission failed: %s',s);
      end
      jobinfo.jobid = str2num(items{3}) ;%p{1};
      fclose(fd) ;

      rproc_register('submit', jobinfo) ;

      %elems=separate(s,' ') ;
      %jobinfo.jobid = str2num(elems{3}) ; % third word is the jobid
      %assert(isequal(['("' jobinfo.prefix '")'], elems{4})) ;
    end ;
  else
    jobinfo.jobid = 0 ;
  end ;
else
  jobinfo.jobid = 0 ;
end ;
