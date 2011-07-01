function [fn_predictor, iter] = train_path_caller(PAR)
% [fn_predictor, iter] = train_path_caller(PAR)

tmpname = [tempname '.mat'];
[engine,tmp,tmp,mccdir] =  determine_engine;

PAR.save_fname = tmpname ;

inventory{1}='PAR' ;
save(tmpname, '-V7', 'inventory', 'PAR') ;
clear PAR
  
%% use precompiled version if exist and we are not on matlab

home=getenv('HOME');
home = deblank(home);
cmp_fct = sprintf('%s/mgene_galaxy/train_path/run_galaxy_train_path_wrapper_lib.sh', home);

if 0%fexist(cmp_fct) && ~isequal(engine, 'matlab')
  fprintf('using precompiled code\n')
  [ret] = system(sprintf('LD_LIBRARY_PATH=/home/galaxy/shogun.matlab_new/trunk/src/libshogun/:/home/galaxy/shogun.matlab_new/trunk/src/libshogunui %s %s %s', cmp_fct, mccdir, tmpname)) ;
  %assert(ret==0) ;
 else
  fprintf('not using precompiled code\n')
  train_path_wrapper(tmpname) ;
end

load(tmpname, 'fn_predictor', 'iter') ;
%assert(exist('fn_predictor','var')) ;
%assert(exist('iter','var')) ;

ret=system(sprintf('rm %s', tmpname)) ;
