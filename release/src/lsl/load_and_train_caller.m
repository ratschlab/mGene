function x = load_and_train_caller(P, Train_Data);
% x = load_and_train_caller(P, Train_Data);

tmpname = tempname ;
[engine,tmp,tmp,mccdir] =  determine_engine;

inventory{1}='P' ;
if nargin==1,
  save(tmpname, '-V7', 'inventory', 'P') ;
else
  inventory{2} = 'Train_Data' ;
  save(tmpname, '-V7', 'inventory', 'P', 'Train_Data') ;
end ;
  
%% use precompiled version if exist and we are not on matlab

home=getenv('HOME');
home = deblank(home);
cmp_fct = sprintf('%s/mgene_galaxy/load_and_train/run_galaxy_load_and_train_wrapper_lib.sh', home);

if 0%fexist(cmp_fct) && ~isequal(engine, 'matlab')
  fprintf('using precompiled code\n')
[ret] = system(sprintf('LD_LIBRARY_PATH=/home/galaxy/shogun.matlab_new/trunk/src/libshogun/:/home/galaxy/shogun.matlab_new/trunk/src/libshogunui %s %s %s', cmp_fct, mccdir, tmpname)) ;
  assert(ret==0) ;
 else
  fprintf('not using precompiled code\n')
  galaxy_load_and_train_wrapper(tmpname) ;
end

load(tmpname, 'x') ;
assert(exist('x','var')==1) ;
ret=system(sprintf('rm %s', tmpname)) ;
