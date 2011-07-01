function gff2anno_helper(Data)
% gff2anno_helper(Data)
paths

[engine,tmp,tmp,mccdir] =  determine_engine;

gff_fname  = Data.gff_fname;
save_fname = Data.save_fname;
source     = Data.source; 
cname      = Data.cname;
TYPE       = Data.TYPE;
log_fname  = Data.log_fname;
genome_config = Data.genome_config;
do_correct_tis_stop = Data.do_correct_tis_stop ;

level3_coding_names= Data.level3_coding_names;
  
%% use precompiled version if exist

home=getenv('HOME');
home = deblank(home);
cmp_fct = sprintf('%s/mgene_galaxy/gff2anno/run_galaxy_gff2anno_helper2.sh',home);
par_name = [save_fname 'par.mat'];
save(par_name, '-v7', 'TYPE', 'cname', 'source','level3_coding_names', 'do_correct_tis_stop');

if 0%fexist(cmp_fct) && ~isequal(engine, 'matlab')
  fprintf('using precompiled code\n')
  [ret] = system(sprintf('%s %s %s %s %s %s %s %s %s %s',...
			      cmp_fct, mccdir, gff_fname, save_fname, par_name, log_fname, genome_config)) ;
  %[ret,output] = unix(sprintf('%s %s %s %s %s %s %s %s %s %s',...
%			      cmp_fct, mccdir, gff_fname, save_fname, par_name, log_fname, genome_config), '-echo') ;
  if ~isequal(ret, 0),
    ret
  end ;
 else
  fprintf('not using precompiled code\n')
  gff2anno_helper2(gff_fname, save_fname, par_name, log_fname, genome_config) ;
end
