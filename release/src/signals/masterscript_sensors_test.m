function [jobinfo,fn_pred]=masterscript_sensors_test(P, check_files, run_locally, do_gzip, save_as, keep_all_files)
% [jobinfo ,all_done]=masterscript_sensors_test(P, check_files, run_locally, do_gzip, save_as, keep_all_files)

num_splits = P.SETs.num_splits ; 
Signal_name = P.Signal_name;

FN = P.FN;  
fn_genome_config = FN.input.fn_genome_config;
fn_training_blocks = FN.output.fn_training_blocks;  
fn_test_blocks = FN.output.fn_test_blocks;

fn_pos =  FN.input_sig.(Signal_name).fn_pos;
fn_pred = FN.output_sig.(Signal_name).fn_pred;
fn_svms = FN.output_sig.(Signal_name).fn_SVMs ;
% fn_svms = sprintf('%sbest_partition=', FN.output_sig.(Signal.name).fn_SVMs) ;

L = load([fn_svms '1.mat'], 'Train');
Signal =  L.Train.PAR.Signal;
clear L 

if nargin<2 
  check_files = 1;
end

if nargin<3 
  run_locally.pos = 1 ;
  run_locally.predict = 1 ;
  run_locally.conf = 1 ;
  run_locally.concat = 1 ;
  run_locally.save = 1 ;
end

if nargin<4 
  do_gzip = 0;
end

if nargin<5 
  save_as.spf_ascii = 1 ;
  save_as.spf_binary = 1 ;
  save_as.wiggle = 0 ;
end

if nargin<6 
  keep_all_files = 0 ;
end

%if ~run_locally  
RPROC = P.RPROC;
%else
%  RPROC = struct;
%end

fid =1;
fprintf(fid,'\n\n-----------\n');
fprintf(fid,'SENSOR :  %s\n',Signal.name);
fprintf(fid,'-----------\n');

all_done = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Generate Positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jobinfo=struct ;

%[temp, num_files] = unix(sprintf('ls %scontig_*.pos | wc -l',fn_pos));
d=dir(sprintf('%scontig_*.pos',fn_pos));
num_files = length(d) ;
%[temp, num_files1] = unix(sprintf('ls %scontig_*.svm | wc -l',fn_pos));
d=dir(sprintf('%scontig_*.svm',fn_pos));
num_files1=length(d) ;

genome_info = init_genome(fn_genome_config);
if isempty(num_files)||(num_files<length(genome_info.contig_names)*2)||...
      isempty(num_files1)||(num_files1<length(genome_info.contig_names)*2)
  jobinfo = gen_all_positions(fn_genome_config,fn_pos,Signal,fn_training_blocks,fn_test_blocks,num_splits,run_locally.pos,RPROC);
  if ~isempty(jobinfo),
    fprintf('generating genome position files\n')
    jobinfo = rproc_resubmit(jobinfo);
    jobinfo = rproc_wait(jobinfo, 20, 1, -1); 
  end ; 
else
  fprintf('genome position files already exist\n')
end

avg_all_svms = 0;
subtract_diff = 0;
jobinfo = genomewide_predictions(fn_genome_config,fn_pos,fn_pred,fn_svms,num_splits,avg_all_svms,run_locally.predict, RPROC);
fprintf('all %s predictions done\n',Signal.name)

jobinfo = convert_predictions2confs(fn_genome_config,fn_pred,fn_svms,num_splits, run_locally.conf, RPROC);
fprintf('conversion done\n')


jobinfo = concat_genomewide_predictons(fn_genome_config,fn_pos,fn_pred, num_splits,avg_all_svms,subtract_diff,run_locally.concat, RPROC);
fprintf('concatination done\n')


resolution = Signal.export_settings.resolution ;


conf_cum_thresh = Signal.export_settings.conf_cum_thresh;

inventory = {'output', 'pos', 'Conf', 'Conf_cum'} ;
[jobinfo] = save_predictions(fn_genome_config, fn_pred, ...
                             Signal_name, resolution, ...
                             conf_cum_thresh, '+-', ...
                             run_locally.save, RPROC, ...
                             do_gzip, save_as, inventory);
fprintf('saving done\n')


if check_files
  % fprintf('checking files \n')
  % check_genomewide_output(P) % this fct compares prediction results with test results from training (cross_validation)
  % all_good = check_predictions(fn_genome_config,fn_pred,Signal_name,resolution)
else
  fprintf('WARNING: files not checked\n')
end


fprintf('%s COMPLETE !\n\n',Signal.name)


if ~keep_all_files
 %  unix(sprintf('rm -rf %s*_split=*', FN.output_sig.(Signal_name).fn_pred));
end
