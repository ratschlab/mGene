function mgene_predict(fasta_fname, fasta_dir, train_fname, train_dir, fn_out, output_dir, run_locally, use_old_files)
%  mGenePredict(fasta_fname, fasta_dir, train_fname, train_dir, fn_out, output_dir, run_locally, use_old_files)
%    

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started mgene_predict at %s\n', timedate) ;

fprintf('Results are written to: %s\n', fn_out) ;

if ~exist('use_old_files','var')
  use_old_files = 0;
end

% find internal files
[train_dir,train_fname] = galaxy_find_internal_files(train_dir,train_fname) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
config_dir = sprintf('%s/genome_dir',output_dir) ;
info_name = sprintf('%s/genome_info',output_dir) ;

unix(sprintf('rm -rf %s', config_dir)) ;

nice_mkdir(config_dir)
tic; galaxy_genometool(fasta_fname, fasta_dir, info_name, config_dir) ;toc;

load(sprintf('%s/mgene_predictor.mat', train_dir), 'train_dirs', 'train_files');

genome_config = [config_dir '/genome.config'] ;
genome_info = init_genome(genome_config);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIGNALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if ~exist('run_locally', 'var')
  run_locally = -1 ;
%end ;

[engine,env]=determine_engine;
if isequal(env, 'internal'),
  bdir = sprintf('/fml/ag-raetsch/home/%s/svn/projects/genefinding/',whoami);
  
  
  jobinfo = rproc_empty(0) ;
  options.start_dir = [bdir 'galaxy'];
  options.addpaths = {bdir, [bdir 'data_processing_sensors'],[bdir 'sensors'],[bdir 'contents'],[bdir 'data_processing_lsl'],[bdir '/lsl_paper'],[bdir 'galaxy'],...
                      [bdir 'utils'],[bdir 'template_files'],[bdir 'template_files/Contents'],[bdir 'template_files/Organisms'],...
                      [bdir 'template_files/Methods'],[bdir 'template_files/Signals'],[bdir 'template_files/lsl/'], [bdir 'template_files/exp']...
                      [bdir '../splicing/splicegraphs/'], [bdir '../splicing/splicegraphs/detect_altsplice/'], [bdir '../../tools/utils'],[bdir '../../tools/genomes']...
                      [bdir '../../tools/prob_plif'],...
                      '/fml/ag-raetsch/home/raetsch/tmp/src'};
  % '/fml/ag-raetsch/home/raetsch/svn/projects/genefinding/shogun'};
  %'/fml/ag-raetsch/home/jonas/shogun/trunk/src/matlab/'};
  
  options.envstr = 'export LD_LIBRARY_PATH=~/mgene_galaxy/shogun.matlab_new/trunk/src/libshogun:~/mgene_galaxy/shogun.matlab_new/trunk/src/libshogunui';
  
  % options.identifier = sprintf('%spred%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;
else
  options = [] ;
  options.force_octave = 1 ;
  options.envstr = 'export LD_LIBRARY_PATH=/home/galaxy/shogun.octave_new/trunk/src/libshogun/:/home/galaxy/shogun.octave_new/trunk/src/libshogunui';
end ;

%[train_dirs.tss, tis_dir, train_dirs.acc, train_dirs.don, train_dirs.cdsStop, train_dirs.cleave] = ...
%    sort_signals(train_dirs.tss, tis_dir, train_dirs.acc, train_dirs.don, train_dirs.cdsStop, train_dirs.cleave) ;

signal_names = {'acc', 'don', 'tis', 'cdsStop', 'tss', 'cleave'} ;

P.rproc_options = options;
P.config_dir = config_dir; 
P.train_on = 0;
P.use_old_files = use_old_files;
P.run_locally = run_locally;

num_jobs = 0;
if 1
for signal_idx = 1:length(signal_names)
  num_jobs = num_jobs + 1;
  P.signal_name = sprintf('signal_%i', signal_idx) ;%signal_names{signal_idx} ;
  
  % training files
  P.tsp_train_dir = train_dirs.(signal_names{signal_idx}) ;%P.signal_name);
  
  % prediction
  P.spf_pred_fname = sprintf('%s/signal_%s_predictions.spf',output_dir, P.signal_name) ;
  P.spf_pred_dir = sprintf('%s/signal_%s_predictions_files',output_dir, P.signal_name) ;
  
  % evaluation
  P.eval_fname = sprintf('%s/signal_%s_eval.txt',output_dir, P.signal_name) ;
 
  if run_locally==-1,
    run_locally_signal_wrapper = rproc_policy('signal_wrapper', genome_info) ;
  else
    run_locally_signal_wrapper = run_locally ;
  end ;

  if ~run_locally_signal_wrapper,
    
    [mem_req, time_req, opts] = rproc_memtime_policy(['signal_wrapper:' signal_names{signal_idx}], [], options) ;

    jobinfo(num_jobs) = rproc('signal_wrapper', P, mem_req, opts, time_req);
  else
    signal_wrapper(P);
  end
  pred_dirs.(signal_names{signal_idx}) =  P.spf_pred_dir;
  pred_files.(signal_names{signal_idx}) = P.spf_pred_fname;
end
end
clear P 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRAIN, PREDICT AND EVAL ALL CONTENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
P.config_dir = config_dir; 
P.use_old_files = use_old_files;
P.run_locally = run_locally;

content_names = {'frame0','intron', 'cds_exon', 'utr3exon', 'utr5exon', 'intergenic'} ;
P.train_on = 0;
for content_idx = 1:length(content_names) ;
  P.content_name = sprintf('content_%i', content_idx) ;%content_names{content_idx} ;
  
  % training
  try
    P.tcp_train_dir = train_dirs.(content_names{content_idx});
  catch
    warning(sprintf('content %s skipped', content_names{content_idx}));
    continue
  end
  num_jobs = num_jobs + 1;
  
  % prediction
  P.cpf_pred_fname = sprintf('%s/content_%s_predictions.spf',output_dir, P.content_name) ;
  P.cpf_pred_dir = sprintf('%s/content_%s_predictions_files',output_dir, P.content_name) ;
  
  % evaluation
  P.eval_fname = sprintf('%s/content_%s_eval.txt',output_dir, P.content_name) ;

  if run_locally==-1,
    run_locally_content_wrapper = rproc_policy('content_wrapper', genome_info) ;
  else
    run_locally_content_wrapper = run_locally ;
  end ;
  
  if ~run_locally_content_wrapper,
    P.rproc_options = options;
    [mem_req, time_req, opts] = rproc_memtime_policy(['content_wrapper:' content_names{content_idx}], [], options) ;
    
    jobinfo(num_jobs) = rproc('content_wrapper', P, mem_req, opts, time_req);
  else
    content_wrapper(P);
  end
  pred_dirs.(content_names{content_idx}) =  P.cpf_pred_dir;
  pred_files.(content_names{content_idx}) =  P.cpf_pred_fname;
end

save([output_dir '/prediction_files.mat'], 'pred_dirs', '-v7');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRAIN, PREDICT AND EVAL Gene Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(jobinfo),
  [jobinfo,num_crashed] = rproc_wait(jobinfo, 60, 1, -1);
end


logfile = [output_dir '/log_train'];

% output file of Gene_predict
% fn_predictions_gff3 = [output_dir '/gene_predictions.gff3'];


galaxy_gene_predict(	config_dir, train_files.gene_predictor, train_dirs.gene_predictor, logfile, ...
	 	pred_files.tss, pred_dirs.tss, pred_files.tis, pred_dirs.tis, pred_files.acc, pred_dirs.acc, ...
            	pred_files.don, pred_dirs.don, pred_files.cdsStop, pred_dirs.cdsStop, pred_files.cleave, pred_dirs.cleave, ...
		pred_files.intergenic, pred_dirs.intergenic, pred_files.utr5exon, pred_dirs.utr5exon, ...
		pred_files.cds_exon, pred_dirs.cds_exon, pred_files.intron, pred_dirs.intron, ...
           	pred_files.utr3exon, pred_dirs.utr3exon, fn_out, output_dir, run_locally)






[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('finished mgene_predict at %s\n', timedate) ;
