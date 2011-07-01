function mgene_train(fasta_fname, fasta_dir, gff_fname, version, seqid, source, level, fn_mGene_predictor, mGene_predictor_dir, run_locally, use_old_files, train_options, converters)
%  mgene_train(fasta_fname, fasta_dir, gff_fname, fn_mGene_predictor, mGene_predictor_dir, run_locally, version, seqid, source, level,  use_old_files, train_options, converters)

paths 
addpath(shogun_settings())

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started mgene_train at %s\n', timedate) ;

fprintf('Results are written to: %s\n', fn_mGene_predictor) ;

run_locally = -1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_dir = mGene_predictor_dir ;

config_dir = sprintf('%s/genome_dir/',output_dir)
info_name = sprintf('%s/genome_info',output_dir) ;

annotation_fname = sprintf('%s/anno_info',output_dir);
annotation_dir = sprintf('%s/',output_dir);

ignore_missing_referenced_entries = 1 ;
coding_flag = 1 ;
nice_mkdir(config_dir);
fd=fopen(fn_mGene_predictor, 'w+') ; 
if fd<1, error('could not open file %s', fn_mGene_predictor); end
fprintf(fd,'------------- \n');
fprintf(fd,'mGeneTrain started %s\n', timedate) ;
fprintf(fd,'------------- \n\n');

fprintf(fd, 'Training Based on FASTA File : %s\n',fasta_fname) ; 
fprintf(fd, '                    GFF File : %s\n',gff_fname) ; 
fclose(fd) ;

%galaxy_genometool(fasta_fname, fasta_dir, info_name, config_dir) ;
%unix(sprintf('cat %s >> %s', info_name, fn_mGene_predictor)) ;
%fd=fopen(fn_mGene_predictor, 'a+') ; fprintf(fd, '\n\n') ; fclose(fd) ;
genome_config = [config_dir '/genome.config'] ;
genome_info = init_genome(genome_config);

%if run_locally==-1,
%  run_locally_gff2anno = rproc_policy('gff2anno', genome_info);
%else
%  run_locally_gff2anno=run_locally ;
%end ;
%galaxy_gff2anno(config_dir, gff_fname, annotation_fname, annotation_dir, version, seqid, source, level, ignore_missing_referenced_entries, run_locally_gff2anno, converters) ;
%
%unix(sprintf('cat %s >> %s',annotation_fname, fn_mGene_predictor)) ;
%fd=fopen(fn_mGene_predictor, 'a+') ; fprintf(fd, '\n\n') ; fclose(fd) ;

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fd=fopen(fn_mGene_predictor, 'a+') ; 

fprintf(fd, 'mGeneTrain (started %s)\n',  timedate);
fprintf(fd, '----------------------------\n\n') ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIGNALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = [] ;
if 1%~run_locally
   
  [engine,env]=determine_engine;
  if isequal(env, 'internal'),
    bdir = get_base_dir();

	[extra_path, envstr] = shogun_settings();
    jobinfo = rproc_empty(0) ;
	options.start_dir = get_base_dir();
    options.envstr = envstr;
	options.addpaths =  {extra_path, fileparts(which('signal_wrapper')), fileparts(which('content_wrapper'))};

    % options.identifier = sprintf('%spred%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;
    MEMREQ = 10000;
    time = 60*10;
  else
    
    options.envstr = 'export LD_LIBRARY_PATH=/home/galaxy/shogun.octave_new/trunk/src/libshogun/:/home/galaxy/shogun.octave_new/trunk/src/libshogunui';

    MEMREQ = 3500;
    time = 60*10;  
  end ;

end

%signal_names = {'acc', 'don', 'tis','cdsStop', 'cleave', 'tss'} ;
signal_names = {'tis','cdsStop', 'cleave', 'tss'} ;
%signal_names = {} ;

P.covered_genic_positions_only = 1 ;
P.config_dir = config_dir ; 
P.annotation_fname = annotation_fname ;
P.annotation_dir = annotation_dir ;
P.train_on = 1 ;

% P.run_locally = run_locally;
% determine the number of contigs and 
% if there are less than ten then submit the subjobs
num_contigs = length(genome_info.contig_names); 
%clear genome_info genome_config


if exist('use_old_files')
  P.use_old_files = use_old_files;
else
  P.use_old_files = 0;
end
P.rproc_options = options;
P.run_locally = run_locally ;

num_jobs = 0;

fprintf(fd,'------------- \n')
fprintf(fd,'start signal and content training and prediction\n') ;
fprintf(fd,'--- this may take a while ---  \n') ;
fprintf(fd,'------------- \n\n')
fclose(fd) ;


for signal_idx = 1 :length(signal_names)
  num_jobs = num_jobs + 1;
  P.signal_name = signal_names{signal_idx} ;
  
  % label generation
  P.spf_label_fname = sprintf('%s/signal_%s_label.spf',output_dir, P.signal_name) ;
  P.spf_label_dir = sprintf('%s/signal_%s_label_files',output_dir, P.signal_name) ;

  % training
  P.tsp_train_info = sprintf('%s/signal_%s_classifier_info',output_dir, P.signal_name) ;
  P.tsp_train_dir = sprintf('%s/signal_%s_classifier_files',output_dir, P.signal_name) ;
  
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
	options.identifier = sprintf('lay1_%s', P.signal_name);
    [mem_req, time_req, opts] = rproc_memtime_policy(['signal_wrapper:' signal_names{signal_idx}], [], options) ;

    jobinfo(num_jobs) = rproc('signal_wrapper', P, mem_req, opts, time_req);
  else
    signal_wrapper(P);
  end
  train_dirs.(P.signal_name) =  P.tsp_train_dir;
  train_files.(P.signal_name) = P.tsp_train_info;
  pred_dirs.(P.signal_name) =  P.spf_pred_dir;
  pred_files.(P.signal_name) = P.spf_pred_fname;
end
clear P 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRAIN, PREDICT AND EVAL ALL CONTENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rmpath ~raetsch/svn/projects/genefinding/shogun
% addpath /fml/ag-raetsch/home/jonas/shogun/trunk/src/matlab

P.config_dir = config_dir; 
P.annotation_fname = annotation_fname;
P.annotation_dir = annotation_dir;
P.vs_content_names = [];
P.intergenic_dist_to_genes = [];
P.covered_genic_positions_only = 1 ;

if exist('use_old_files')
  P.use_old_files = use_old_files;
else
  P.use_old_files = 0;
end

P.run_locally = run_locally ;
% P.run_locally = run_locally;
%if num_contigs>10
%  fprintf('MORE than 10 contigs => subjobs will NOT be submitted\n');
%  P.run_locally = 1;
%else
%  fprintf('LESS than 10 contigs => subjobs will be submitted\n');
%  P.run_locally = 0;
%end
P.rproc_options = options;

content_names = {'cds_exon', 'intron', 'intergenic', 'utr5exon', 'utr3exon'} ;
%content_names = {} ;
P.train_on = 1;
for content_idx = 1:length(content_names) ;
  num_jobs = num_jobs + 1;
  P.content_name = content_names{content_idx} ;
  
  % label generation 
  P.cpf_label_fname = sprintf('%s/content_%s_label.cpf',output_dir, P.content_name) ;
  P.cpf_label_dir = sprintf('%s/content_%s_label_files',output_dir, P.content_name) ;
  
  % training
  P.tcp_train_info = sprintf('%s/content_%s_classifier_info',output_dir, P.content_name) ;
  P.tcp_train_dir = sprintf('%s/content_%s_classifier_files',output_dir, P.content_name) ;
  
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
	options.identifier = sprintf('lay1_%s', P.content_name);
    [mem_req, time_req, opts] = rproc_memtime_policy(['content_wrapper:' content_names{content_idx}], [], options) ;

    jobinfo(num_jobs) = rproc('content_wrapper', P, mem_req, opts,time_req);
  else
    content_wrapper(P);
  end
  train_dirs.(P.content_name) =  P.tcp_train_dir;
  train_files.(P.content_name) = P.tcp_train_info;
  pred_dirs.(P.content_name) =  P.cpf_pred_dir;
  pred_files.(P.content_name) = P.cpf_pred_fname;
end

if ~isempty(jobinfo),
  [jobinfo,num_crashed] = rproc_wait(jobinfo, 60, 1, -1);
end
return
save([output_dir '/prediction_files.mat'],'pred_dirs', '-v7');

fd=fopen(fn_mGene_predictor, 'a+') ; fprintf(fd, '\nUnderlying signal predictors:\n-----------------------------\n\n') ; fclose(fd) ;
for signal_idx = 1:length(signal_names)
  unix(sprintf('cat %s >> %s', train_files.(signal_names{signal_idx}), fn_mGene_predictor)) ;
  fd=fopen(fn_mGene_predictor, 'a+') ; fprintf(fd, '\n\n') ; fclose(fd) ;
end ;

fd=fopen(fn_mGene_predictor, 'a+') ; fprintf(fd, '\n\nUnderlying content predictors:\n------------------------------\n\n') ; fclose(fd) ;
for content_idx = 1:length(content_names) ;
  unix(sprintf('cat %s >> %s', train_files.(content_names{content_idx}), fn_mGene_predictor)) ;
  fd=fopen(fn_mGene_predictor, 'a+') ; fprintf(fd, '\n\n') ; fclose(fd) ;
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRAIN, PREDICT AND EVAL Gene Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


logfile = [output_dir '/log_train'];
pred_dirs.gene_predictor = sprintf('%s/gene_predictor/',output_dir);
fn_gene_predictor = sprintf('%s/gene_predictor.dat',output_dir);

save([output_dir '/trained_model.mat'],'train_dirs', '-v7');

signal_fns = ''; % spf-files if no matlab files 
content_fns = ''; % cpf-files if no matlab files 

save('workspace_mgenetrain', '-v7')
gene_train(config_dir, annotation_dir, logfile, ...
           signal_fns, pred_dirs.tss, signal_fns, pred_dirs.tis,...  
	   signal_fns, pred_dirs.acc, signal_fns, pred_dirs.don,... 
           signal_fns, pred_dirs.cdsStop, signal_fns, pred_dirs.cleave,...
           content_fns, pred_dirs.intergenic, content_fns, pred_dirs.utr5exon, ...
           content_fns, pred_dirs.cds_exon, content_fns, pred_dirs.intron,...
           content_fns, pred_dirs.utr3exon, fn_gene_predictor, pred_dirs.gene_predictor, train_options);


make_mgene_predictor(fn_gene_predictor, pred_dirs.gene_predictor, ...
			train_files.tss, train_dirs.tss, train_files.tis, train_dirs.tis,...
			train_files.acc, train_dirs.acc, train_files.don, train_dirs.don,...
			train_files.cdsStop, train_dirs.cdsStop, train_files.cleave, train_dirs.cleave, ... 
			train_files.intergenic, train_dirs.intergenic, train_files.utr5exon, train_dirs.utr5exon,...
			train_files.cds_exon, train_dirs.cds_exon, train_files.intron, train_dirs.intron,...
			train_files.utr3exon, train_dirs.utr3exon, fn_mGene_predictor, mGene_predictor_dir);


[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf(fd,'finished mgene_train at %s\n', timedate) ;


fd=fopen(fn_mGene_predictor, 'a+') ; 
fprintf(fd, '\n mGene.web Train finished : %s\n', timedate) ;  
fclose(fd) ;
