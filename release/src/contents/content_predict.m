function content_predict(PAR_dir, genome_config_dir, fn_out, fn_out_dir, run_locally_par, save_as,options)

paths

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started content_predict at %s\n', timedate) ;

fprintf('results are written to: %s\n', fn_out) ;

% find internal files
[PAR_dir] = find_internal_files(PAR_dir) ;

[extra_path shogun_envstr] = shogun_settings()
addpath(extra_path)

genome_config = [genome_config_dir '/genome.config'] ;
genome_info = init_genome(genome_config) ;

%run_locally_par==0 ;

if ~exist('run_locally_par', 'var'),
  run_locally_par = -1 ;
end ;

if run_locally_par==-1,
  run_locally.pos = rproc_policy('content_predict:pos',  genome_info) ; % not implemented
  run_locally.predict = rproc_policy('content_predict:predict',  genome_info) ;
  run_locally.concat = rproc_policy('content_predict:concat',  genome_info) ;
  run_locally.save = rproc_policy('content_save:save',  genome_info) ;
else
  run_locally.pos = run_locally_par ; % not implemented 
  run_locally.predict = run_locally_par ;
  run_locally.concat = run_locally_par ;
  run_locally.save = run_locally_par ;
end ;

engine = determine_engine() ;
if isequal(engine, 'octave'),
  struct_levels_to_print(1) ;
end ;
run_locally
if isequal(engine, 'octave'),
  struct_levels_to_print(0) ;
end ;


if ~exist('save_as', 'var'),
  save_as.spf_ascii = 1 ;
  save_as.spf_binary = 1 ;
  save_as.wiggle = 0 ;
end ;
if ~exist('options', 'var'),
  options = [] ;
end ;


disp('------------------------------------') ;
disp('Step 1: Setting up datastructures...') ;
disp('------------------------------------') ;

L=load([PAR_dir '/PAR.mat'], 'PAR') ;
PAR=L.PAR ;
content = PAR.Content_name ;
P.Content_name = content ;
P.RPROC = PAR.RPROC ;
P.RPROC.options = options;
if ~isfield(P.RPROC.options, 'addpaths')
  P.RPROC.options.addpaths={} ;
end ;
P.RPROC.options.addpaths{end+1}=extra_path ;
P.RPROC.options.addpaths{end+1}=fileparts(which('load_and_predict')) ;
P.RPROC.options.envstr = shogun_envstr ;

genome_config = [genome_config_dir '/genome.config'] ;

% make sure there is no file with the same name
%system(['rm -f ' fn_out]) ;

% create the working directory
%dir_orig = fn_out_dir ;
%system(['mkdir -p ' dir_]) ;

fprintf(1,'\n\nstart testing for content %s\n\n', content);

% prepare data structures
fid = 1; % stdout

P.SETs.num_splits = PAR.SETs.num_splits ; 
P.FN.output_cont.(content).fn_SVMs = sprintf('%sbest_partition=', PAR.FN.output_cont.(content).fn_SVMs);
if isequal(genome_config, PAR.FN.input.fn_genome_config)
  P.FN.output.fn_training_blocks = PAR.FN.output.fn_training_blocks ;
  P.FN.output.fn_test_blocks = PAR.FN.output.fn_test_blocks ;
else
  P.FN.output.fn_training_blocks = [];
  P.FN.output.fn_test_blocks = [];
end

P.FN.input.fn_genome_config = genome_config;
P.FN.input_cont.(content).fn_pos =  sprintf('%s/pos/',fn_out_dir);
P.FN.output_cont.(content).fn_pred =  sprintf('%s/pred/',fn_out_dir) ;

% create output directories
system(sprintf('mkdir -p %s',fn_out_dir));
system(sprintf('mkdir -p %s/pos',fn_out_dir));
system(sprintf('mkdir -p %s/pred',fn_out_dir));

fprintf('Done.\n\n') ;

disp('------------------------------') ;
disp('Step 2: Starting Prediction...') ;
disp('------------------------------') ;

% start the actual work

[jobinfo,fn_pred] = masterscript_contents_test(P, 0, run_locally, 0, save_as) ;

fprintf('Done.\n\n') ;

disp('-------------------------------') ;
disp('Step 3: Writing output files...') ;
disp('-------------------------------') ;

% collect the output
genome_info = init_genome(genome_config) ;
fclose(fopen(fn_out, 'w+')) ;
cnt = 0 ;
spf_info = [] ;
for c = 1:length(genome_info.contig_names)
  for s = '+-'
    fprintf('processing contig_%i%s\r', c, s);

    P.chrom = c;
    P.strand = s;
    
    filename = sprintf('%scontig_%i%s', fn_pred, c, s);
    if save_as.spf_ascii,
      fn_in=[filename '_output.spf'];
      for j=1:20
        try
          tmp = fopen(fn_in, 'r');
          assert(tmp>-1) ;
          fclose(tmp);
        catch
  	  fprintf('unable to load file %s\nwaiting 30 secs\n', fn_in)
          pause(30)
        end
      end
      system(['cat ' fn_in ' >> ' fn_out]) ;
      %assert(ret==0) ;
      system(['rm -f ' filename '_output.spf ']) ;
      %assert(ret==0) ;
    end ;

    if save_as.spf_binary,
      fn_mat = sprintf('%s/contig_%i%s_output_spf.mat', fn_pred, c, s) ;
      fnall_mat = sprintf('%s/output_spf.mat', fn_out_dir) ;

      try %if file exists, but writing is not completed
        SPF=load(fn_mat);
      catch
        fprintf('unable to load %s\nwaiting 300sec\n', fn_mat) ;
	pause(10)
        SPF=load(fn_mat);
      end
      spf_info(end+1).contig_name = SPF.contig_name ;
      spf_info(end).strand = SPF.strand ;
      spf_info(end).score_name = SPF.score_name ;
      spf_info(end).signal_name = SPF.signal_name ;

      save_append(fnall_mat, (cnt>0), sprintf('SPF_%i', cnt), SPF) ;
      cnt = cnt + 1 ;

      % cleanup
      system(['rm -f ' filename '_output_spf.mat ']) ;
      %assert(ret==0) ;
    end ;
    system(['rm -f ' filename '_split=?.mat ' filename '_all.mat']) ;
    %assert(ret==0) ;
  end ;
end ;

if save_as.spf_binary,
  save_append(fnall_mat, (cnt>0), 'SPF_info', spf_info) ;
end ;

fprintf('Done.\n\n') ;


[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('finished content_predict at %s\n', timedate) ;
