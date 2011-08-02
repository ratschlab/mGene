function signal_predict(PAR_dir, genome_config_dir, fn_out, fn_out_dir, run_locally_par, save_as, options)
% signal_predict(PAR_dir, genome_config_dir, fn_out, fn_out_dir, run_locally, save_as,options)
%
% predicts signals and writes the predictions to the file fn_out
% temporary files are stored in a directory fn_out which is deleted

paths

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started signal_predict at %s\n', timedate) ;

fprintf('Results are written to: %s\n', fn_out) ;

% find internal files
[PAR_dir] = find_internal_files(PAR_dir) ;

if ~exist('run_locally_par', 'var'),
  run_locally_par = -1 ;
end ;

genome_config = [genome_config_dir '/genome.config'] ;
genome_info = init_genome(genome_config) ;

run_locally_par

run_locally=struct ;
if run_locally_par==-1,
  run_locally.pos = rproc_policy('signal_predict:pos',  genome_info) ;
  run_locally.predict = rproc_policy('signal_predict:predict',  genome_info) ;
  run_locally.conf = rproc_policy('signal_predict:conf',  genome_info) ;
  run_locally.concat = rproc_policy('signal_predict:concat',  genome_info) ;
  run_locally.save = rproc_policy('signal_predict:save',  genome_info) ;
else
  run_locally.pos = run_locally_par ;
  run_locally.predict = run_locally_par ;
  run_locally.conf = run_locally_par ;
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

 
disp('------------------------------------') ;
disp('Step 1: Setting up datastructures...') ;
disp('------------------------------------') ;

L=load([PAR_dir '/PAR.mat'], 'PAR') ;
PAR=L.PAR ;
signal = PAR.Signal_name ;
P.Signal_name = signal ;
P.RPROC = PAR.RPROC;
if exist('options','var')
  P.RPROC.options = options;
else
  %P.RPROC.options.addpaths = {};
  %P.RPROC.options.envstr = 'export LD_LIBRARY_PATH=/fml/ag-raetsch/home/jonas/shogun/trunk/src/libshogun:/fml/ag-raetsch/home/jonas/shogun/trunk/src/libshogunui';
end



% make sure there is no file with the same name
%unix(['rm -f ' fn_out]) ;

% create the working directory
%dir_orig = fn_out_dir ;
%dir_ = [fn_out_dir '/galaxy_signal_test'] ;
%unix(['mkdir -p ' dir_]) ;

fprintf(1,'\n\nstart testing for sensor %s\n\n', signal);

% prepare data structures
fid = 1; % stdout

P.SETs.num_splits = PAR.SETs.num_splits ; 
P.FN.output_sig.(signal).fn_SVMs = sprintf('%sbest_partition=', PAR.FN.output_sig.(signal).fn_SVMs);
if isfield(PAR.FN, 'base_dir')&&~fexist([P.FN.output_sig.(signal).fn_SVMs '1.mat'])
	P.FN.output_sig.(signal).fn_SVMs = [PAR.FN.base_dir '/' P.FN.output_sig.(signal).fn_SVMs];
end
if isequal(genome_config, PAR.FN.input.fn_genome_config)
  P.FN.output.fn_training_blocks = PAR.FN.output.fn_training_blocks ;
  P.FN.output.fn_test_blocks = PAR.FN.output.fn_test_blocks ;
else
  P.FN.output.fn_training_blocks = [];
  P.FN.output.fn_test_blocks = [];
end

P.FN.input.fn_genome_config = genome_config;
P.FN.input_sig.(signal).fn_pos =  sprintf('%s/pos/',fn_out_dir);
P.FN.output_sig.(signal).fn_pred =  sprintf('%s/pred/',fn_out_dir) ;

% create output directories
system(sprintf('mkdir -p %s',fn_out_dir));
system(sprintf('mkdir -p %s/pos',fn_out_dir));
system(sprintf('mkdir -p %s/pred',fn_out_dir));

fprintf('Done.\n\n') ;

disp('------------------------------') ;
disp('Step 2: Starting Prediction...') ;
disp('------------------------------') ;

% start the actual work
[jobinfo, fn_pred] = masterscript_sensors_test(P, 1, run_locally, 0, save_as) ;

fprintf('Done.\n\n') ;


disp('-------------------------------') ;
disp('Step 3: Writing output files...') ;
disp('-------------------------------') ;

% collect the output
%fn_out_temp = tempname ;
genome_info = init_genome(genome_config) ;
fd = fopen(fn_out, 'w+') ;
if fd>0
	fclose(fd);
else
	error('could not open output file: %s\n', fn_out);
end
cnt = 0 ;
spf_info = [] ;
for c = 1:length(genome_info.contig_names)
  for s = '+-'
    fprintf('processing contig_%i%s\r', c, s);

    P.chrom = c;
    P.strand = s;
  
    filename = sprintf('%scontig_%i%s', fn_pred, c, s);
    if save_as.spf_ascii,
      fn_in=[filename '_Conf_cum.spf'];
      %fn_in=[filename '_output.spf'];
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

      % cleanup
      system(['rm -f ' filename '_Conf_cum.spf ' filename '_Conf.spf ' filename '_output.spf ']) ;
      %assert(ret==0) ;
    end ;

    if save_as.spf_binary,
      fn_mat = sprintf('%s/contig_%i%s_Conf_cum_spf.mat', fn_pred, c, s);
      fnall_mat = sprintf('%s/Conf_cum_spf.mat', fn_out_dir);

      try
        SPF = load(fn_mat);
      catch
        fprintf('unable to load %s\nwaiting 120sec\n', fn_mat) ;
        pause(120)
        SPF = load(fn_mat);
      end

      spf_info(end+1).contig_name = SPF.contig_name ;
      spf_info(end).strand = SPF.strand ;
      spf_info(end).score_name = SPF.score_name ;
      spf_info(end).signal_name = SPF.signal_name ;

      save_append(fnall_mat, (cnt>0), sprintf('SPF_%i', cnt), SPF) ;
      cnt = cnt + 1 ;

      system(['rm -f ' filename '_Conf_cum_spf.mat ' filename '_Conf_spf.mat ' filename '_output_spf.mat ']) ;
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
fprintf('finished signal_predict at %s\n', timedate) ;

return
