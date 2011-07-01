function gene_predict(	genome_info_dir, fn_gene_predictor, dir_gene_predictor, logfile, ...
			tss_fn, tss_dir, tis_fn, tis_dir, acc_fn, acc_dir, don_fn, don_dir, stop_fn, stop_dir, cleave_fn, cleave_dir,...
			intergenic_fn, intergenic_dir, utr5_fn, utr5_dir, exon_fn, exon_dir, intron_fn, intron_dir, utr3_fn, utr3_dir,... 
			fn_out_gff, output_dir, run_locally, train_options, reg_num)



[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started gene_predict at %s\n', timedate) ;

fprintf('Results are written to: %s\n', fn_out_gff) ;


[extra_path shogun_envstr] = shogun_settings()

% handle training options encoded in a single string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = set_default_train_opts();

opts =  read_train_opts_from_string(train_options, opts);

print_train_opts(opts, 1)

fn_train_best = find_best_model(dir_gene_predictor);
[a b c] = fileparts(fn_train_best);
fn_train_best = fullfile(pwd, a, [b c]);
unix(sprintf('ln -sf %s %s/predictor.mat', fn_train_best, dir_gene_predictor))

%which('sg')
%rmpath(fileparts(which('sg')))
addpath(extra_path)
%clear('sg')
which('sg')

%[ret, timedate] = unix('date') ;
%timedate(timedate==sprintf('\n')) = [] ;
%fprintf('started gene_predict at %s', timedate) ;

% find internal files
[dir_gene_predictor] = find_internal_files(dir_gene_predictor) ;

% sort inputs
%[tss_dir, tis_dir, acc_dir, don_dir, stop_dir, cleave_dir] = ...
%    sort_signals(tss_dir, tis_dir, acc_dir, don_dir, stop_dir, cleave_dir) ;
%[exon_dir, intron_dir, intergenic_dir, utr5_dir, utr3_dir] = ...
%    sort_contents(exon_dir, intron_dir, intergenic_dir, utr5_dir, utr3_dir) ;


fid = 1;
% set some parameters
%---------------------------
PAR.fn_log = logfile;
PAR.submit=1;
PAR.LSL.method.add_lin_feat = 0;
PAR.RPROC.MEMREQ = 4000;
PAR.RPROC.options = struct;
PAR.RPROC.exm_per_batch=1;
PAR.RPROC.options.envstr = shogun_envstr ;

PAR.cleave_offset = opts.cleave_offset;

predictor_file = [dir_gene_predictor '/predictor.mat']
load(predictor_file, 'TRAIN_PAR');
PAR.model = TRAIN_PAR.model;
PAR.fn_gene_predictor = predictor_file;

fn_config =  [genome_info_dir '/genome.config'] ;
PAR.FN.input.fn_genome_config = fn_config ;
PAR.FN.genome.fn_genome_config = fn_config ;

PAR.FN.output_sig.tss.fn_pred = [tss_dir '/pred/'] ;
PAR.FN.output_sig.tis.fn_pred = [tis_dir '/pred/'] ;
PAR.FN.output_sig.acc.fn_pred = [acc_dir '/pred/'] ;
PAR.FN.output_sig.don.fn_pred = [don_dir '/pred/'] ;
PAR.FN.output_sig.cdsStop.fn_pred = [stop_dir '/pred/'] ;
if isfield(PAR.FN.output_sig, 'polya') ;
  PAR.FN.output_sig=rmfield(PAR.FN.output_sig, 'polya') ;
end 
PAR.FN.output_sig.cleave.fn_pred = [cleave_dir '/pred/'] ;

PAR.fn_pred = [output_dir '/genome_wide_predictions/'];
nice_mkdir(output_dir);
nice_mkdir(PAR.fn_pred);

write_sig_log([PAR.fn_pred 'signal.log'], PAR);

PAR.RPROC.options.addpaths = { fileparts(which('sg'))};



% set up blocks for prediction
%------------------------------
blocks=init_regions(fn_config);

blocks = make_chunks(blocks); 

fn_pred_blocks = sprintf('%s/prediction_blocks.mat', PAR.fn_pred);
if ~fexist(fn_pred_blocks)
	save(fn_pred_blocks, 'blocks');
else
	fprintf(1,'File %s exists, will not be overwritten\n', fn_pred_blocks);
end

if exist('reg_num', 'var')
  blocks = blocks(reg_num);
end


for j=1:length(blocks), 
%  blocks(j).stop = blocks(j).stop-1;, 
  blocks(j).split = [1];
end

genome_info = init_genome(fn_config);

%%%%%%%%%%%%%%%%%%%%
% add additional tracks to blocks
%%%%%%%%%%%%%%%%%%%%
for j = 1:length(opts.track_names)
  fprintf(1,'Adding track %s from file %s\n', opts.track_names{j}, opts.track_files{j});
  blocks = feval(opts.track_functions{j}, blocks, opts.track_files{j}, opts.track_params{j});
end
for j = 1:length(opts.segment_feature_names)
  fprintf(1,'Adding segment_feature %s from file %s\n', opts.segment_feature_names{j}, opts.segment_feature_files{j});
  blocks = feval(opts.segment_feature_functions{j}, blocks, opts.segment_feature_files{j}, genome_info, opts.segment_feature_params{j});
end
%for j = 1:length(PAR.model.track_names)
%  blocks = feval(PAR.model.track_functions{j}, blocks, PAR.model.track_files{j}, PAR.model.track_params{j});
%end
%for j = 1:length(PAR.model.segment_feature_names)
%  blocks = feval(PAR.model.segment_feature_functions{j}, blocks, PAR.model.segment_feature_files{j}, genome_info, PAR.model.segment_feature_params{j});
%end


% add signals to blocks 
%------------------------------------------------------------- 
fprintf(fid, 'load signal predictions...');
blocks = make_prediction_blocks(PAR, blocks, PAR.model, 0);
fprintf(fid, 'done\n');

% add content predictions to blocks
%------------------------------------------------------------- 
fprintf(fid, 'load predictions of content sensors...');
cont_pred_dirs = {exon_dir, intron_dir, intergenic_dir, utr5_dir}; 
blocks = add_content_predictions_to_blocks(blocks, cont_pred_dirs);
fprintf(fid, 'done\n'); 

% construct DP matrix with signal features 
%-------------------------------------------------------
fprintf(fid, 'create feature matrix for dynamic programm...');
fprintf(1,'generate feature matrix for blocks \n')
tic; blocks = gen_features(blocks, PAR.model); toc
fprintf(fid, 'done\n');

%%%%%%%%%%%%%%%%%%%%
% subsample to block positions
%%%%%%%%%%%%%%%%%%%%
for i=1:length(blocks)
  blocks(i).content_pred = blocks(i).content_pred(:, blocks(i).all_pos) ;
end ;

% do predictions
%-------------------------
pred_path(PAR, blocks);
%blocks = pred_path(PAR, blocks);


% write predictions to gff file
%--------------------------------
%save('~/tmp/workspace_gene_predict', '-v7')
%genome_info = init_genome(fn_config);
%all_empty = 1;
%for j=1:length(blocks)
%  for k=1:length(blocks(j).prediction.genes)
%    if isempty(blocks(j).prediction.genes{k})
%      continue;
%    end
%    all_empty = 0;
%    if blocks(j).strand == '-'
%    %  blocks(j).prediction.genes{k}(:,1:2) = blocks(j).prediction.genes{k}(:,1:2)+1; % fix a problem of the blocks to genes function
%    end
%  end
%end
%
%if all_empty == 1
%  warning('no genes predicted');
%  return
%end
%
%clade = '';
%genes = blocks2genes(blocks, PAR.model, 'prediction', genome_info, clade);
%
%fn_genes = [fileparts(fn_out_gff) '/genes.mat'];
%fprintf(fid, 'Saving predicted genes to %s\n', fn_genes);
%save(fn_genes, 'genes')
%
%write_gff3(genes,fn_out_gff,'mGene', 'Predicted');


[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('finished gene_train at %s', timedate) ;

return

function write_sig_log(logfile, PAR)

	fd = fopen(logfile, 'w+'); 
	if fd<1
		error('could not open file %s for writing\n', logfile);
	end
	fprintf(fd, '%s\n', PAR.FN.output_sig.tss.fn_pred);
	fprintf(fd, '%s\n', PAR.FN.output_sig.tis.fn_pred);
	fprintf(fd, '%s\n', PAR.FN.output_sig.don.fn_pred);
	fprintf(fd, '%s\n', PAR.FN.output_sig.acc.fn_pred);
	fprintf(fd, '%s\n', PAR.FN.output_sig.cdsStop.fn_pred);
	fprintf(fd, '%s\n', PAR.FN.output_sig.cleave.fn_pred);
	fclose(fd)
return
