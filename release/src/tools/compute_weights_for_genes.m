function weights = compute_weights_for_genes(organism, genes) 
%% script => names of signal prediction directories
rgasp_file_names

[tmp fn_config] = rgasp_genome_config_dir(organism);
[extra_path shogun_envstr] = shogun_settings();

addpath(extra_path)
%clear('sg')
%which('sg')


fd_out = 1;

% handle training options encoded in a single string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%configuration.bam_exon_track = 1;
%configuration.bam_intron_track = 1;
%configuration.bam_intron_list = 1; 

train_options = '';
%train_options = define_tracks(configuration, '', organism);

opts = set_default_train_opts();
opts =  read_train_opts_from_string(train_options, opts);
%print_train_opts(opts, fd_out)


galaxy_PAR.dir = '/fml/ag-raetsch/home/jonas/tmp/features2/';
galaxy_PAR.genome_config = fn_config;

galaxy_PAR.track_names = opts.track_names;
galaxy_PAR.track_functions = opts.track_functions;
galaxy_PAR.track_files = opts.track_files;
galaxy_PAR.track_params = opts.track_params;
galaxy_PAR.track_monoton_functions = opts.track_monoton_functions;

galaxy_PAR.segment_feature_names = opts.segment_feature_names;
galaxy_PAR.segment_feature_functions = opts.segment_feature_functions;
galaxy_PAR.segment_feature_files = opts.segment_feature_files;
galaxy_PAR.segment_feature_params = opts.segment_feature_params;
galaxy_PAR.segment_feature_monoton_functions = opts.segment_feature_monoton_functions;

PAR = create_PAR('exp_galaxy', galaxy_PAR);
PAR.tasks.run_locally=1 ;

create_directories(PAR.FN, PAR.organism) ;

PAR.FN.input.fn_genome_confiG = fn_config ;
PAR.FN.genome.fn_genome_config = fn_config ;
PAR.model.long_trans_thresh = opts.long_trans_thresh;

C_lin = 0;
C_sq = opts.C;

PAR.FN.output_sig.tss.fn_pred = [tss_dir '/pred/'] ;
PAR.FN.output_sig.tis.fn_pred = [tis_dir '/pred/'] ;
PAR.FN.output_sig.acc.fn_pred = [acc_dir '/pred/'] ;
PAR.FN.output_sig.don.fn_pred = [don_dir '/pred/'] ;
PAR.FN.output_sig.cdsStop.fn_pred = [stop_dir '/pred/'] ;
PAR.FN.output_sig.cleave.fn_pred = [cleave_dir '/pred/'] ;
if isfield(PAR.model.signals, 'polya')
  PAR.model.signals = rmfield(  PAR.model.signals, 'polya') ;
end ;

genome_info = init_genome(fn_config);
genes = get_chr_num(genes, genome_info);

for i=1:length(genes),
  for j=1:length(genes(i).transcripts)
    genes(i).start = min(genes(i).start, min(min(genes(i).exons{j}(:,1:2)))) ;
    genes(i).stop  = max(genes(i).stop, max(max(genes(i).exons{j}(:,1:2)))) ;
  end ;
  genes(i).id = i;
end ;

for i=1:length(genes),
	genes(i).transcript_status = ones(1,length(genes(i).transcripts)) ;
	genes(i).transcript_valid = ones(1,length(genes(i).transcripts)) ;
	genes(i).transcript_complete = ones(1,length(genes(i).transcripts)) ;
	genes(i).is_complete = 1;
	if isfield(genes, 'utr_5prime')
		genes(i).utr5_exons = genes(i).utr_5prime;
		genes(i).utr3_exons = genes(i).utr_3prime;
	end
end

maxlen = 1000;
min_len = 500;
rm_dist = 500;
 
  [blocks rm_idx] = create_blocks_from_genes(genes, maxlen, min_len, rm_dist, fn_config);
  blocks(unique([rm_idx])) = [];

[blocks, genes, split] = gen_lsl_blocks(genes, blocks, PAR.FN.genome.fn_genome_config, PAR.model, '');

for j = 1:length(opts.track_names)
  blocks = feval(opts.track_functions{j}, blocks, opts.track_files{j}, opts.track_params{j});
end
for j = 1:length(opts.segment_feature_names)
  blocks = feval(opts.segment_feature_functions{j}, blocks, opts.segment_feature_files{j}, genome_info, opts.segment_feature_params{j});
end

%%%%%%%%%%%%%%%%%%%%
% add signals to blocks 
%%%%%%%%%%%%%%%%%%%%
blocks = make_prediction_blocks(PAR, blocks, PAR.model, 0);




%%%%%%%%%%%%%%%%%%%%
% add content predictions to blocks
%%%%%%%%%%%%%%%%%%%%
cont_pred_dirs = {exon_dir, intron_dir, intergenic_dir, utr5_dir, utr3_dir}; 
blocks = add_content_predictions_to_blocks(blocks, cont_pred_dirs);

blocks = gen_features(blocks, PAR.model);

%%%%%%%%%%%%%%%%%%%%
% subsample to block positions
%%%%%%%%%%%%%%%%%%%%
for i=1:length(blocks)
  blocks(i).content_pred = blocks(i).content_pred(:, blocks(i).all_pos) ;
end ;

%%%%%%%%%%%%%%%%%%%%
% move tss and cleave signal positions in true segmentation to nearest 
% position with tss or cleave prediction
%%%%%%%%%%%%%%%%%%%%
rm_idx = [];
for b = 1:length(blocks)
  for t=1:length(blocks(b).truth)
   	new_truth = correct_cleave_and_tss_pos(blocks(b).truth(t), blocks(b), PAR.model);
    if isempty(new_truth), 
      rm_idx = [rm_idx b];
      continue ;
    end
    blocks(b).truth(t) = new_truth ;
  end
end
if length(rm_idx)>0
	Tfprintf([1 fd_out], 'Remove %i of %i blocks because no tss or cleave position was found\n', length(rm_idx), length(blocks));
	blocks(rm_idx) = [];
end
if length(rm_idx)==length(genes)
	weights.no_signal_cons = 1;
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%
% create paths from segmentation
%%%%%%%%%%%%%%%%%%%%%%%%
take_map = logical(ones(1,length(blocks))) ;

for b = 1:length(blocks)

  if mod(b, 10)==0
  	fprintf('\r %i (%i) blocks done', b, length(blocks))
  end
  for t=1:length(blocks(b).truth)
    [blocks(b).truth(t).path,blocks(b).truth(t).pos_idx,blocks(b).truth(t).pos]=...
         segmentation2path(blocks(b).truth(t),blocks(b),PAR.model);
    % check that the inverse function produces again the original segmentation
    [segments,tmp,ok] = path2segmentation(blocks(b).truth(t).path, blocks(b).truth(t).pos_idx, blocks(b), PAR.model);
    if ok,
      assert(isequal(segments,blocks(b).truth(t).segments(:,1:3)))
    else
      take_map(b)=false ;
    end ;
    clear segments
  end
end
blocks = blocks(take_map) ;
%%%%%%%%%%%%%%%%%%%%%%%%
% determine length distributions for segment types
%%%%%%%%%%%%%%%%%%%%%%%%
%PAR.model = prepare_model(PAR,blocks) ;
%
if ~isempty(strmatch(organism, {'Col_0', 'Bur_0', 'Can_0', 'Ct_1', 'Edi_0', 'Hi_0', 'Kn_0', 'Ler_0', 'Mt_0', 'No_0', 'Oy_0', 'Po_0', 'Rsch_4', 'Sf_2', 'Tsu_0', 'Wil_2', 'Ws_0', 'Wu_0', 'Zu_0'}, 'exact'))
	load('/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/Col_0/lsl/RNA_seq_ns/output/lsl/data/boundary_model.mat', 'boundaries')
	PAR.model.boundaries = boundaries;
else
	error('specify boundary file')
end 

%%%%%%%%%%%%%%%%%%%%%%%%
% compute weights for the true paths
%%%%%%%%%%%%%%%%%%%%%%%%


blocks = initialize_blocks(blocks, PAR.model, 1, PAR.FN);
%Tfprintf([1 fd_out], 'done\n');


weights = blocks(1).truth.weights;
%for j = 2:length(blocks)
%	names = fieldnames(weights);
%	for k = 1:length(names)
%		weights.(names{k}) = weights.(names{k})+blocks(j).truth.weights.(names{k});
%	end
%end

%save(fn_average_weight, 'weights')
%keyboard
%save_struct(fn_average_weight, blocks, 'blocks')
