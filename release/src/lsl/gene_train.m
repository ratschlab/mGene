function gene_train(	genome_info_dir, anno_dir, logfile, tss_fn, tss_dir, tis_fn, tis_dir, ...
			acc_fn, acc_dir, don_fn, don_dir, stop_fn, stop_dir, cleave_fn, cleave_dir,...
			intergenic_fn, intergenic_dir, utr5_fn, utr5_dir, exon_fn, exon_dir, intron_fn, intron_dir, ...
			utr3_fn, utr3_dir, fn_out, output_dir, train_options, chr_num)
% gene_train(	genome_info_dir, anno_dir, logfile, tss_fn, tss_dir, tis_fn, tis_dir, ...
%		acc_fn, acc_dir, don_fn, don_dir, stop_fn, stop_dir, cleave_fn, cleave_dir,...
%		intergenic_fn, intergenic_dir, utr5_fn, utr5_dir, exon_fn, exon_dir, intron_fn, intron_dir, ...
%               utr3_fn, utr3_dir, fn_out, output_dir, maxNumIter)


[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started gene_train at %s\n', timedate) ;

fprintf('Logfile is: %s\n', fn_out) ;
nice_mkdir(fileparts(fn_out))
[fd_out msg] = fopen(fn_out, 'w+') ;
if fd_out<1, error('could not open file %s: %s',fn_out, msg); end
Tfprintf([1 fd_out], 'Trained gene predictor (training started %s)\n\n', timedate) ;

disp('------------------------------------') ;
disp('START GENE TRAIN');
disp('------------------------------------') ;


[extra_path shogun_envstr] = shogun_settings()

ld_lib_path = getenv('LD_LIBRARY_PATH')

%which('sg')
%rmpath(fileparts(which('sg')))
addpath(extra_path)
clear('sg')
which('sg')

if ~exist('train_options', 'var'),
  train_options = '' ;
end ;


% handle training options encoded in a single string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = set_default_train_opts();

opts =  read_train_opts_from_string(train_options, opts);

print_train_opts(opts, fd_out)


disp('------------------------------------') ;
disp('Obtain parameter structure');
disp('------------------------------------') ;


fn_config =  [genome_info_dir '/genome.config'] ;

galaxy_PAR.dir = output_dir ;
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
if isfield(opts, 'num_expression_levels')&&opts.num_expression_levels>1
	PAR.model = extent_model_to_expression_levels(PAR.model, opts.num_expression_levels)
end
PAR.tasks.run_locally=1 ;

create_directories(PAR.FN, PAR.organism) ;

PAR.FN.input.fn_genome_config = fn_config ;
PAR.FN.genome.fn_genome_config = fn_config ;
organism = PAR.organism.name;
PAR.model.long_trans_thresh = opts.long_trans_thresh;

%struct_levels_to_print(2)
%PAR.LSL.method.par_ms.C_regul
C_lin = 0;
C_sq = opts.C;

PAR.LSL.method.max_num_iterations = opts.maxNumIter;
PAR.LSL.method.par_ms.C_regul.contents = C_lin ;
PAR.LSL.method.par_ms.C_regul.signals = C_lin ;
PAR.LSL.method.par_ms.C_regul.lengths = C_lin ;
PAR.LSL.method.par_ms.C_regul.tracks = C_lin ;
PAR.LSL.method.par_ms.C_regul.segment_features = C_lin ;
PAR.LSL.method.par_ms.C_regul.rna_seq_intron_list = C_lin ;
PAR.LSL.method.par_ms.C_regul.rna_seq_intron_quality = C_lin ;
PAR.LSL.method.par_ms.C_regul.transitions = C_lin ;
PAR.LSL.method.par_ms.C_regul.transitions_sq = C_sq ;
PAR.LSL.method.par_ms.C_regul.plif_ys_sq = C_sq ;
PAR.LSL.method.par_ms.C_regul.smoothness_sq = C_sq*10;

%[tss_dir, tis_dir, acc_dir, don_dir, stop_dir, cleave_dir] = ...
%    sort_signals(tss_dir, tis_dir, acc_dir, don_dir, stop_dir, cleave_dir) ;
%
%ok = check_signal(tss_dir, 'tss') ; assert(ok) ;
%ok = check_signal(tis_dir, 'tis') ; assert(ok) ;
%ok = check_signal(acc_dir, 'acc') ; assert(ok) ;
%ok = check_signal(don_dir, 'don') ; assert(ok) ;
%ok = check_signal(stop_dir, 'cdsStop') ; assert(ok) ;
%ok = check_signal(cleave_dir, 'cleave') ; assert(ok) ;

PAR.FN.output_sig.tss.fn_pred = [tss_dir '/pred/'] ;
PAR.FN.output_sig.tis.fn_pred = [tis_dir '/pred/'] ;
PAR.FN.output_sig.acc.fn_pred = [acc_dir '/pred/'] ;
PAR.FN.output_sig.don.fn_pred = [don_dir '/pred/'] ;
PAR.FN.output_sig.cdsStop.fn_pred = [stop_dir '/pred/'] ;
if isfield(PAR.FN.output_sig, 'polya') ;
  PAR.FN.output_sig=rmfield(PAR.FN.output_sig, 'polya') ;
end 
PAR.FN.output_sig.cleave.fn_pred = [cleave_dir '/pred/'] ;

if isfield(PAR.model.signals, 'polya')
  PAR.model.signals = rmfield(  PAR.model.signals, 'polya') ;
end ;

disp('------------------------------------') ;
disp('Loading genes and fixing boundaries/splicegraphs');
disp('------------------------------------') ;

Tfprintf([1 fd_out], 'Start preprocessing of genes for training:\n');
genes = load_genes(anno_dir);

if 0%~isfield(genes, 'transcript_weights')
	Tfprintf([1 fd_out], 'Getting quantification values for transcripts\n');
	genes = get_quantification_values(genes, PAR.FN.genome.fn_genome_config)
end

if exist('chr_num', 'var')
	genes = genes([genes.chr_num]==chr_num);
end

genome_info = init_genome(fn_config);
genes = get_chr_num(genes, genome_info);
gg = filter_invalid_genes(genes, genome_info);

%% I consider genes suspicious where the minimal transcript start or 
%% maximal transcript stop are not equal to the gene start or gene 
%% stop, respectively. It is likely that there were more transcripts 
%% that were filtered out by some previous steps or couldn't 
%% be parsed from the gff
%rm_idx = find_suspicious_genes(gg);
%Tfprintf([1 fd_out], 'remove %i (%i) genes because their start/end do not match the transcript boundaries\n', length(rm_idx), length(gg));
%gg(rm_idx) = [];

%% 
rm_idx = find_genes_with_nonoverlaping_transcripts(gg);
Tfprintf([1 fd_out], 'remove %i (%i) genes because they have nonoverlapping transcripts\n', length(rm_idx), length(gg));
gg(rm_idx) = [];

%rm_idx = find([gg.chr_num]>10);
%Tfprintf([1 fd_out], 'remove %i (%i) genes with chr_num >10\n', length(rm_idx), length(gg));
%gg(rm_idx) = [];

%max_num_genes = 10000;
%if length(gg)>max_num_genes
%	perm = randperm(length(gg));
%	rm_idx = perm(max_num_genes+1:end);
%	Tfprintf([1 fd_out], 'keep %i (%i) \n', max_num_genes, length(gg));
%	gg(rm_idx) = [];
%end


if opts.use_train_region
	l = load(opts.train_region_file);
	train_genes = [];
	for j = 1:length(l.blocks)
		region_genes_idx = find([gg.chr_num]==l.blocks(j).chr_num&[gg.start]>l.blocks(j).start&[gg.stop]<l.blocks(j).stop);
		%region_genes_idx = find([gg.chr_num]==l.blocks(j).chr_num&[gg.strand]==l.blocks(j).strand&[gg.start]>l.blocks(j).start&[gg.stop]<l.blocks(j).stop);
		train_genes = [train_genes, gg(region_genes_idx)];
	end
	Tfprintf([1 fd_out], 'remove %i (%i) genes because they live outside of training regions\n', length(gg)-length(train_genes), length(gg));
	gg = train_genes;
	clear train_genes region_genes_idx l 
end


keep_idx = find(ismember([genes.id], [gg.id]));
remove_idx = find(~ismember([genes.id], [gg.id]));
clear gg
if isempty(keep_idx)
	return
end

%genes = genes(~[genes.is_alt]) ;
%genes = genes(1:3) ;
for i=1:length(genes),
  for j=1:length(genes(i).transcripts)
    genes(i).start = min(genes(i).start, min(min(genes(i).exons{j}(:,1:2)))) ;
    genes(i).stop  = max(genes(i).stop, max(max(genes(i).exons{j}(:,1:2)))) ;
  end ;
  genes(i).id = i;
end ;
%% caused error "Subscripted assignment between dissimilar structures"
%% because field splicegraph was removed and then added again in a different place
%genes(keep_idx) = update_splicegraph(genes(keep_idx), PAR.tasks.run_locally) ;
Tfprintf([1 fd_out], 'Done.\n');

Tfprintf([1 fd_out], 'Generated %i usable gene structures from genome annotation\n\n', length(keep_idx)) ;

disp('------------------------------------') ;
disp('generate blocks from genes');
disp('------------------------------------') ;

if ~isfield(genes, 'transcript_status') || length([genes.transcript_status])==0
	for i=1:length(genes),
		genes(i).transcript_status = ones(1,length(genes(i).transcripts)) ;
	end
else
	% append transcript status 0 where nesseccary
	for i=1:length(genes)
		if length(genes(i).transcript_status)<length(genes(i).transcripts)
			genes(i).transcript_status(end+1:length(genes(i).transcripts))=0;
		end
	end
end

if isequal(opts.block_design, 'regionauto')
  %% create blocks from a genes structure that contains consecutive streches of 
  %% genes as comprehensive as possible 
  %% worst case is if only every second gene is contained in the gene structure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %maxlen = 4000;
  %min_len = 1000;
  %rm_dist = 100;
  if ~isempty(strfind(fn_config, 'human'))
    maxlen = 30000;
    min_len = 3000;
    rm_dist = 900;
  else
	% elegans and drosophila
    maxlen = 15000;% elegans: about 1.4% are longer  
    min_len = 300;% elegans: about 1.1% are shorter
    rm_dist = 50;
  end

  Tfprintf([1 fd_out], 'create blocks assuming that regions between genes \n')
  Tfprintf([1 fd_out], 'are true intergenic regions:\n')
  Tfprintf([1 fd_out], '\t\t use intergenic regions up to a length of %i nt\n', maxlen)
  Tfprintf([1 fd_out], '\t\t if available use at least  %i nt of intergenic region\n', min_len)

  [blocks rm_idx] = create_blocks_from_genes(genes, maxlen, min_len, rm_dist, fn_config);
  Tfprintf([1 fd_out], 'remove %i (%i) blocks because their intergenic region is shorter than %i\n', length(rm_idx), length(blocks), rm_dist);
  Tfprintf([1 fd_out], 'remove %i (%i) blocks because the corresponding genes were filtered out\n', length(remove_idx), length(blocks)); 
  blocks(unique([remove_idx rm_idx])) = [];

elseif isequal(opts.block_design, 'regionlist')
  %% create blocks from a list or regions stored in a ascii file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  Tfprintf([1 fd_out], 'parse training blocks from file %s \n', opts.block_design_regionlist)
  blocks = parse_region_list(opts.block_design_regionlist);
  [blocks rm_idx] = assign_genes_to_blocks(blocks, genes);
  Tfprintf([1 fd_out], 'remove %i (%i) blocks because they cut annotated genes\n', length(rm_idx), length(blocks));
  blocks(rm_idx) = [];
  rm_idx = find_empty_blocks(blocks);
  Tfprintf([1 fd_out], 'remove %i (%i) blocks because they contain no annotated genes\n', length(rm_idx), length(blocks));
  blocks(rm_idx) = [];

  for i=1:length(blocks)
    blocks(i).chr_num = strmatch(blocks(i).chr, genome_info.contig_names, 'exact') ;
    assert(length(blocks(i).chr_num)==1) ;
  end ;
  rm_idx = find_invalid_blocks(blocks, genome_info);
  Tfprintf([1 fd_out], 'remove %i (%i) blocks because their boundaries are invalid\n', length(rm_idx), length(blocks));
  blocks(rm_idx) = [];

	if 0
		% create blocks on the oposite strand to get good negative examples for genes
		blocks_rev = blocks;
		for j = 1:length(blocks)
			blocks_rev(j).gene_ids = [];
			if blocks(j).strand=='+'
				blocks_rev(j).strand='-';
			else
				blocks_rev(j).strand='+';
			end
		end
		blocks = [blocks, blocks_rev];
		blocks = blocks(randperm(length(blocks)));
	end

else
  %% create blocks from a genes structure that are patchy with
  %% many genes missing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  Tfprintf([1 fd_out], 'create blocks as regions arround genes\n')
  blocks = make_blocks_from_genes(genes, fn_config, opts.block_design_merge_sep, 0) ;
  for i=1:length(blocks), blocks(i).split = [1] ; end ;
end

Tfprintf([1 fd_out], 'Saving original blocks\n') ;
if exist('chr_num', 'var')
	save('-v7', strrep(PAR.FN.input_lsl.fn_training_blocks, 'training_blocks', sprintf('blocks_orig_chr%i',chr_num )), 'blocks') ;
else
	save('-v7', strrep(PAR.FN.input_lsl.fn_training_blocks, 'training_blocks', 'blocks_orig'), 'blocks') ;
end

% infer UTRs if not annotated
%--------------------------------------------------------------
if opts.use_rna_seq_for_label_gen==-1
	%% do nothing
elseif opts.use_rna_seq_for_label_gen==0
	Tfprintf([1 fd_out],'Infer utrs using tss and cleave signals\n')
	genes = infer_utrs(blocks, genes, PAR);
else
	Tfprintf([1 fd_out],'Infer utrs using tss and cleave signals and tracks\n')
	j = opts.use_rna_seq_for_label_gen;% track number 
	genes = infer_utrs(blocks, genes, PAR, opts.track_functions{j}, opts.track_files{j});
end

Tfprintf([1 fd_out],'Generate segmentation and filter\n')
[blocks, genes, split] = gen_lsl_blocks(genes, blocks, PAR.FN.genome.fn_genome_config, PAR.model, '');


Tfprintf([1 fd_out], 'Saving genes\n') ;
save('-v7', strrep(PAR.FN.input_lsl.fn_training_blocks, 'blocks', 'genes'), 'genes', 'PAR') ;
clear genes ;



%%%%%%%%%%%%%%%%%%%%
% reduce number of blocks
%%%%%%%%%%%%%%%%%%%%
train_test_cnt = subsample_policy(length(blocks), 'gene');
train_test_cnt(1) = min(train_test_cnt(1), opts.maxNumBlock) ;
maxval = sum(train_test_cnt);
merge_num = opts.block_design_merge_num; 
if exist('chr_num', 'var')
  some_more = 1;
else
  some_more = 1.5;
end
if length(blocks)/merge_num>maxval*some_more
  Tfprintf([1 fd_out],'Reduce total number of blocks to some reasonable number (from %i to ', length(blocks));
  blocks = blocks(1:floor(maxval*merge_num*some_more)) ;
  Tfprintf([1 fd_out], '%i)\n', length(blocks)) ;
end



%%%%%%%%%%%%%%%%%%%%
% add additional tracks to blocks
%%%%%%%%%%%%%%%%%%%%
if isfield(opts, 'subsample_reads')%&&opts.subsample_reads
	Tfprintf([1 fd_out],'Determine subsampling factors for blocks\n');
	for j = 1:length(blocks)
		x = rand(1);
		blocks(j).subsample_reads = x^opts.subsample_reads;
		%blocks(j).subsample_reads = x^2;
		%blocks(j).subsample_reads = x^3;
	end
end

for j = 1:length(opts.track_names)
  Tfprintf([1 fd_out], 'Adding track %s from file %s\n', opts.track_names{j}, opts.track_files{j});
  blocks = feval(opts.track_functions{j}, blocks, opts.track_files{j}, opts.track_params{j});
end
for j = 1:length(opts.segment_feature_names)
  Tfprintf([1 fd_out], 'Adding segment_feature %s from file %s\n', opts.segment_feature_names{j}, opts.segment_feature_files{j});
  blocks = feval(opts.segment_feature_functions{j}, blocks, opts.segment_feature_files{j}, genome_info, opts.segment_feature_params{j});
end

%if exist('chr_num', 'var')&&chr_num==16
%  bb = blocks(1:min(50, length(blocks)));
%  save('~/tmp/human_rep_mask_blocks', 'bb')
%  clear bb
%end

%%%%%%%%%%%%%%%%%%%%
% generate mask for loss such that 
% no loss will be incured at positions 
% where transcripts disagree
%%%%%%%%%%%%%%%%%%%%
if opts.use_loss_mask==1
  blocks = generate_loss_mask(blocks, PAR.model);
  rm_idx = [];
  for j = 1:length(blocks)
	if mean(blocks(j).loss_mask)<0.4
		rm_idx = [rm_idx j];
	end
  end
  Tfprintf([1 fd_out], 'Filter blocks according to loss mask: remove %i (%i)\n', length(rm_idx), length(blocks));
  blocks(rm_idx) = [];
end
if opts.scale_loss_mask==1
  blocks = scale_loss_according_to_confirmation(blocks);
end

if isempty(blocks)
  return
end

%%%%%%%%%%%%%%%%%%%%
% filter blocks according to intron lists data
% if available
%%%%%%%%%%%%%%%%%%%%
if 1
	% reorder the truth field such that the first field 
	% has maximal expression
	if opts.use_rna_seq_for_label_gen>0
		[blocks max_cov] = sort_truth(blocks);

		%cutoff = 1e-2;%prctile(max_cov, 25);
		cutoff = prctile(max_cov, 50);
		%cutoff = my_prctile(max_cov, 80);
		Tfprintf([1 fd_out], 'Removing %i (%i) blocks with coverage less than %i\n', sum(max_cov<cutoff), length(blocks), cutoff);
		blocks = blocks(max_cov>=cutoff);	
	end
else
	idx = filter_blocks_according_to_intron_confirmation(blocks);
	Tfprintf([1 fd_out], 'Reduce total number of blocks from %i to %i according to intron confirmation\n', length(blocks), length(idx));
	blocks = blocks(idx);
end

if isfield(opts, 'exon_map_file')
	rm_idx = filter_according_to_exon_maps(blocks, PAR.model, opts.exon_map_file);
	%rm_idx = filter_according_to_exon_maps(blocks, PAR.model, PAR.FN.genome.fn_genome_config);
	Tfprintf([1 fd_out], 'remove %i (%i) blocks because there are exons (possibly noncoding) within intergenic regions\n', length(rm_idx), length(blocks))
	blocks(rm_idx) = [];
end

if length(blocks)/merge_num>maxval*1.1
  Tfprintf([1 fd_out], 'Reduce total number of blocks to the number really needed (from %i to ', length(blocks));
  blocks = blocks(1:floor(maxval*merge_num*1.1)) ;
  Tfprintf([1 fd_out], '%i)\n', length(blocks)) ;
end

%%%%%%%%%%%%%%%%%%%%
% add signals to blocks 
%%%%%%%%%%%%%%%%%%%%
Tfprintf([1 fd_out], 'Load signal predictions:\n');
blocks = make_prediction_blocks(PAR, blocks, PAR.model, 0);
%Tfprintf([1 fd_out], 'done\n');




%%%%%%%%%%%%%%%%%%%%
% add content predictions to blocks
%%%%%%%%%%%%%%%%%%%%
Tfprintf([1 fd_out], 'Load predictions of content sensors:\n');
%[exon_dir, intron_dir, intergenic_dir, utr5_dir, utr3_dir] = ...
%    sort_contents(exon_dir, intron_dir, intergenic_dir, utr5_dir, utr3_dir) ;
cont_pred_dirs = {exon_dir, intron_dir, intergenic_dir, utr5_dir, utr3_dir}; 
blocks = add_content_predictions_to_blocks(blocks, cont_pred_dirs);
%Tfprintf([1 fd_out], 'done\n'); 

%%%%%%%%%%%%%%%%%%%%
% merge blocks for lsl training
%%%%%%%%%%%%%%%%%%%%

%blocks_orig = blocks ;
if merge_num>1
	Tfprintf([1 fd_out], 'Merge single-gene blocks \n')
	blocks = merge_blocks(blocks, PAR.model, merge_num);
	for i=1:length(blocks), blocks(i).split = [1] ; end ;
end

%%%%%%%%%%%%%%%%%%%%
% obtain signal predictions for lsl training
%%%%%%%%%%%%%%%%%%%%
Tfprintf([1 fd_out], 'Create feature matrix for dynamic programs:\n');
%Tfprintf([1 fd_out], 'generate feature matrix for blocks \n')
blocks = gen_features(blocks, PAR.model);
%Tfprintf([1 fd_out], 'done\n');


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
Tfprintf([1 fd_out], 'Find nearest predicted tss and cleave position for the annotated ones and move them\n')
rm_idx = [];
for b = 1:length(blocks)
  for t=1:length(blocks(b).truth)
   	new_truth = correct_cleave_and_tss_pos(blocks(b).truth(t), blocks(b), PAR.model);
	if ~isempty(new_truth)
		new_truth = shift_tss_and_cleave_to_best_prediction(new_truth, blocks(b), PAR.model);
	end
	if isfield(PAR.model.use.segments, 'rna_seq_polya')&&PAR.model.use.segments.rna_seq_polya
    	new_truth = move_cleave_pos(new_truth, blocks(b), PAR.model);
	end
    if isempty(new_truth), 
      rm_idx = [rm_idx b];
      continue ;
    end
    blocks(b).truth(t) = new_truth ;
  end
end
Tfprintf([1 fd_out], 'Remove %i of %i blocks because no tss or cleave position was found\n', length(rm_idx), length(blocks));
blocks(rm_idx) = [];

for xx = 1:length(blocks), emp(xx) = size(blocks(xx).truth(1).segments, 1)==1; end
if any(emp)
	fprintf('remove %i blocks because they contain no genes\n', sum(emp));
	blocks(find(emp)) = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%
% create paths from segmentation
%%%%%%%%%%%%%%%%%%%%%%%%
Tfprintf([1 fd_out], 'Create label sequence for dynamic programm:\n');
%Tfprintf([1 fd_out], 'create paths from segmentations\n')
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
PAR.model = prepare_model(PAR,blocks) ;

f = fieldnames(PAR.model.boundaries.lengths); 
for j=1:length(f), 
  %if strcmp(f{j},'intergenic_long'), continue, end
  if PAR.model.lengths_range.(f{j})(2) < exp(PAR.model.boundaries.lengths.(f{j})(end-1))+1
    PAR.model.lengths_range.(f{j})(2) = exp(PAR.model.boundaries.lengths.(f{j})(end-1))+1; 
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%
% check consistency of 'true' path with model
%%%%%%%%%%%%%%%%%%%%%%%%
Tfprintf([1 fd_out], 'Check consistency of label sequence with gene model: ')
ok_idx = check_block_model_consistency(blocks,PAR.model);
Tfprintf([1 fd_out], '\n%i out of %i blocks have label sequences that are not consistent with our gene model\n', length(blocks)-length(ok_idx),length(blocks));
blocks = blocks(ok_idx);

for i=1:length(blocks),
  blocks(i).id = i ;
end ;
Tfprintf([1 fd_out], 'Done.\n');


%%%%%%%%%%%%%%%%%%%%%%%%
% compute weights for the true paths
%%%%%%%%%%%%%%%%%%%%%%%%

Tfprintf([1 fd_out], 'Obtain features for the label sequences:\n');
blocks = initialize_blocks(blocks, PAR.model, 1, PAR.FN);
%Tfprintf([1 fd_out], 'done\n');



%%%%%%%%%%%%%%%%%%%%%%%%
% split blocks in training and validation examples
%%%%%%%%%%%%%%%%%%%%%%%%

if exist('chr_num', 'var')
	xx = whos('blocks')
	maxsize = 1.9e9;
	if 0%xx.bytes>maxsize
		new_num = floor(length(blocks)*maxsize/xx.bytes);
		Tfprintf([1 fd_out], 'reduce num of blocks from %i to %i to be able to save them\n', length(blocks), new_num);
		blocks = blocks(1:new_num);
	end
end

if isfield(blocks, 'train_val_split')
	idx_train = find([blocks.train_val_split]==1);
	idx_val = find([blocks.train_val_split]==-1);
	num_train = length(idx_train);
	num_test = length(idx_val);
	blocks = [blocks(idx_train(randperm(num_train))), blocks(idx_val(randperm(num_test)))];
else
	Tfprintf([1 fd_out], 'Generated %i usable blocks (each with 3 genes) from genome annotation\n\n', length(blocks)) ;
	train_test_cnt = subsample_policy(length(blocks), 'gene');
	num_train = min(train_test_cnt(1), opts.maxNumBlock) ;
	num_test = train_test_cnt(2);
end

Tfprintf([1 fd_out], '\n\n-------------------------------------------------------------------------------------\n') ;
blocks_train = blocks(1:num_train) ;
Tfprintf([1 fd_out], 'training on %i blocks\n', num_train) ;
if num_test>0
  Tfprintf([1 fd_out], 'performing validation on %i different blocks\n', num_test) ;
  blocks_val = blocks(num_train+1:num_train+num_test) ;
else
  Tfprintf([1 fd_out], 'performing validation on %i training blocks (not enough examples)\n', length(blocks)) ;
  blocks_val = blocks;
end
Tfprintf([1 fd_out], '-------------------------------------------------------------------------------------\n\n') ;

Tfprintf([1 fd_out], 'saving blocks') ;
if exist('chr_num', 'var')
  fn_train_blocks = strrep(PAR.FN.input_lsl.fn_training_blocks, 'blocks', sprintf('blocks_%i', chr_num));
  save_struct(fn_train_blocks, blocks_train, 'blocks');
  fn_val_blocks = strrep(PAR.FN.input_lsl.fn_val_blocks, 'blocks', sprintf('blocks_%i', chr_num));
  save_struct(fn_val_blocks, blocks_val, 'blocks');
else
  save_struct(PAR.FN.input_lsl.fn_training_blocks, blocks_train, 'blocks') ;
  save_struct(PAR.FN.input_lsl.fn_val_blocks, blocks_val, 'blocks') ;

  %save_append(PAR.FN.input_lsl.fn_training_blocks, 0, 'blocks', blocks_train) ;
  %save_append(PAR.FN.input_lsl.fn_val_blocks, 0, 'blocks', blocks_val) ;
end
blocks_split{1} = 1:length(blocks_train) ;
save('-v7', PAR.FN.input_lsl.fn_training_split, 'blocks_split') ;

%PAR_predict = PAR ;
%PAR_predict.fn_pred = [output_dir '/predict/'] ;
%unix(sprintf('rm -rf %s; mkdir -p %s', PAR_predict.fn_pred, PAR_predict.fn_pred)) ;
%PAR_predict.fn_blocks = [PAR_predict.fn_pred '/blocks.mat'] ;
%save_append(PAR_predict.fn_blocks, 0, 'blocks', blocks_val) ;


PAR.LSL.method.add_lin_feat = 0 ;
PAR.LSL.method.exm_per_solve = 500 ;
%PAR.RPROC.options.immediately_bg = 1 ;
if length(blocks_train)>500,
  PAR.RPROC.exm_per_batch = 5 ;
else
  PAR.RPROC.exm_per_batch = 1 ;
end ;

PAR.RPROC.options.addpaths = { fileparts(which('sg'))};
PAR.RPROC.options.envstr = shogun_envstr ;

clear blocks blocks_val blocks_train ;

PAR.submit = 1;

PAR.debug_mode = 0 ;

save('-v7', strrep(PAR.FN.input_lsl.fn_training_blocks, 'blocks', 'PAR'),'PAR') ;

fclose(fd_out) ;


 
return
function rm_idx = find_empty_blocks(blocks)
	rm_idx = [];
	for j = 1:length(blocks),
		if isempty(blocks(j).gene_ids)
			rm_idx = [rm_idx j];
		end
	end
return

function rm_idx = find_suspicious_genes(genes);
%% I consider genes suspicious where the minimal transcript start or 
%% maximal transcript stop are not equal to the gene start or gene 
%% stop, respectively. It is likely that there were more transcripts 
%% that were filtered out by some previous steps or couldn't 
%% be parsed from the gff
	rm_idx = [];
	for j = 1:length(genes)
    	min_start = inf;
	    max_stop = -inf;
    	for k = 1:length(genes(j).exons)
        	min_start = min([min_start min(genes(j).exons{k})]);
	        max_stop = max([max_stop max(genes(j).exons{k})]);
    	end
	    diff_start = min_start-genes(j).start;
    	diff_stop = genes(j).stop-max_stop;
		if diff_start>0|diff_stop>0
			rm_idx = [rm_idx j];        	
		end
	end


return
