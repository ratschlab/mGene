function [PAR, jobinfo, all_done] = signal_train(spf_label_fname, spf_label_dir, signal_name, genome_config_dir, classifier_info, classifier_dir, run_locally_par, options)
% PAR = signal_train(spf_label_fname, spf_label_dir, signal, genome_config_dir, classifier_info, classifier_dir, run_locally,options)

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started signal_train at %s\n', timedate) ;

fprintf('Results are written to: %s\n', classifier_info) ;

% make sure no function uses cplex
global NO_CPLEX
NO_CPLEX = 0 ;

genome_config = [genome_config_dir '/genome.config'] ;
genome_info = init_genome(genome_config) ;

if ~exist('run_locally_par', 'var'),
  run_locally_par = -1 ;
end ;

if run_locally_par==-1,
  run_locally.pos = rproc_policy('signal_train:pos',  genome_info) ;
  run_locally.train = rproc_policy('signal_train:train',  genome_info) ;
  run_locally.eval = rproc_policy('signal_train:eval',  genome_info) ;
  run_locally.predict = rproc_policy('signal_train:predict',  genome_info) ;
  run_locally.conf = rproc_policy('signal_train:conf',  genome_info) ;
  run_locally.save = rproc_policy('signal_train:save',  genome_info) ;
else
  run_locally.pos = run_locally_par ;
  run_locally.train = run_locally_par ;
  run_locally.eval = run_locally_par ;
  run_locally.predict = run_locally_par ;
  run_locally.conf = run_locally_par ;
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

if exist('options','var')
  PAR.RPROC.options = options;
  fprintf('signal_train:::run_locally = %i', run_locally_par)
end

keep_all_files = 1;

fid = 1;
  
disp('------------------------------------') ;
fprintf('Processing signal: %s\n', signal_name );
disp('------------------------------------') ;


disp('------------------------------------') ;
disp('Step 1: Setting up datastructures...') ;
disp('------------------------------------') ;



fprintf(fid,'\n\nstart training for sensor %s\n\n',signal_name);

galaxy_PAR.dir = classifier_dir ;
galaxy_PAR.genome_config = genome_config ;
PAR = create_PAR('exp_galaxy', galaxy_PAR);

PAR.tasks.signals.prepare_labels = 1;
PAR.tasks.signals.train = 1;
PAR.tasks.signals.learn_plifs = 1;
PAR.tasks.signals.get_train_error = 1;
PAR.tasks.signals.run_locally = run_locally ;


create_directories(PAR.FN,PAR.organism,0,signal_name) ;

genome_config = PAR.FN.input.fn_genome_config;
organism = PAR.organism.name;

fprintf(fid,'initializing regions\n') ;
training_regions = init_regions(genome_config);
save(PAR.FN.input.fn_trainings_regions,'training_regions', '-v7');

test_regions=training_regions(1);
test_regions.stop = 0;
test_regions.start = 0;
save(PAR.FN.input.fn_test_regions, 'test_regions', '-v7');

trivial_regions = training_regions;
save(PAR.FN.output.fn_trivial_regions,'trivial_regions', '-v7');
clear trivial_regions

blocks=training_regions([]) ;
block_length=1000 ;
num_blocks = 0 ;
for i=1:length(training_regions)
  for p=1:block_length:training_regions(i).stop,
    num_blocks = num_blocks + 1;
  end ;
end ;
if num_blocks>10000,
  block_length= block_length*10 ; % 1000*round(num_blocks/100000) 
  if num_blocks>100000,
    block_length= block_length*10 ; % 1000*round(num_blocks/100000) 
  end ;
  num_blocks = 0 ;
  for i=1:length(training_regions)
    for p=1:block_length:training_regions(i).stop,
      num_blocks = num_blocks + 1;
    end ;
  end ;
end ;

fprintf('creating %i blocks\n', num_blocks) ;
blocks(num_blocks).start=[] ;
num_blocks=0 ;
for i=1:length(training_regions)
  for p=1:block_length:training_regions(i).stop,
    num_blocks = num_blocks + 1 ;
    blocks(num_blocks).start = p ;
    blocks(num_blocks).stop = min(training_regions(i).stop, p+block_length-1) ;
    blocks(num_blocks).chr = training_regions(i).chr ;
    blocks(num_blocks).chr_num = training_regions(i).chr_num ;
    blocks(num_blocks).strand = training_regions(i).strand ;
    blocks(num_blocks).id = num_blocks ;
    if blocks(num_blocks).strand=='+',
      blocks(num_blocks).offset = blocks(num_blocks).start - 1 ;
    else
      blocks(num_blocks).offset = blocks(num_blocks).stop + 1 ;
    end ;
    blocks(num_blocks).config = genome_config ;
    blocks(num_blocks).organism = organism ;
  end ;
end ;
rand('seed', 87543875848) ;
sp = ceil(rand(1,length(blocks))*PAR.SETs.num_splits) ;
split = struct ;
for i=1:PAR.SETs.num_splits,
  split.(sprintf('split_%i', i)) = find(sp==i) ;
end ;

fprintf('saving blocks\n') ;
save_struct(PAR.FN.output.fn_training_blocks, blocks, 'blocks') ;
save_append(PAR.FN.output.fn_training_blocks, 1, 'PAR', PAR, 'split', split);

fprintf('Done.\n\n') ;


disp('-----------------------------------------') ;
disp('Step 2: Reading and filter label files...') ;
disp('-----------------------------------------') ;

genome_info =  init_genome(genome_config) ;
spf_label_mat_fname = sprintf('%s/label_spf.mat', spf_label_dir) ;
if fexist(spf_label_mat_fname),
  disp('Reading binary SPF label files') ;
  [chr, strand, pos, score, signal_names, score_names, chr_dic, signal_names_dic, score_names_dic, success] = read_sigpred_bin(spf_label_mat_fname) ;
  assert(success) ;
elseif fexist(sprintf('%s/contig_1+.label', spf_label_dir))
  disp('Reading interval query label files') ;
  [chr, strand, pos, score] = read_sigpred_interval(spf_label_dir, genome_info) ;
  signal_names_dic = {signal_name};
  signal_names = ones(1, length(pos));
  score_names_dic = {'label'}; 
  score_names = ones(1, length(pos));
  chr_dic = genome_info.contig_names;
  
else
  disp('Reading SPF label files') ;
  [chr, strand, pos, score, signal_names, score_names, chr_dic, signal_names_dic, score_names_dic]=read_sigpred(spf_label_fname) ;
end ;

if isequal(signal_name, 'auto') && length(signal_names_dic)==1,
  signal_name = signal_names_dic{1} ;
end ;

% filter out spurious entries
idx_ = strmatch(signal_name, signal_names_dic, 'exact') ;
idx = find(signal_names == idx_) ;
if length(idx)~=length(signal_names),
  warning('dropping %i/%i entries from other signals', length(signal_names)-length(idx), length(signal_names)) ;
  chr = chr(idx) ;
  strand = strand(idx) ;
  pos = pos(idx) ;
  score = score(idx) ;
  signal_names = signal_names(idx) ;
  score_names = score_names(idx) ;
end 

idx_ = strmatch('label', score_names_dic, 'exact') ;
idx = find(score_names == idx_) ;
if length(idx)~=length(score_names),
  warning('dropping %i/%i non-label entries', length(score_names)-length(idx), length(score_names)) ;
  chr = chr(idx) ;
  strand = strand(idx) ;
  pos = pos(idx) ;
  score = score(idx) ;
  signal_names = signal_names(idx) ;
  score_names = score_names(idx) ;
end 

idx = find(abs(abs(score)-1)<=1e-3) ;
if length(idx)~=length(score),
  warning('dropping %i/%i entries with label not +/-1', length(score)-length(idx), length(score)) ;
  chr = chr(idx) ;
  strand = strand(idx) ;
  pos = pos(idx) ;
  score = score(idx) ;
  signal_names = signal_names(idx) ;
  score_names = score_names(idx) ;
end 

unix(sprintf('mkdir -p %s', PAR.FN.input_sig.(signal_name).fn_candsites)) ;
num_pos = 0 ;
num_neg = 0 ;
for i=1:length(genome_info.contig_names),
  for s='+-',
    idx_ = strmatch(lower(genome_info.contig_names{i}), lower(chr_dic), 'exact') ;
    if isempty(idx_),
      idx = [] ;
    else
      idx = find(chr==idx_) ;
    end ;
    idx2 = find(strand(idx)==s) ;
    idx = idx(idx2) ;
    save_append(sprintf('%s/contig_%i%c.mat', PAR.FN.input_sig.(signal_name).fn_candsites, i, s), 0, 'pos', pos(idx), 'label', score(idx)) ;
    num_pos = num_pos + sum(sign(score(idx))==1) ;
    num_neg = num_neg + sum(sign(score(idx))==-1) ;
  end ;
end 
num_pos
num_neg
if num_pos<50
   warning('number of positive examples is low!!')
end
if num_neg<50
   warning('number of negative examples is low!!')
end

fprintf('Done.\n\n') ;

unix(sprintf('mkdir -p %s', PAR.FN.output_sig.(signal_name).fn_SVMs)) ;
unix(sprintf('mkdir -p %s', PAR.FN.output_sig.(signal_name).fn_pred)) ;

disp('----------------------------') ;
disp('Step 3: Starting Training...') ;
disp('----------------------------') ;
PAR.Signal_name = signal_name ;
all_done = 0 ; iter = 0 ;

[PAR, jobinfo, all_done] = masterscript_sensors(PAR, 1,keep_all_files) ;

if all_done,
  fprintf('Done.\n\n') ;

  disp('-------------------------------') ;
  disp('Step 4: Writing output files...') ;
  disp('-------------------------------') ;
  
  PAR.Signal_name = signal_name ;
  save([classifier_dir '/PAR.mat'], 'PAR', '-v7') ;
  
  L=load(PAR.FN.output_sig.(signal_name).fn_models) ;
  num_pos = sum(score>0) ;
  num_neg = sum(score<0) ;
  
  fi = fopen(classifier_info, 'w+') ;
  fprintf(fi, 'Trained "%s" classifier\n', signal_name) ;
  fprintf(fi, ' * based on %i labeled examples (%i positive, %i negative)\n', num_pos+num_neg, num_pos, num_neg) ;
  fprintf(fi, ' * using %i-fold cross-validation for model-selection from %i models (inner cv loop)\n', PAR.SETs.num_splits, size(L.MODELS,2)) ;
  fprintf(fi, ' * using %i-fold cross-validation for obtaining unbiased predictions (outer cv loop)\n\n', PAR.SETs.num_splits) ;
  fprintf(fi, 'Performance\n') ;
  fprintf(fi, ' * Average area under ROC curve on test splits: %1.3f\n', mean(L.RSV_test)) ;
  fprintf(fi, ' * Average area under PRC curve on test splits: %1.3f\n', mean(L.RFV_test)) ;

  fclose(fi) ;
  
  fprintf('Done.\n\n') ;
end

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('finished signal_train at %s\n', timedate) ;

return
