function PAR = content_train(cpf_label_fname, cpf_label_dir, content_name, genome_config_dir, anno_dir, classifier_info, classifier_dir, run_locally_par, options)
% PAR = content_train(cpf_label_fname, cpf_label_dir, content_name, genome_config_dir, classifier_info, classifier_dir, run_locally)

paths

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started content_train at %s\n', timedate) ;

fprintf('Results are written to: %s\n\n', classifier_info) ;

% make sure no function uses cplex
global NO_CPLEX
NO_CPLEX = 1 ;

genome_config = [genome_config_dir '/genome.config'] ;
genome_info = init_genome(genome_config) ;

if ~exist('run_locally_par', 'var'),
  run_locally_par = -1 ;
end ;

%run_locally_par=0 ;

if run_locally_par==-1,
  run_locally.pos = 1 ; % rproc_policy('content_train:pos',  genome_info) ; % not implemented 
  run_locally.train = rproc_policy('content_train:train',  genome_info) ;
  run_locally.eval = rproc_policy('content_train:eval',  genome_info) ;
  run_locally.predict = rproc_policy('content_train:predict',  genome_info) ; 
  run_locally.save = 1;%rproc_policy('content_train:save',  genome_info) ; % not implemented
else
  run_locally.pos = 1 ;%run_locally_par ; % not implemented 
  run_locally.train = run_locally_par ;
  run_locally.eval = run_locally_par ; 
  run_locally.predict = run_locally_par ; 
  run_locally.save = 1 ; % run_locally_par ; % not implementede
end ;

engine = determine_engine() ;
if isequal(engine, 'octave'),
  struct_levels_to_print(1) ;
end ;
run_locally
if isequal(engine, 'octave'),
  struct_levels_to_print(0) ;
end ;


fid = 1;
  
disp('------------------------------------') ;
fprintf('Processing content: %s\n', content_name );
disp('------------------------------------') ;


disp('------------------------------------') ;
disp('Step 1: Setting up datastructures...') ;
disp('------------------------------------') ;



fprintf(fid,'\n\nstart training for sensor %s\n\n',content_name);

galaxy_PAR.dir = classifier_dir ;
galaxy_PAR.genome_config = genome_config ;
PAR = create_PAR('exp_galaxy', galaxy_PAR);

PAR.tasks.contents.prepare_labels = 1;
PAR.tasks.contents.train = 1;
PAR.tasks.contents.modelsel = 1;
PAR.tasks.contents.learn_plifs = 0;
PAR.tasks.contents.get_train_error = 1;
PAR.tasks.contents.pred_on_genome = 0;
PAR.tasks.contents.run_locally = run_locally ;

%PAR.RPROC.MEMREQ = 2000;

% keyboard
create_directories(PAR.FN,PAR.organism,0,content_name) ;

genome_config = PAR.FN.input.fn_genome_config;
organism = PAR.organism.name;

if ~isempty(anno_dir),

  genes = load_genes(sprintf('%s/genes.mat', anno_dir)) ;

  fprintf(fid,'\n\nloaded %i annotated genes \n\n', length(genes));
else
  genes = [] ;
end ;

fprintf(fid,'init_regions\n')
trivial_regions = init_regions(genome_config);

genome_info =  init_genome(genome_config) ;

% it is necessary to do this again here since the genes are not saved 
% by anno2contentlabel, where the filtering was already done
			  genes = filter_invalid_genes(genes, genome_info);

[training, test] = assign_subsets(trivial_regions, PAR.FN.input.fn_trainings_regions,PAR.FN.input.fn_test_regions,genes,PAR.SETs,PAR.regions,PAR.tasks.contents) ;

% training_regions = init_regions(genome_config);
% save(PAR.FN.input.fn_trainings_regions,'training_regions');

% test_regions=training_regions(1);
% test_regions.stop = 0;
% test_regions.start = 0;
% save(PAR.FN.input.fn_test_regions, 'test_regions');

% trivial_regions = training_regions;
% save(PAR.FN.output.fn_trivial_regions,'trivial_regions');
% clear trivial_regions

fprintf('Saving blocks:\n\n') ;
blocks = training.blocks;
split = training.split;
save_struct(PAR.FN.output.fn_training_blocks, blocks, 'blocks') ;
save_append(PAR.FN.output.fn_training_blocks, 1, 'PAR', PAR, 'split', split);
fprintf('Done.\n\n') ;


disp('-----------------------------------------') ;
disp('Step 2: Reading and filter label files...') ;
disp('-----------------------------------------') ;

cpf_label_mat_fname = sprintf('%s/label_cpf.mat', cpf_label_dir) ;
if fexist(cpf_label_mat_fname),
  disp('Reading binary CPF label files') ;
  [chr, strand, pos, score, content_names, score_names, chr_dic, content_names_dic, score_names_dic, success] = read_contpred_bin(cpf_label_mat_fname) ;
  assert(success) ;

  %[chr_, strand_, pos_, score_, content_names_, score_names_, chr_dic_, content_names_dic_, score_names_dic_]=read_contpred(cpf_label_fname) ;
  %keyboard
  %assert(isequal(chr,chr_)) ;
  %assert(isequal(strand,strand_)) ;
  %assert(isequal(pos,pos_)) ;
  %assert(isequal(score,score_)) ;
  %assert(isequal(score_names,score_names_)) ;
  %assert(isequal(content_names_dic,content_names_dic)) ;
  %assert(isequal(score_names_dic,score_names_dic)) ;
else
  disp('Reading CPF label files') ;
  [chr, strand, pos, score, content_names, score_names, chr_dic, content_names_dic, score_names_dic]=read_sigpred(cpf_label_fname) ;
end ;

if isequal(content_name, 'auto') && length(content_names_dic)==1,
  content_name = content_names_dic{1} ;
end ;

% filter out spurious entries
idx_ = strmatch(content_name, content_names_dic, 'exact') ;
idx = find(content_names == idx_) ;
if length(idx)~=length(content_names),
  warning('dropping %i/%i entries from other contents', length(content_names)-length(idx), length(content_names)) ;
  chr = chr(idx) ;
  strand = strand(idx) ;
  pos = pos(idx,:) ;
  score = score(idx) ;
  content_names = content_names(idx) ;
  score_names = score_names(idx) ;
end 

idx_ = strmatch('label', score_names_dic, 'exact') ;
idx = find(score_names == idx_) ;
if length(idx)~=length(score_names),
  warning('dropping %i/%i non-label entries', length(score_names)-length(idx), length(score_names)) ;
  chr = chr(idx) ;
  strand = strand(idx) ;
  pos = pos(idx,:) ;
  score = score(idx) ;
  content_names = content_names(idx) ;
  score_names = score_names(idx) ;
end 

idx = find(abs(abs(score)-1)<=1e-3) ;
if length(idx)~=length(score),
  warning('dropping %i/%i entries with label not +/-1', length(score)-length(idx), length(score)) ;
  chr = chr(idx) ;
  strand = strand(idx) ;
  pos = pos(idx,:) ;
  score = score(idx) ;
  content_names = content_names(idx) ;
  score_names = score_names(idx) ;
end 

unix(sprintf('mkdir -p %s', PAR.FN.input_cont.(content_name).fn_candsites)) ;
num_pos = 0 ;
num_neg = 0 ;
for i=1:length(genome_info.contig_names),
  for s='+-',
    idx_ = strmatch(lower(genome_info.contig_names{i}), lower(chr_dic), 'exact') ;
    if ~isempty(idx_)
      idx = find(chr==idx_) ;
      idx2 = find(strand(idx)==s) ;
      idx = idx(idx2) ;
    else
      idx = [];
    end
    save_append(sprintf('%s/contig_%i%c.mat', PAR.FN.input_cont.(content_name).fn_candsites, i, s), 0, 'pos', pos(idx,:), 'label', score(idx)) ;
    num_pos = num_pos + sum(sign(score(idx))==1) ;
    num_neg = num_neg + sum(sign(score(idx))==-1) ;
  end ;
end 
fprintf('Number of pos: %i \n',num_pos) ;
fprintf('Number of neg: %i \n',num_neg) ;

fprintf('Done.\n\n') ;

unix(sprintf('mkdir -p %s', PAR.FN.output_cont.(content_name).fn_SVMs)) ;
unix(sprintf('mkdir -p %s', PAR.FN.output_cont.(content_name).fn_pred)) ;

disp('----------------------------') ;
disp('Step 3: Starting Training...') ;
disp('----------------------------') ;
PAR.Content_name = content_name ;
all_done = 0 ; iter = 0 ;
while ~all_done,
  fprintf('\n\n start masterscript_contents: iter:%i \n\n', iter)
  [PAR, jobinfo, all_done] = masterscript_contents(PAR, 1) ;
  iter = iter + 1 ;
  if iter>10, % this should not normally happen
    error('Computations aborted') ;
  end ;
end ;
fprintf('Done.\n\n') ;

disp('-------------------------------') ;
disp('Step 4: Writing output files...') ;
disp('-------------------------------') ;

PAR.Content_name = content_name ;
save([classifier_dir '/PAR.mat'], 'PAR', '-v7') ;

L=load(PAR.FN.output_cont.(content_name).fn_models) ;
num_pos = sum(score>0) ;
num_neg = sum(score<0) ;

fi = fopen(classifier_info, 'w+') ;
fprintf(fi, 'Trained "%s" classifier\n', content_name) ;
fprintf(fi, ' * based on %i labeled examples (%i positive, %i negative)\n', num_pos+num_neg, num_pos, num_neg) ;
fprintf(fi, ' * using %i-fold cross-validation for model-selection from %i models (inner cv loop)\n', PAR.SETs.num_splits, size(L.MODELS,2)) ;
fprintf(fi, ' * using %i-fold cross-validation for obtaining unbiased predictions (outer cv loop)\n\n', PAR.SETs.num_splits) ;
fprintf(fi, 'Performance\n') ;
fprintf(fi, ' * Average area under ROC curve on test splits: %1.3f\n', mean(L.RSV_test)) ;
fprintf(fi, ' * Average area under PRC curve on test splits: %1.3f\n', mean(L.RFV_test)) ;

fclose(fi) ;

fprintf('Done.\n\n') ;


[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('finished content_train at %s\n', timedate) ;

return
