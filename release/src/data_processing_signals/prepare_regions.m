function [trivial_regions,blocks,training,test] = prepare_regions(PAR);
%  [trivial_regions,training_split,training_blocks,test_split,test_blocks] = prepare_regions(PAR,genes);


% generate trivial regions
if fexist(PAR.FN.output.fn_trivial_regions)
  fprintf(1,'loading trivial regions from file %s\n',PAR.FN.output.fn_trivial_regions)
  load(PAR.FN.output.fn_trivial_regions,'trivial_regions');
  assert(isequal(trivial_regions(1).config,PAR.FN.input.fn_genome_config));
else
  fprintf(1,'initializing trivial regions\n')
  trivial_regions = init_regions(PAR);
  save(PAR.FN.output.fn_trivial_regions,'trivial_regions');
end

% generate training and test blocks and splits, 
% considering training and test region
if fexist(PAR.FN.output.fn_training_blocks)& fexist(PAR.FN.output.fn_test_blocks)
  Old = load(PAR.FN.output.fn_training_blocks, 'PAR');
  assert(isequal(Old.PAR.FN.output.fn_trivial_regions,PAR.FN.output.fn_trivial_regions));
  assert(isequal(Old.PAR.FN.input.fn_genome_config,PAR.FN.input.fn_genome_config));
  diff_fields = compare_PAR(PAR,Old.PAR) ;
  if ismember('regions',diff_fields)
    warning('regions exist, Parameters are different ')
    %fprintf('regions exist, Parameters are different ')
  else
    fprintf(1,'loading split and blocks from file %s \n',PAR.FN.output.fn_training_blocks)
    training = load(PAR.FN.output.fn_training_blocks)
    assert(length(fieldnames(training.split))==PAR.SETs.num_splits)
  end
  Old = load(PAR.FN.output.fn_test_blocks, 'PAR');
  assert(isequal(Old.PAR.FN.output.fn_trivial_regions,PAR.FN.output.fn_trivial_regions));
  assert(isequal(Old.PAR.FN.input.fn_genome_config,PAR.FN.input.fn_genome_config));
  diff_fields = compare_PAR(PAR,Old.PAR) ;
  if ismember('regions',diff_fields)
    warning('regions exist, Parameters are different ')
    %fprintf('regions exist, Parameters are different ')
  else
    fprintf(1,'loading split and blocks from file %s \n',PAR.FN.output.fn_test_blocks)
    test = load(PAR.FN.output.fn_test_blocks)
    assert(length(fieldnames(test.split))==PAR.SETs.num_splits)
  end
else
  if isfield(PAR.FN.output,'fn_genes_merged') && fexist(PAR.FN.output.fn_genes_merged)
    load(PAR.FN.output.fn_genes_merged,'genes');
  else
    genes=[];
  end
  [training, test] = assign_subsets(trivial_regions,PAR.FN.input.fn_trainings_regions,PAR.FN.input.fn_test_regions,genes,PAR.SETs,PAR.regions,PAR.tasks)
end

if 1
  if isfield(PAR.FN.output,'fn_genes_merged') && fexist(PAR.FN.output.fn_genes_merged)
    load(PAR.FN.output.fn_genes_merged,'genes');
  else
    genes=[];
  end
  load(PAR.FN.input.fn_trainings_regions,'training_regions');
  load(PAR.FN.input.fn_test_regions,'test_regions');
  load(PAR.FN.output.fn_trivial_regions,'trivial_regions');
  if  ~exist('training')
    training = load(PAR.FN.output.fn_training_blocks);
  end
  if  ~exist('test')
    test = load(PAR.FN.output.fn_test_blocks);
  end
  % load(PAR.FN.output.fn_blocks,'blocks');
  blocks =test.blocks;
  check_and_plot_splits(trivial_regions(1:min(10,length(trivial_regions))),training_regions,test_regions,genes,blocks,training, test)
end
  
