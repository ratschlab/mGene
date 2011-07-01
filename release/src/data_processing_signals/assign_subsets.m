function [train,test] = assign_subsets(trivial_regions,fn_trainings_regions,fn_test_regions,genes,SETs,regions_par,tasks)
% [train,test] = assign_subsets(trivial_regions,fn_trainings_regions,fn_test_regions,genes,SETs,regions_par,tasks)


% label_source
% if fexist(fn_genes)&(label_source.use_confirmed|label_source.use_anno)
  
if ~isempty(fn_trainings_regions)&fexist(fn_trainings_regions)
  assert(ismember('training_regions',who('-file',fn_trainings_regions)))
  load(fn_trainings_regions,'training_regions');
else
  training_regions = trivial_regions;
  fprintf('no training region specified, training on whole genome except test region\n')
end

if ~isempty(fn_test_regions)&fexist(fn_test_regions)
  assert(ismember('test_regions',who('-file',fn_test_regions)))
  load(fn_test_regions,'test_regions');
else
  test_regions = trivial_regions(1);
  test_regions.stop = 0;
  test_regions.start = 0;
  fprintf('no test region specified, training on whole genome\n')
end
fprintf('make sure there is no overlap between training and test regions\n')  
if ~(tasks.train)
  training_regions = assert_noRegionOverlap(test_regions,training_regions);
else
  warning('region overlap not checked!! Training and Test region could overlap!\n')
end
rand('seed',1234);
idx = [];

fprintf('Splitting testing regions\n') ;
regions_par.offset = regions_par.offset_test;
test.blocks = split_regions(test_regions, [], regions_par, 0);
fprintf('Done\n') ;

fprintf('Splitting training regions\n') ;
regions_par.offset = regions_par.offset_train;
train.blocks = split_regions(training_regions, genes, regions_par,0);
fprintf('Done\n') ;

if ~isempty(test.blocks)
  test.blocks = rmfield(test.blocks, 'left_dist');
  test.blocks = rmfield(test.blocks, 'right_dist');
  % test.blocks = sort_regions(test.blocks) ;
  % for j=1:length(test.blocks)
  %   test.blocks(j).id = j;
  % end
  [test.blocks, test.split] = make_split(test.blocks,SETS)
else
  test.split = [];
end

% merge training blocks and generate training split
if ~isempty(train.blocks)
  % train.blocks = sort_regions(train.blocks) ; 
  % for j=1:length(train.blocks)
  %   train.blocks(j).id = j;
  % end
  if regions_par.merge_confirmed
    train.blocks = construct_merged_blocks(train.blocks,PAR);
  else
    train.blocks = rmfield(train.blocks,'left_dist');
    train.blocks = rmfield(train.blocks,'right_dist');
  end
  if regions_par.cluster_paralogs
    split_idx = splits_by_paralogs(train.blocks,SETs.num_splits);
    for s = 1:SETs.num_splits
      idx = setfield(idx,sprintf('split_%i',s),split_idx{s});
      error('something wron here; check implementation')
      train.blocks(idx.(sprintf('split_%i',s))) = s;
    end
  else
    [train.blocks, train.split] = make_split(train.blocks,SETs);
  end
else
  train.split = [];
end

return

function [blocks, split] = make_split(blocks,SETs)

blocks_per_split = ceil(length(blocks)/SETs.num_splits);
block_idx  = randperm(length(blocks));
for s = 1:SETs.num_splits
  IDX = (s-1)*blocks_per_split+[1:blocks_per_split];
  IDX(IDX>length(blocks)) = [];
  IDX = block_idx(IDX);
  split.(sprintf('split_%i',s)) = IDX;
% [blocks(IDX).split] = deal(s); % problematic piece of code
  for i=IDX
    blocks(i).split = s;
  end ;
end

y = [];
for s = 1:SETs.num_splits
  y = [y getfield(split,sprintf('split_%i',s))];
end
assert(length(y)==length(unique(y)))
