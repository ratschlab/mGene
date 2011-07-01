




starts1 = [100 101 10000]; 
stops1 = [101 101 10003];

starts2 = [10 99 101 10004]; 
stops2 = [1000 100 102 10007];

[a,b] = interval_overlap(starts1, stops1, starts2, stops2)

load('/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/mouse/label_gen_alt_long/blocks_all_short.mat', 'blocks');

blocks = blocks([blocks.chr_num]==10);

[tmp, sort_idx] = sort([blocks.start]);
blocks = blocks(sort_idx);
%blocks = blocks(76:80);

tic; idx1 = find_overlapping_regions_slow(blocks, []); toc
tic; idx2 = find_overlapping_regions(blocks, [] ); toc
setdiff(idx1, idx2, 'rows')
setdiff(idx2, idx1, 'rows')
tic; idx1 = find_overlapping_regions_slow(blocks, blocks); toc
tic; idx2 = find_overlapping_regions(blocks, blocks ); toc
setdiff(idx1, idx2, 'rows')
setdiff(idx2, idx1, 'rows')

