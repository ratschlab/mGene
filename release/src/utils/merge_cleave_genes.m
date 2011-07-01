function merge_cleave_genes(organism, experiment, offsets)

fname = 'genes.mat';

cnt = 0;
for c = offsets
	cnt = cnt+1;
	bdir = get_pred_dir(organism, experiment, c);
	load(sprintf('%s/%s', bdir, fname), 'genes');
	%genes = genes([genes.chr_num]==1);
	%genes = closed_to_half_open(genes);
	if cnt==1
		g.genes = genes;
	else
		g(end+1).genes = genes;
	end
	g(end).new = zeros(1, length(g(end).genes));
end

other = length(g);
for j = 2:length(g)
	idx = find_overlapping_regions(g(j).genes, g(1).genes);
	g(j).new = ~ismember(1:length(g(j).genes), idx(:, 1));
	[l1, l2] = compute_overlap_lists(idx, length(g(j).genes), length(g(1).genes));
	g(j).map_split1 = zeros(1, length(l1));
	for k = 1:length(l1)
		if length(l1{k})==1
			g(j).map_split1(k) = length(l2{l1{k}});
		end
	end	
	idx = find_overlapping_regions(g(j).genes, g(j-1).genes);
	[l1, l2] = compute_overlap_lists(idx, length(g(j).genes), length(g(j-1).genes));
	g(j).map_split2 = zeros(1, length(l1));
	for k = 1:length(l1)
		if length(l1{k})==1 
			g(j).map_split2(k) = length(l2{l1{k}});
		end
	end	
end
for j = 1:length(g)-1
	idx = find_overlapping_regions(g(j).genes, g(other).genes);
	[l1, l2] = compute_overlap_lists(idx, length(g(j).genes), length(g(other).genes));
	g(j).map_merge = zeros(1, length(l1));
	for k = 1:length(l1)
		g(j).map_merge(k) = length(l1{k});
	end	
end

%%% merge
genes = g(1).genes;
idx = find_overlapping_regions(genes, g(other).genes);
[l1, l2] = compute_overlap_lists(idx, length(genes), length(g(other).genes));
genes_add = [];
rm_idx = [];
for j = 1:length(genes)
	if length(l1{j})>2
		rm_idx = [rm_idx j];
		split_genes = g(other).genes(l1{j});
		take_idx = find(g(other).map_split2(l1{j})==1);
		genes_add = [genes_add split_genes(take_idx)];
	end
end

exp_merge = [experiment '_cleave_merge'];
save_dir = get_pred_dir(organism, exp_merge, 0)

genes(rm_idx) = [];
nice_mkdir(save_dir)
save(sprintf('%sgenes.mat', save_dir), 'genes')

%genes = [genes genes_add];

return

