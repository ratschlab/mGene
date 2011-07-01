function genes_from_predictions(organism, experiment, cleave_off, base_exp, base_organism)

if nargin<5
	base_organism = organism;
end
if nargin<4
	base_exp = experiment;
end

if exist(organism, 'dir')==7
	% reusing this function for mGeneToolbox 
	genome_info_dir = experiment;
	fn_config = fullfile(genome_info_dir , 'genome.config');
else
	[genome_info_dir, fn_config] = rgasp_genome_config_dir(organism);
end



basedir = get_pred_dir(organism, experiment, cleave_off)

fn_predicted_genes = sprintf('%s/genes.mat', basedir)
fn_predicted_genes_rquant = sprintf('%s/genes_half_open.mat', basedir)
fn_predicted_genes_date = sprintf('%s/genes_%s.mat', basedir, date)
fn_pred_blocks = sprintf('%s/prediction_blocks.mat', basedir)
have_predictions = 0;
if 0%fexist(fn_predicted_genes)
  return
end

if fexist(fn_pred_blocks)
	load(fn_pred_blocks, 'blocks');
else
	blocks = init_regions(fn_config);
end

genome_info = init_genome(fn_config);

basedir2 = get_pred_dir(base_organism, base_exp, 0);
%load(sprintf('/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/%s/lsl/%s/output/lsl/train/init.mat', base_organism, base_exp), 'PAR');
load(sprintf('%s/../lsl/train/init.mat', basedir2), 'PAR');

for j = 1:length(blocks)
	fn_pred = sprintf('%s/_viterbi_block%i.mat', basedir, j);
	if fexist(fn_pred)
		try
			load(fn_pred, 'pred');
			blocks(j).prediction = pred;
			have_predictions(j) = 1;
		catch
			fprintf('file %s not completely written\n', fn_pred);
			blocks(j).prediction.genes = {};
			have_predictions(j) = 0;
		end
	else
		blocks(j).prediction.genes = {};
		have_predictions(j) = 0;
	end
	clear pred
end


% determine number of genes
all_genes = {}; 
for j = 1:length(blocks), 
	all_genes = {all_genes{:}, blocks(j).prediction.genes{:}}; 
end

fprintf('%i (%i) blocks done; number of predicted genes: %i estimated num of genes for all blocks: %i\n', sum(have_predictions), length(blocks), length(all_genes) , length(all_genes)/sum(have_predictions)*length(blocks));
clear all_genes
%keyboard

if 1%all(have_predictions==1)

	for j = 1:length(blocks)-1
		if length(blocks(j).prediction.genes)==1&&isempty(blocks(j).prediction.genes{1})
			continue;
		end
		if blocks(j).chr_num==blocks(j+1).chr_num&&blocks(j).strand==blocks(j+1).strand
			if blocks(j).strand=='+'
				assert(blocks(j).stop>blocks(j+1).start)
				gene_map1 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
				gene_map2 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
				overlap = floor(length(gene_map1)/2);
	
				gidx = length(blocks(j).prediction.genes);
				while gidx>0
					gene = blocks(j).prediction.genes{gidx};
					global_gene = gene(:, 1:2) + blocks(j).start-1;
						if global_gene(end, 2)<blocks(j+1).start
							break;
						end
					if global_gene(1,1)>blocks(j+1).start
						gene_map1(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
					else
						gene_map1(1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
					end
					gidx = gidx-1;
				end
				gidx = 1;
				while gidx<length(blocks(j+1).prediction.genes);
					gene = blocks(j+1).prediction.genes{gidx};
					global_gene = gene(:, 1:2) + blocks(j+1).start-1;
					if global_gene(1, 1)>blocks(j).stop
						break;
					end
				if global_gene(end,2)<blocks(j).stop
					gene_map2(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				else
					gene_map2(global_gene(end, 2)-blocks(j+1).start+1:end) = gidx;
				end	
				gidx = gidx+1;
			end
			%figure;
			%plot(gene_map1>0), hold on , plot((gene_map2>0)+0.1, 'r');
	
			%find candidate split points 
			% (intergenic nucleotide close to the middle of the overlapping region)
			candidates = find((gene_map1+gene_map2)==0);
			idx1 = find(candidates<overlap, 1, 'last');
			idx2 = find(candidates>overlap, 1, 'first');

			candidate1 = candidates(idx1);
			candidate2 = candidates(idx2);

			%choose the candidate that is closes to the middle of the overlap region
			if abs(candidate1-overlap)<abs(candidate2-overlap)
				split_point = candidate1;
			else
				split_point = candidate2;
			end

			% find last gene left of the split
			map1_right_part = gene_map1(split_point:end);
			gene_idx = min(map1_right_part(map1_right_part>0))-1;
			if isempty(gene_idx)
				% no gene in this part of the overlap
				% => take all predicted genes from this block
				gene_idx = length(blocks(j).prediction.genes);
			end

			% remove predicted genes
			if gene_idx<length(blocks(j).prediction.genes)
				%assert that removed gene is in overlapping region
				gene = blocks(j).prediction.genes{gene_idx+1} + blocks(j).start-1;
				assert(all(gene(1,1)>blocks(j+1).start))
			end

			blocks(j).prediction.genes = blocks(j).prediction.genes(1:gene_idx);
			clear gene_idx 

			gene_idx = max(gene_map2(1:split_point))+1; % first gene to the right of the split

			if gene_idx>1
				%assert that removed gene is in overlapping region
				gene = blocks(j+1).prediction.genes{gene_idx-1}+blocks(j+1).start-1;
				assert(all(gene(end,2)<blocks(j).stop))
			end
			blocks(j+1).prediction.genes = blocks(j+1).prediction.genes(gene_idx:end);
			clear gene_idx

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% fill the maps again to see the plot
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			gene_map1 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
			gene_map2 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
			overlap = floor(length(gene_map1)/2);

			gidx = length(blocks(j).prediction.genes);
			while gidx>0
				gene = blocks(j).prediction.genes{gidx};
				global_gene = gene(:, 1:2) + blocks(j).start-1;
				if global_gene(end, 2)<blocks(j+1).start
					break;
				end
				if global_gene(1,1)>blocks(j+1).start
					gene_map1(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				else
					gene_map1(1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				end
				gidx = gidx-1;
			end
			gidx = 1;
			while gidx<length(blocks(j+1).prediction.genes);
				gene = blocks(j+1).prediction.genes{gidx};
				global_gene = gene(:, 1:2) + blocks(j+1).start-1;
				if global_gene(1, 1)>blocks(j).stop
					break;
			end
			if global_gene(end,2)<blocks(j).stop
				gene_map2(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
			else
				gene_map2(global_gene(end, 2)-blocks(j+1).start+1:end) = gidx;
			end
			gidx = gidx+1;
		end
		if rand>0.9
			assert(all(gene_map1.*gene_map2==0))
		end

		%figure;
		%plot((gene_map1>0)+0.02, 'b--'), hold on , plot((gene_map2>0)+0.12, 'r--');
		%keyboard
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		else %strand=='-'
			assert(blocks(j).stop>blocks(j+1).start)
			gene_map1 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
			gene_map2 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
			overlap = floor(length(gene_map1)/2);

			genes1 = blocks(j).prediction.genes(end:-1:1);
			genes2 = blocks(j+1).prediction.genes(end:-1:1);

			gidx = length(genes1);
			while gidx>0
				gene = genes1{gidx};
				global_gene = blocks(j).stop - gene(end:-1:1, 1:2)+1;
				if global_gene(end, 2)<blocks(j+1).start
					break;
				end
				if global_gene(1,1)>blocks(j+1).start
					gene_map1(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				else
					gene_map1(1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				end
				gidx = gidx-1;
			end
			gidx = 1;
			while gidx<length(genes2);
				gene = genes2{gidx};
				global_gene = blocks(j+1).stop - gene(end:-1:1, 1:2)+1;
				if global_gene(1, 1)>blocks(j).stop
					break;
				end
				if global_gene(end,2)<blocks(j).stop
					gene_map2(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				else
					gene_map2(global_gene(end, 2)-blocks(j+1).start+1:end) = gidx;
				end
				gidx = gidx+1;
			end
			%figure;
			%plot(gene_map1>0), hold on , plot((gene_map2>0)+0.1, 'r');


			candidates = find((gene_map1+gene_map2)==0);
			idx1 = find(candidates<overlap, 1, 'last');
			idx2 = find(candidates>overlap, 1, 'first');

			candidate1 = candidates(idx1);
			candidate2 = candidates(idx2);

			%choose the candidate that is closes to the middle of the overlap region
			if abs(candidate1-overlap)<abs(candidate2-overlap)
				split_point = candidate1;
			else
				split_point = candidate2;
			end

			% find last gene left of the split
			map1_right_part = gene_map1(split_point:end);
			gene_idx = min(map1_right_part(map1_right_part>0))-1;
			if isempty(gene_idx)
				% no gene in this part of the overlap
				% => take all predicted genes from this block
				gene_idx = length(genes1);
			end

			if gene_idx<length(genes1)
				%assert that removed gene is in overlapping region
				gene = genes1{gene_idx+1};
				gene = blocks(j).stop - gene(end:-1:1, 1:2)+1;
				assert(all(gene(1,1)>blocks(j+1).start))
			end
			blocks(j).prediction.genes = genes1(gene_idx:-1:1);
			clear gene_idx 

			gene_idx = max(gene_map2(1:split_point))+1; % first gene to the right of the split
			if gene_idx>1
				%assert that removed gene is in overlapping region
				gene = genes2{gene_idx-1};
				gene = blocks(j+1).stop - gene(end:-1:1, 1:2)+1;
				assert(all(gene(end,2)<blocks(j).stop))
			end
			blocks(j+1).prediction.genes = genes2(end:-1:gene_idx);
			clear gene_idx
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			gene_map1 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
			gene_map2 = zeros(1, blocks(j).stop-blocks(j+1).start+1);
			overlap = floor(length(gene_map1)/2);

			genes1 = blocks(j).prediction.genes(end:-1:1);
			genes2 = blocks(j+1).prediction.genes(end:-1:1);


			gidx = length(genes1);
			while gidx>0
				gene = genes1{gidx};
				global_gene = blocks(j).stop - gene(end:-1:1, 1:2)+1;
				if global_gene(end, 2)<blocks(j+1).start
					break;
				end
				if global_gene(1,1)>blocks(j+1).start
					gene_map1(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				else
					gene_map1(1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				end
				gidx = gidx-1;
			end
			gidx = 1;
			while gidx<length(genes2);
				gene = genes2{gidx};
				global_gene = blocks(j+1).stop - gene(end:-1:1, 1:2)+1;
				if global_gene(1, 1)>blocks(j).stop
					break;
				end
				if global_gene(end,2)<blocks(j).stop
					gene_map2(global_gene(1,1)-blocks(j+1).start+1:global_gene(end, 2)-blocks(j+1).start+1) = gidx;
				else
					gene_map2(global_gene(end, 2)-blocks(j+1).start+1:end) = gidx;
				end
				gidx = gidx+1;
			end
			if rand>0.9
				assert(all(gene_map1.*gene_map2==0))
			end
			%figure;
			%plot((gene_map1>0)+0.02, 'b--'), hold on , plot((gene_map2>0)+0.12, 'r--');
			%keyboard
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		end
		if isempty(blocks(j).prediction.genes)
			blocks(j).prediction.genes = {[]};
		end
		if isempty(blocks(j+1).prediction.genes)
			blocks(j+1).prediction.genes = {[]};
		end
	end
end		

genes = blocks2genes(blocks(logical(have_predictions)), PAR.model, 'prediction', genome_info, '');

%% filter out non coding genes
rm_idx = [];
for j = 1:length(genes)
	if length(genes(j).cds_exons)<1||isempty(genes(j).cds_exons{1})
		rm_idx = [rm_idx j];
	end
end
nc_genes = genes(rm_idx); 
genes(rm_idx) = [];
fprintf('found %i coding and  %i non coding genes. Stored in variable nc_genes in the same file\n', length(genes), length(nc_genes));

% rquant assumes half open intervals as input
% and outputs closed intervals
save(fn_predicted_genes_rquant, 'genes');

genes = convert_to_common_gene_structure(genes); 
save(fn_predicted_genes, 'genes', 'nc_genes');
save(fn_predicted_genes_date, 'genes', 'nc_genes');
end

