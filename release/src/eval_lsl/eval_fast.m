function [correct SN SP F] = eval_fast(genes, agenes)

if ischar(agenes)
	organism = agenes;
	fn_genes = rgasp_label_dir(organism, 'eval')
	l = load(fn_genes, 'genes');
	agenes = l.genes;
end

genes_bk = genes;

for j = 1:3 % check if half open or closed interval

	if j == 2
		genes = closed_to_half_open(genes_bk);
	elseif j == 3
		genes = half_open_to_closed(genes_bk);
	end

	% collect exons from annotated genes
	[acoding, aexons] = get_features(agenes);
	
	% collect exons from predicted genes
	[coding, exons] = get_features(genes);
	
	num_wrong = size(setdiff(exons, aexons, 'rows'), 1);
	num_correct = size(exons, 1)-num_wrong;
	SN = num_correct/size(aexons,1);
	SP = num_correct/size(exons, 1);
	F = fscore(SN, SP);
	if SP>0.1
		fprintf('coding exon level: SN: %.4f, SP: %.4f, F: %.4f\n', SN, SP, F)
		break
	end
end

% eval transcripts
correct = zeros(1, length(genes));
idx = find_overlapping_regions(genes, agenes);
for j = 1:size(idx, 1)
	correct(idx(j, 1)) = correct(idx(j, 1)) + compare(genes(idx(j, 1)), agenes(idx(j, 2)));
end
SN = sum(correct)/sum(acoding);
SP = mean(correct(coding>0));
F = fscore(SN, SP);
fprintf('coding transcript level: SN: %.4f, SP: %.4f, F: %.4f\n', SN, SP, F)

return

function [coding, exons] = get_features(genes)
	coding = zeros(1,length(genes));
	exons = zeros(50*length(genes), 4); % chr strand start stop
	cnt = 0;
	for j = 1:length(genes)
		coding(j) = sum(~cellfun('isempty', genes(j).cds_exons));
		tmp = collect_exons(genes(j));
		exons(cnt+1:cnt+size(tmp, 1), :) = tmp;
		cnt = cnt+size(tmp, 1);
	end
	exons(cnt+1:end, :) = [];
	exons = unique(exons, 'rows');
return

function exons = collect_exons(gene)

	exons = [];
	for j = 1:length(gene.cds_exons)
		if isempty(gene.cds_exons{j})
			continue
		end
		num_exons = size(gene.cds_exons{j}, 1);
		chr = repmat(gene.chr_num, num_exons, 1);	
		strand = repmat(double(gene.strand), num_exons, 1);	
		exons = [exons; [chr strand gene.cds_exons{j}(:, 1:2)]];
	end
	exons = unique(exons, 'rows');

return

function res = compare(gene1, gene2)

	res = 0;
	for j = 1:length(gene1.cds_exons)
		for k = 1:length(gene2.cds_exons)
			if isempty(gene1.cds_exons{j}) || isempty(gene2.cds_exons{k})
				continue
			end
			if isequal(gene1.cds_exons{j}(:, 1:2), gene2.cds_exons{k}(:, 1:2))
				res = res+1;
			end
		end
	end

return
