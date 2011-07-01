function idx = find_genes_with_nonoverlaping_transcripts(genes)

idx = [];
for j = 1:length(genes)
	exon_map = zeros(1, genes(j).stop-genes(j).start+1);
	minstart = inf;
	maxstop = -inf;
	for k=1:length(genes(j).exons)
		start = min(min(genes(j).exons{k}))-genes(j).start+1;
		minstart = min(minstart, start);
		stop = max(max(genes(j).exons{k}))-genes(j).start+1;
		maxstop = max(maxstop, stop);
		exon_map(start:stop) = 1;
	end
	if mean(exon_map(minstart:maxstop))<1
		idx = [idx j];
	end
end
