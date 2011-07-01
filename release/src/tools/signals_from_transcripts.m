function genes = signals_from_transcripts(genes)

for j = 1:length(genes)
	for k = 1:length(genes(j).exons)
		if isfield(genes, 'cds_exons') && length(genes(j).cds_exons)>=k && ~isempty(genes(j).cds_exons{k})
			if genes(j).strand=='+'
				genes(j).tis{k} = genes(j).cds_exons{k}(1, 1);
				genes(j).cdsStop{k} = genes(j).cds_exons{k}(end, 2)-3;
			else
				genes(j).tis{k} = genes(j).cds_exons{k}(end, 2);
				genes(j).cdsStop{k} = genes(j).cds_exons{k}(1, 1)+3;
			end
		end
		if isfield(genes, 'exons') && length(genes(j).exons)>=k && ~isempty(genes(j).exons{k})
			if genes(j).strand=='+'
				genes(j).tss{k} = genes(j).exons{k}(1, 1);
				genes(j).cleave{k} = genes(j).exons{k}(end, 2);
			else
				genes(j).tss{k} = genes(j).exons{k}(end, 2);
				genes(j).cleave{k} = genes(j).exons{k}(1, 1);
			end
		end
	end
end
