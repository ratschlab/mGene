function genes = sort_transcripts(genes)
	
	
	for j = 1:length(genes)
		closed_interval_convention = 0;
		if isfield(genes, 'cds_exons')
			len = trans_len(genes(j).cds_exons, closed_interval_convention);
		else
			len = trans_len(genes(j).exons, closed_interval_convention);
		end
		[tmp idx]  = sort(len, 'descend');
		genes(j) = apply_transcript_idx(genes(j), idx);
	end	
return

