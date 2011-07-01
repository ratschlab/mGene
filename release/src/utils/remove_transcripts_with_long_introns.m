function genes = remove_transcripts_with_long_introns(genes, cutoff);
	
	rm_idx = [];
	for j = 1:length(genes)
		idx = [];
		for k=1:length(genes(j).exons)
			exons = genes(j).exons{k};
			intron_len = exons(2:end, 1)-exons(1:end-1, 2);
			if ~any(intron_len>cutoff)
				idx = [idx k]; 
			end
		end
		if isempty(idx)
			rm_idx = [rm_idx j];
			continue
		end	
		genes(j) = apply_transcript_idx(genes(j), idx);
	end	
	fprintf('removing %i genes with intron length>%i for all transcripts\n', length(rm_idx), cutoff)
	genes(rm_idx) = [];
return

