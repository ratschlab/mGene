function genes = remove_nc_transcripts(genes);
	
	rm_idx = [];
	for j = 1:length(genes)
		idx = find(~cellfun('isempty' , genes(j).cds_exons));
		if ~isempty(idx)
			genes(j) = apply_transcript_idx(genes(j), idx);
		else
			rm_idx = [rm_idx j];
		end
	end	
	if length(rm_idx)>0
		fprintf('removing %i(%i) nc genes\n', length(rm_idx), length(genes));
		genes(rm_idx) = [];
	end
return

