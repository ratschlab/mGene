function genes = compute_max_weight(genes)

	for j = 1:length(genes)
		max_weight = max(genes(j).transcript_weights);
		genes(j).max_transcript_weight = max_weight;
	
		% set low random weight if there are no quantification values given
		no_weights=isempty(max_weight);
		no_weights=no_weights||max_weight==0;% for asigning all genes to specific bins
		no_weights=no_weights||max_weight==-1;
		no_weights=no_weights||isnan(max_weight);
		if no_weights
			genes(j).max_transcript_weight = rand(1)/1e18;
		end
	end
return


