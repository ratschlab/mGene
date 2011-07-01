function len = trans_len(exons, close_interval_convention, transcripts)

if nargin<2
	close_interval_convention = 0;
end

if isstruct(exons)
	genes = exons;
	len = zeros(1, length(genes)*10);
	cnt = 0;
	for j = 1:length(genes)
		for k = 1:length(genes(j).exons)
			cnt = cnt+1;
			len(cnt) = trans_len(genes(j).exons, close_interval_convention, k);
		end
	end
	len(cnt+1:end) = []; 
	return
end

	if nargin<3
		len = zeros(1, length(exons));
	else
		len = zeros(1, length(transcripts));
	end
	if nargin<3
		transcripts = 1:length(exons);
	end
	cnt = 1;
	for k = transcripts
		if ~isempty(exons{k})
			if close_interval_convention
				len(cnt) = sum(exons{k}(:, 2)-exons{k}(:, 1)+1);
			else
				len(cnt) = sum(exons{k}(:, 2)-exons{k}(:, 1));
			end
		end
		cnt = cnt+1;
	end

return

