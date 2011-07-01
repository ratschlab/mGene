function blocks = generate_loss_mask(blocks, model);

%warning('generate_loss_mask not yet active')
%return


fields = fieldnames(model.segments);

for j = 1:length(blocks)
	loss_mask = ones(1, blocks(j).stop-blocks(j).start+1);
	for f = 1:length(fields)
		seg_type = fields{f};
		exon_map = zeros(1, blocks(j).stop-blocks(j).start+1);
		transcnt = 0;
		for k=1:length(blocks(j).truth)
			segments = blocks(j).truth(k).segments;
			if isempty(segments), continue ; end ;
			transcnt = transcnt+1;
			for l = 1:size(segments, 1)
				if ismember(segments(l, 3),  model.segments.(seg_type))
					exon_map(segments(l, 1):segments(l, 2)-1) = exon_map(segments(l, 1):segments(l, 2)-1)+1;
				end
			end
		end
		loss_mask = loss_mask.*(exon_map==0|exon_map==transcnt);
	end
	% for nucleotides that are exonic in all transcripts the 
	% value in the map should be equal to the number of transcripts;
	% for unambigious intergenic or intronic nucleotides the 
	% value should be equal to zero
	blocks(j).loss_mask = loss_mask;

%	rand('seed', 12345)
%	if transcnt>1&& rand>0.8
%		figure; hold on
%		for k=1:length(blocks(j).truth)
%			plot_segments(segments, model, k)
%		end
%		plot(blocks(j).loss_mask)
%		plot(exon_map, 'r')
%    end
end ;

