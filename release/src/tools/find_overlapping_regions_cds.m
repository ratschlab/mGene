function idx = find_overlapping_regions_cds(genes1, genes2, consider_strands)
% idx = find_overlapping_regions_cds(genes1, genes2, consider_strands)

genes1 = set_start_stop(genes1);
genes2 = set_start_stop(genes2);

if nargin==3
	idx = find_overlapping_regions(genes1, genes2, consider_strands);
else
	idx = find_overlapping_regions(genes1, genes2);
end

return

function genes = set_start_stop(genes)

for j = 1:length(genes)
	min_start = inf;
	max_stop = 0;
	for k = 1:length(genes(j).exons)
        if (~isfield(genes(j), 'transcript_coding') || (isfield(genes(j), 'transcript_coding') && genes(j).transcript_coding(k))) && ~isempty(genes(j).cds_exons{k})
            if genes(j).cds_exons{k}(1,1)<min_start
                min_start = genes(j).cds_exons{k}(1,1);
            end
            if genes(j).cds_exons{k}(end,2)>max_stop
                max_stop = genes(j).cds_exons{k}(end,2);
            end
        else
            if ~isempty(genes(j).exons{k}),
                if genes(j).exons{k}(1,1)<min_start
                    min_start = genes(j).exons{k}(1,1);
                end
                if genes(j).exons{k}(end,2)>max_stop
                    max_stop = genes(j).exons{k}(end,2);
                end
            end ;
        end ;
	end
    if isinf(min_start), min_start=0 ; end ;

	genes(j).start = min_start;
	genes(j).stop = max_stop;
end		
return
