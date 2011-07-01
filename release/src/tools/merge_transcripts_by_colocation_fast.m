function genes=merge_transcripts_by_colocation_fast(genes, fields, merge_by_cds, verb)
%genes=merge_transcripts_by_colocation_fast(genes, fields, merge_by_cds, verbosity)

if nargin==1||isempty(fields)
	fields = {'transcripts', 'exons', 'cds_exons', 'utr3_exons', 'utr5_exons', 'tss', 'tis', 'cdsStop', 'cleave'} ;
end
if nargin<3
	merge_by_cds = 0;
end

idx = 1;
while ~isempty(idx)
	nofgenes = length(genes);
	[genes] = find_start_stop(genes, merge_by_cds);
	[genes.remove] = deal(0);
	consider_strands = 1;
	idx = find_overlapping_regions(genes, [], consider_strands);
	assert(issorted(idx(:, 1)))
	for row = 1:size(idx, 1)
		if genes(idx(row, 1)).remove==1 || genes(idx(row, 2)).remove==1
			% avoid duplicating transcripts
			continue
		end
		genes(end+1) = merge_genes(genes(idx(row, 1)), genes(idx(row,2)), fields);
		genes(idx(row, 1)).remove=1;
		genes(idx(row, 2)).remove=1;
	end
	rm_idx = find([genes.remove]);
	if nargin<4||verb
		fprintf('reducing number of genes from %i to %i\n', nofgenes, nofgenes-length(rm_idx)/2);
	end
	genes(rm_idx) = [];
	if isempty(rm_idx), 
		break; 
	end
end
genes = find_start_stop(genes, 0);
genes = rmfield(genes, 'remove');

return

function genes = find_start_stop(genes, merge_by_cds)

	for i=1:length(genes) ;
		start=inf; stop=0 ;
		for j=1:length(genes(i).cds_exons)
			if ~isempty(genes(i).cds_exons{j}) && merge_by_cds,
				start=min(start, min(genes(i).cds_exons{j}(:)));
				stop=max(stop, max(genes(i).cds_exons{j}(:)));
			else
				start=min(start, min(genes(i).exons{j}(:)));
				stop=max(stop, max(genes(i).exons{j}(:)));
			end ;
		end ;
		assert(~isinf(start) & stop~=0) ;
		genes(i).start=start ;
		genes(i).stop=stop ;
	end
return

function gene1 = merge_genes(gene1, gene2, fields)
	assert(gene1.strand==gene2.strand) ;
	for f= 1:length(fields),
		if isfield(gene1, fields{f}),
			if iscell(gene1.(fields{f})) && iscell(gene2.(fields{f})),
				gene1.(fields{f}) = {gene1.(fields{f}){:} gene2.(fields{f}){:}} ;
			elseif iscell(gene1.(fields{f}))
				gene1.(fields{f}) = {gene1.(fields{f}){:} {}} ;
			elseif iscell(gene1.(fields{f}))
				gene1.(fields{f}) = {{} gene2.(fields{f}){:}} ;
			else
				gene1.(fields{f}) = [gene1.(fields{f}) gene2.(fields{f})] ;
			end
		end
	end 
return

