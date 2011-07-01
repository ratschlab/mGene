function [blocks bad_idx]= assign_genes_to_blocks(blocks, genes)


block_chrs = unique({blocks.chr})
gene_chrs = unique({genes.chr})
chrs = intersect(block_chrs, gene_chrs);

if isempty(chrs)
	error('gene chromosome names do not match block chromosome names\n');
end

bad_idx = [];
for j=1:length(chrs)
	for s = '+-'
		block_ids = find(strcmp({blocks.chr}, chrs{j}) & [blocks.strand]==s);
		gene_ids = find(strcmp({genes.chr}, chrs{j}) & [genes.strand]==s);

		for b = 1:length(block_ids)
			[in_region, ok] = find_overlapping_genes(blocks(block_ids(b)), genes(gene_ids));
			if ok
				if length(in_region)>1
					fprintf('%i genes in region found\n', length(in_region));
					%keyboard
				end
				blocks(block_ids(b)).gene_ids = [genes(gene_ids(in_region)).id];
				if isfield(genes, 'is_complete')
					blocks(block_ids(b)).gene_complete = [genes(gene_ids(in_region)).is_complete];
				else
					blocks(block_ids(b)).gene_complete = 1;
				end
				if isfield(genes, 'transcript_status')
					blocks(block_ids(b)).gene_status = max([genes(gene_ids(in_region)).transcript_status]);
				else
					blocks(block_ids(b)).gene_status = 1;
				end
			else
				bad_idx = [bad_idx block_ids(b)];
			end
		end
	end
end


return

function [gene_ids, ok] = find_overlapping_genes(block, genes)
% all blocks and genes are on the same chr and strand

	start_in = find([genes.start]>=block.start & [genes.start]<=block.stop);
	stop_in	= find([genes.stop ]>=block.start & [genes.stop ]<=block.stop);

	complete_in = intersect(start_in, stop_in);

	if length(start_in)>length(complete_in)||length(stop_in)>length(complete_in)
		% discard block
		ok = 0;
		gene_ids = [];
		return
	end
	ok = 1; 
	gene_ids = complete_in;
return
