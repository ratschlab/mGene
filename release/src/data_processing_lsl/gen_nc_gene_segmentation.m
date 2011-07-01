function gene = gen_nc_gene_segmentation(gene, segment_ids)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('gen_nc_gene_segmentation\n')
for j = 1:length(gene.exons)
	isoform = [];
		
	if gene.strand=='+'
		for k=1:size(gene.exons{j},1)-1
			isoform = [isoform; [gene.exons{j}(k,1:2) segment_ids.nc_exon]];
			isoform = [isoform; [gene.exons{j}(k,2) gene.exons{j}(k+1,1) segment_ids.intron]];
		end
		isoform = [isoform; [gene.exons{j}(end,1:2) segment_ids.nc_exon]];
	else
		isoform = [isoform; [gene.exons{j}(1,1:2) segment_ids.nc_exon]];
		for k=2:size(gene.exons{j},1)
			isoform = [isoform; [gene.exons{j}(k-1,2) gene.exons{j}(k,1) segment_ids.intron]];
			isoform = [isoform; [gene.exons{j}(k,1:2) segment_ids.nc_exon]];
		end
	end

	gene.segments{j} = isoform;
	assert(isequal(isoform(2:end, 1), isoform(1:end-1, 2)));
end
%fprintf('isoform size: %i\n', size(isoform, 1)) 
