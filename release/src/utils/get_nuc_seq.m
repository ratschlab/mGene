function [nuc valid] = get_nuc_seq(gene, k, genome_info)

	valid = 1;
	if isfield(gene, 'gene2.nuc_seq') && iscell(gene.nuc_seq)
		assert(length(gene.nuc_seq)==length(gene.exons))
		nuc = gene.nuc_seq{k};
		valid = ~isempty(nuc);
		return
	end
	if isempty(gene.exons{k})
		nuc = '';
		valid = 0;
		return
	end
	if gene.strand=='+'
		nuc = load_genomic(gene.chr, char(gene.strand), gene.exons{k}(:,1), gene.exons{k}(:,2)-1, genome_info,1);
	else
		nuc = load_genomic(gene.chr, char(gene.strand), gene.exons{k}(:,1)+1, gene.exons{k}(:,2), genome_info,1);
	end
return

