function [aa valid] = get_aa_seq(gene, k, genome_info)

	valid = 1;
	if isfield(gene, 'aa_seq') && iscell(gene.aa_seq) && length(gene.aa_seq)>=k
		aa = gene.aa_seq{k};
		valid = ~isempty(aa);
		return
	end
	if isempty(gene.cds_exons{k})
		aa = '';
		valid = 0;
		return
	end
	if gene.strand=='+'
		cds = load_genomic(gene.chr, char(gene.strand), gene.cds_exons{k}(:,1), gene.cds_exons{k}(:,2)-1, genome_info,1);
	else
		cds = load_genomic(gene.chr, char(gene.strand), gene.cds_exons{k}(:,1)+1, gene.cds_exons{k}(:,2), genome_info,1);
	end
	aa = translate(cds);
	if ~isequal(find(aa=='*'), length(aa))
		valid = 0;
	end
	if (aa(1)~='M')
		valid = 0;
	end
return

