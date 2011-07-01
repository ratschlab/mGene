function genes_out = convert_to_common_gene_structure(genes)

genes_out = struct;
genes_out.name = '';
genes_out.chr = '';
genes_out.chr_num = [];
genes_out.strand = '';
genes_out.start = [];
genes_out.stop 	= [];
genes_out.transcripts = {};
genes_out.exons = {};
genes_out.cds_exons = {};
genes_out.utr_5prime = {};
genes_out.utr_3prime = {};


for j = 1:length(genes)
	genes_out(j).name = ['mGene_' genes(j).name]; 
	genes_out(j).start = genes(j).start; 
	genes_out(j).stop = genes(j).stop; 
	genes_out(j).chr = genes(j).chr;
	genes_out(j).chr_num = genes(j).chr_num;
	genes_out(j).strand = genes(j).strand;
	for k = 1:length(genes(j).transcripts)
		genes_out(j).transcripts{k} = ['mGene_' genes(j).transcripts{k}];
		exons = genes(j).exons{k};
		cds_exons = genes(j).cds_exons{k};
		utr5_exons = genes(j).utr5_exons{k};
		utr3_exons = genes(j).utr3_exons{k};
		if genes(j).strand=='+'
			genes_out(j).exons{k} = [exons(:,1), exons(:,2)-1];
			if ~isempty(cds_exons)
				genes_out(j).cds_exons{k} = [cds_exons(:,1), cds_exons(:,2)-1];
			else
				genes_out(j).cds_exons{k} = [];
			end
			if ~isempty(utr5_exons)
				genes_out(j).utr_5prime{k} = [utr5_exons(:,1), utr5_exons(:,2)-1];
			else
				genes_out(j).utr_5prime{k} = [];
			end
			if ~isempty(utr3_exons)
				genes_out(j).utr_3prime{k} = [utr3_exons(:,1), utr3_exons(:,2)-1];
			else
				genes_out(j).utr_3prime{k} = [];
			end
		else
			genes_out(j).exons{k} = [exons(:,1)+1, exons(:,2)];
			if ~isempty(cds_exons)
				genes_out(j).cds_exons{k} = [cds_exons(:,1)+1, cds_exons(:,2)];
			else
				genes_out(j).cds_exons{k} = [];
			end
			if ~isempty(utr5_exons)
				genes_out(j).utr_5prime{k} = [utr5_exons(:,1)+1, utr5_exons(:,2)];
			else
				genes_out(j).utr_5prime{k} = [];
			end
			if ~isempty(utr3_exons)
				genes_out(j).utr_3prime{k} = [utr3_exons(:,1)+1, utr3_exons(:,2)];
			else
				genes_out(j).utr_3prime{k} = []; 
			end
		end
	end
end
