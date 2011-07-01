function write_gff3(organism, experiment, cleave, genes_prefix, genes)

base_dir = get_pred_dir(organism, experiment, cleave);

if nargin<4
	genes_prefix = 'genes'
end
fn_genes = sprintf('%s/%s.mat', base_dir, genes_prefix);
gff_fname = sprintf('%s/%s.gff3', base_dir, genes_prefix);

if nargin<5
	fprintf('load genes from file: %s\n', fn_genes);
	load(fn_genes, 'genes');
else
	fprintf('Do not load genes from file: %s\n', fn_genes);
end

genes = closed_to_half_open(genes);
for j = 1:length(genes)
	genes(j).utr5_exons = genes(j).utr_5prime;
	genes(j).utr3_exons = genes(j).utr_3prime;
	genes(j).is_valid = 1;
	genes(j).transcript_valid = ones(1, length(genes(j).exons));
	genes(j).coding = ~cellfun('isempty',genes(j).cds_exons);
end

addpath ~/svn/projects/rgasp.2/submission/
rgasp_write_gff3(genes,gff_fname, 'mGene')
