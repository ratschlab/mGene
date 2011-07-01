function write_coding_sequence_to_fasta(genes, genome_info, fasta_file, translate_to_AA)

if nargin<4
	translate_to_AA=0;
end
fd = fopen(fasta_file, 'w');
for j = 1:length(genes)
	exons = genes(j).cds_exons{1}; 
	[cds, splice0, ok] = load_genomic(genes(j).chr_num, char(genes(j).strand), exons(:,1), exons(:,2), genome_info, 1);
	fprintf(fd, '>%s\n', genes(j).name);
	if translate_to_AA
		fprintf(fd, '%s\n', translate(cds));
	else
		fprintf(fd, '%s\n', cds);
	end
end
fclose(fd);
