function regions = get_chr_num(regions, genome_info)

if ~isstruct(genome_info)&&exist(genome_info, 'file')
	fn_config = genome_info;
	clear genome_info
	genome_info = init_genome(fn_config);
end

for j = 1:length(regions)
	regions(j).chr_num = strmatch(regions(j).chr,genome_info.contig_names,'exact') ;
end
