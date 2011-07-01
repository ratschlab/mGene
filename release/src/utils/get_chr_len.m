function len = get_chr_len(genome_info, chr_num)

d = dir(genome_info.flat_fnames{chr_num}) ;
if isempty(d)
  error('flat file with genomic sequence not found: %s', genome_info.flat_fnames{chr_num});
end
len = d.bytes;


