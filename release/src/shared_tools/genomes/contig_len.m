function len = contig_len(contig, genome_info)
% len = contig_len(contig, genome_info)

if ischar(contig),
  %contig_idx = find(ismember(upper(genome_info.contig_names),upper(contig))) ;
  contig_idx = strmatch(upper(contig), upper(genome_info.contig_names), 'exact') ;
  assert(all(size(contig_idx)==1) & ~isempty(contig_idx)) ;
else
  contig_idx = contig ;
end ;
  
fname=genome_info.flat_fnames{contig_idx} ;

d=dir(fname) ; 
len=d.bytes;
