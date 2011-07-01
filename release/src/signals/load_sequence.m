function regions=load_sequence(regions);

% regions=load_sequence(regions);
  
genome_info = init_genome(regions(1).config);


for i=1:length(regions)
  if isfield(regions(i),'num_fragments')&&~isempty(regions(i).num_fragments)
    num_fragments = regions(i).num_fragments;
  else
    num_fragments = 1;
  end
  for p=1:num_fragments
    if mod(i,100)==0,
      fprintf(1,'region %i of %i\r',i,length(regions))
    end ;
    chr = regions(i).chr_num(p) ;
    str = regions(i).strand(p) ;
    start = regions(i).start(p) ;
    stop = regions(i).stop(p) ;
    seq = upper(load_genomic(genome_info.contig_names{chr},str,start,stop,genome_info));
    assert(ischar(seq)) ;
    idx = (~ismember(seq,'ACGT')) ;
    seq(idx)='N';
    if p==1
      regions(i).seq = seq;
    else
      regions(i).seq = [regions(i).seq ,seq];
    end
  end
end
