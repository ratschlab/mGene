function regions=init_regions(fn_genome_config, organism)
% regions=init_regions(fn_genome_config,organism)

  
regions(1).chr = [];
regions(1).chr_num = [];
regions(1).strand = [];
regions(1).start = [];
regions(1).stop = [];
regions(1).offset =[];
regions(1).id = [];
regions(1).organism = [];
regions(1).config = [];

if ~exist('fn_genome_config')
  return ;
else
  regions = regions([]) ;
  genome_config = init_genome(fn_genome_config);
  strands = '+-';
  for c=1:length(genome_config.contig_names)
    %+strand
    d = dir(genome_config.flat_fnames{c}) ;
    for s=1:length(strands)
      region = [];
      region.chr=genome_config.contig_names{c};
      region.chr_num = c;
      region.strand=strands(s);
      region.start = 1;
      region.stop = d.bytes ;
      switch strands(s)
       case '+',
        region.offset      = region.start - 1;
       case '-',
        region.offset      = region.stop + 1;
      end
      region.id = (c-1)*2+s;
      region.organism = [];
      if exist('organism')
        region.organism = organism;
      end
      region.config = fn_genome_config;
      regions(region.id) = region;
    end
  end
end

