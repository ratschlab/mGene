function [blocks rm_idx] = create_blocks_from_genes(genes, maxoffset, mindist, rm_dist, fn_genome_config)
%maxoffset: maximal intergenic region created
%mindist:  if intergenic region between two genes is smaler than mindist take the whole region
%rm_dist: remove blocks with less intergenic region

reg_chr 	= {};
reg_strand 	= repmat('.', 1, length(genes));
reg_starts 	= zeros(1,length(genes));
reg_stops 	= zeros(1,length(genes));
reg_gene_ids 	= zeros(1,length(genes));
reg_complete 	= zeros(1,length(genes));
reg_status 	= zeros(1,length(genes));
left_off 	= zeros(1,length(genes));
right_off 	= zeros(1,length(genes));

genome_info = init_genome(fn_genome_config);
for j = 1:length(genome_info.flat_fnames), 
  d = dir(genome_info.flat_fnames{j});
  chr_len(j) = d.bytes;
end



reg_cnt = 1;
chr_nums = unique([genes.chr_num]);

all_gene_chr_idxs = store_gene_chr_idxs(genes, chr_nums);

for chr_num = chr_nums
  for s = '+-'
	fprintf('\rchr_num: %i (%i), strand:%s', chr_num, max(chr_nums), s);
    %chr_genes_idx = find([genes.chr_num]==chr_num&[genes.strand]==s);
	chr_genes_idx = all_gene_chr_idxs{chr_num}; % precomputed to make it faster for many contigs
	str_genes_idx = find([genes(chr_genes_idx).strand]==s);
	chr_genes_idx = chr_genes_idx(str_genes_idx);
    chr_genes_start = sort([genes(chr_genes_idx).start]);
    chr_genes_stop = sort([genes(chr_genes_idx).stop]);
  
    for j = chr_genes_idx

      % determin start of region
      last_stop = chr_genes_stop(find(chr_genes_stop<genes(j).start, 1, 'last'));
      start = max([last_stop, 1, genes(j).start-maxoffset]);
      if genes(j).start-start<mindist
		reg_starts(reg_cnt) = start;
      else
        reg_starts(reg_cnt) = start + floor(rand(1)*(genes(j).start-start-mindist));
      end
      left_off(reg_cnt) = genes(j).start-reg_starts(reg_cnt);

      assert(left_off(reg_cnt)<=maxoffset)

      % determin end of region
      first_start = chr_genes_start(find(chr_genes_start>genes(j).stop, 1, 'first'));
      stop = min([first_start, chr_len(chr_num), genes(j).stop+maxoffset]);
      if stop -genes(j).stop<mindist
		reg_stops(reg_cnt) = stop;
      else
		tmp = floor(rand(1)*(stop -genes(j).stop-mindist));
        reg_stops(reg_cnt) = stop - tmp;
      end
      right_off(reg_cnt) = reg_stops(reg_cnt)- genes(j).stop;

	  assert(reg_stops(reg_cnt)<=chr_len(chr_num))
      assert(right_off(reg_cnt)<=maxoffset)

      reg_chr{reg_cnt} = genes(j).chr;
      reg_strand(reg_cnt) = genes(j).strand;
      reg_gene_ids(reg_cnt) = genes(j).id;
      reg_complete(reg_cnt) = genes(j).is_complete;
      if ~isempty(genes(j).transcript_status)
      	reg_status(reg_cnt) = max(genes(j).transcript_status);
      end
      reg_cnt = reg_cnt+1;
    end
  end
end


blocks = struct( 'chr', reg_chr);

for j=1:length(blocks)
  blocks(j).start 	     = reg_starts(j);
  blocks(j).stop 	     = reg_stops(j);
  blocks(j).strand 	     = reg_strand(j);
  blocks(j).chr_num          = strmatch(upper(blocks(j).chr), upper(genome_info.contig_names), 'exact');
  assert(length(blocks(j).chr_num)==1) ;
  blocks(j).gene_ids 	     = reg_gene_ids(j);
  blocks(j).paralogs         = [];
  blocks(j).is_alt           = [];
  blocks(j).gene_complete    = reg_complete(j);
  blocks(j).gene_status      = reg_status(j);
  blocks(j).config           = [];
  blocks(j).split            = 1;
end

rm_idx = find(right_off<rm_dist | left_off<rm_dist);

assert(all(right_off<=maxoffset))
assert(all(left_off<=maxoffset))

return

x = histc(left_off, linspace(min(left_off), max(left_off),100));
plot(linspace(min(left_off), max(left_off),100), x)
mean(left_off)
min(left_off)
max(left_off)
figure
y = histc(right_off, linspace(min(right_off), max(right_off),100))
plot(linspace(min(right_off), max(right_off),100), y)
mean(right_off)
min(right_off)
max(right_off)

