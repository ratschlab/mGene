function idx = find_overlapping_regions(regions1, regions2, consider_strands)
%%

if isempty(regions1) & isempty(regions2)
  idx = [];
  return;
end
  
if isempty(regions2)&&(nargin==2 || consider_strands)
	%% find overlapping pairs within one set of regions
	chr_nums = unique([regions1.chr_num]);
	all_gene_chr_idxs1 = store_gene_chr_idxs(regions1, chr_nums);%fast if there are many contigs
	idx = [];
	for chr=chr_nums
		for s='+-'
			fprintf('\rchr%i strand%s',chr,s)
			chr_genes_idx1 = all_gene_chr_idxs1{chr}; 
			str_genes_idx1 = find([regions1(chr_genes_idx1).strand]==s);
			r1_chr_idx = chr_genes_idx1(str_genes_idx1);

			if isempty(r1_chr_idx)
				continue
			end 
			len = max([[regions1(r1_chr_idx).stop]]);
			map1 = zeros(1, len, 'uint16');
			for r = r1_chr_idx 
				area = map1(regions1(r).start:regions1(r).stop);

				if any(area>0)
					o_regs = unique(area);
					if o_regs(1)==0
						o_regs = o_regs(2:end);
					end
					%append overlapping regions
					idx = [idx; [ones(length(o_regs), 1)*r, o_regs']];
				end
				map1(regions1(r).start:regions1(r).stop) = r;
			end
		end
	end
	return
end

chr_nums = unique([[regions1.chr_num] [regions2.chr_num]]);
all_gene_chr_idxs1 = store_gene_chr_idxs(regions1, chr_nums);%fast if there are many contigs
all_gene_chr_idxs2 = store_gene_chr_idxs(regions2, chr_nums);

if nargin==2 || consider_strands
	idx = [];
	for chr=chr_nums
		for s='+-'
			fprintf('\rchr%i strand%s',chr,s)
			%r1_chr_idx = find([regions1.chr_num]==chr&[regions1.strand]==s);
			%r2_chr_idx = find([regions2.chr_num]==chr&[regions2.strand]==s);
			chr_genes_idx1 = all_gene_chr_idxs1{chr}; % precomputed to make it faster for many contigs
			str_genes_idx1 = find([regions1(chr_genes_idx1).strand]==s);
			r1_chr_idx = chr_genes_idx1(str_genes_idx1);

			chr_genes_idx2 = all_gene_chr_idxs2{chr}; 
			str_genes_idx2 = find([regions2(chr_genes_idx2).strand]==s);
			r2_chr_idx = chr_genes_idx2(str_genes_idx2);

			if isempty(r1_chr_idx)||isempty(r2_chr_idx)
				continue
			end 
			len = max([[regions1(r1_chr_idx).stop] [regions2(r2_chr_idx).stop]]);
			map1 = zeros(1, len, 'uint16');
			map2 = zeros(1, len, 'uint16');
			for r = r1_chr_idx 
				map1(regions1(r).start:regions1(r).stop) = r;
			end
			for r = r2_chr_idx 
				map2(regions2(r).start:regions2(r).stop) = r;
			end
			overlap_map = uint16((map1.*map2)>0);
			x = unique([[map1.*overlap_map]' [map2.*overlap_map]'], 'rows');
			assert(x(1,1)==0&&x(1,2)==0)	
			idx = [idx;x(2:end, :)];
		end
	end
else
	idx = [];
	for chr=chr_nums
		fprintf('\rchr%i',chr)
		r1_chr_idx = all_gene_chr_idxs1{chr};
		r2_chr_idx = all_gene_chr_idxs2{chr};

		if isempty(r1_chr_idx)||isempty(r2_chr_idx)
			continue
		end 
		len = max([[regions1(r1_chr_idx).stop] [regions2(r2_chr_idx).stop]]);
		map1 = zeros(1, len, 'uint16');
		map2 = zeros(1, len, 'uint16');
		for r = r1_chr_idx 
			map1(regions1(r).start:regions1(r).stop) = r;
		end
		for r = r2_chr_idx 
			map2(regions2(r).start:regions2(r).stop) = r;
		end
		overlap_map = uint16((map1.*map2)>0);
		x = unique([[map1.*overlap_map]' [map2.*overlap_map]'], 'rows');
		assert(x(1,1)==0&&x(1,2)==0)	
		idx = [idx;x(2:end, :)];
	end
end


return

