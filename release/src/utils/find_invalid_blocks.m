function rm_idx = find_invalid_blocks(blocks, genome_info)

rm_idx = [];

for j = 1:length(blocks)
	if blocks(j).start<1
		rm_idx = [rm_idx j];
	elseif blocks(j).stop>get_chr_len(genome_info, blocks(j).chr_num)
		rm_idx = [rm_idx j];
	end
end
