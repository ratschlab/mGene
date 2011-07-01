function num = trans_cnt(genes)

num = zeros(1, length(genes));
for j = 1:length(genes)
	num(j) = length(genes(j).exons);
end
