function bins = myhistc(x, edges)

x = sort(x);

bins = zeros(1, length(edges));
k = 1;
for j = 1:length(edges)-1
	while k<=length(x) && x(k)<=edges(j+1)
		if x(k)>edges(j)
			bins(j) = bins(j)+1;
		end
		k = k+1;
	end
end
