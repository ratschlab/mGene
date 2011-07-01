function gene = update_start_stop(gene)

min_start = inf;
max_stop = 0;
for k = 1:length(gene.exons)
    if ~isempty(gene.exons{k}),
        if gene.exons{k}(1,1)<min_start
            min_start = gene.exons{k}(1,1);
        end
        if gene.exons{k}(end,2)>max_stop
            max_stop = gene.exons{k}(end,2);
        end
    end ;
end
if isinf(min_start), min_start=0 ; end ;
assert(max_stop>=min_start) ;

gene.start = min_start;
gene.stop = max_stop;
