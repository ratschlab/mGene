function issimple=detectsimplegraph(gene)
%  issimple=detectsimplegraph(gene)
  
[dummy,exon_order] = sort(gene.splicegraph{1}(1,:));
gene.splicegraph{1} = gene.splicegraph{1}(:,exon_order);
gene.splicegraph{2} = gene.splicegraph{2}(exon_order,exon_order);

issimple=1 ;
for i=1:size(gene.splicegraph{1},2)-1
  if ~gene.splicegraph{2}(i,i+1)
    issimple=0; 
    return ;
  end ;
end ;
  
  