function genes = merge_genes_by_name(genes, fields)
% genes = merge_genes_by_name(genes, fields)

if nargin<2 || isequal(fields, []),
    fields=default_gene_field_names ;
end ;

gene_names = {genes.name} ;
[tmp,idx]=sort(gene_names) ;
genes=genes(idx) ;
rm_idx=[] ;
for i=1:length(genes)-1,
    j=1 ;
    while i+j<=length(genes) && isequal(genes(i).name, genes(i+j).name),
        genes(i)=merge_genes(genes(i), genes(i+j), fields) ;
        rm_idx(end+1) = i+j ;
        j=j+1 ;
    end ;
end ;
genes(rm_idx) = [] ;
assert(length(unique(gene_names))==length(genes)) ;

for i=1:length(genes),
    genes(i) = update_start_stop(genes(i));
end ;

return ;

function g1 = merge_genes(g1, g2, fields)
% g1 = merge_genes(g1, g2, fields)

for i=1:length(g1.exons),
    rm_idx=[] ;
    for j=1:length(g2.exons),
        id=1 ;
        if isfield(g1, fields{i}),
            if ~isequal(g1.(fields{i})(i), g2.(fields{i})(j)),
                id=0 ;
            end ;
        end
        if id,
            rm_idx(end+1) = j ;
        end ;
    end ;
    g2 = apply_transcript_idx(g2, setdiff(1:length(g2.exons), rm_idx), fields, 1) ;
end ;

for i=1:length(fields)
	if isfield(g1, fields{i}) 
		g1.(fields{i}) 	= [g1.(fields{i}) g2.(fields{i})] ;
    end
end 
