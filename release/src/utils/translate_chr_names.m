function genes = translate_chr_names(genes, names1, names2)

% e.g.: 
% names1 = unique({l.genes.chr})
%
% which makes names1 to be 
% names1 = {'Chr1' 'Chr2' 'Chr3' 'Chr4' 'Chr5' 'ChrC' 'ChrM'}
% names2 = {'1'  '2'  '3'  '4'  '5' 'chloroplast' 'mitochondria'}

% assume a one to one relationship
assert(length(names1)==length(names2))

for j = 1:length(genes)
	matchidx = strmatch(genes(j).chr, names1);	

	if isempty(matchidx)
		error('chr name not found in cell');
	end	
	if length(matchidx)>1
		error('names1 is not unique')
	end

	genes(j).chr = names2{matchidx};
end

