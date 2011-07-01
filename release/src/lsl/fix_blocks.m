function blocks = fix_blocks(blocks, model) ;
% blocks = fix_blocks(blocks, model) ;

assert(isequal(model, fix_model(model))) ;
w_fields = fieldnames(init_weights(model)) ;

for i=1:length(blocks),
	if isfield(blocks, 'truth')
		for j=1:length(blocks(i).truth),
			blocks(i).truth(j).weights = reorder_fields(blocks(i).truth(j).weights, w_fields) ;
        	end ;
	end
	if isfield(blocks, 'pred')
		for j=1:length(blocks(i).pred),
			blocks(i).pred(j).weights = reorder_fields(blocks(i).pred(j).weights, w_fields) ;
        	end ;
	end
end ;
