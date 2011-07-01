function blocks = scale_loss_according_to_confirmation(blocks)


if ~isfield(blocks, 'gene_status')
	return
end
if ~isfield(blocks, 'loss_mask')
	return
end

for j = 1:length(blocks)
	if max(blocks(j).gene_status)==1
		% assume genes is cDNA confirmed
		blocks(j).loss_mask = blocks(j).loss_mask.*3;
	elseif max(blocks(j).gene_status)==0
		% assume genes is EST confirmed
		blocks(j).loss_mask = blocks(j).loss_mask.*1.5;
	elseif max(blocks(j).gene_status)==-1
		% assume gene is not confirmed
		% do nothing
	else
		warning('gene_status has some unexpected value')
	end
end
