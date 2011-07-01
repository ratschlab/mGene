function create_region_list(fn_blocks, outfile, fn_val_blocks)

%blocks = load_struct('/fml/ag-raetsch/share/projects/genefinding/A_thaliana/jonas_new_genes/exp/lsl/tiling_C1/data/blocks_all.mat', 'blocks');
if nargin==3
	fn_train_blocks = fn_blocks;
	fprintf('writing training blocks\n')
	blocks = load_struct(fn_train_blocks, 'blocks');
	write_region_list_as_gff(blocks, outfile, 0, 'training')
	clear blocks

	fprintf('appending validation blocks\n')
	blocks = load_struct(fn_val_blocks, 'blocks');
	write_region_list_as_gff(blocks, outfile, 1, 'validation')

else
	blocks = load_struct(fn_blocks, 'blocks');
	write_region_list_as_gff(blocks, outfile)
end


return


