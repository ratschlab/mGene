function blocks = make_prediction_blocks(PAR, blocks, model, generate_features)
% blocks = make_prediction_blocks(PAR, blocks, model, generate_features)

if nargin<4, 
  generate_features=1 ;
end ;
% load sequence
%-------------------------------------------------------------
if ~isfield(blocks, 'seq')
  blocks = load_sequence(blocks);
end
% add signals to blocks
%-------------------------------------------------------------
fprintf(1,'load signals and append to blocks \n')

signal_names = fieldnames(PAR.model.signals) ;
blocks = add_signals2blocks(blocks,PAR, signal_names, 'pred_dir', 0);

if 0
	fprintf('manipulating signal predictions with masks generated from confirmed sequences\n')
	signal_masks = load_struct('/fml/ag-raetsch/share/databases/genomes/C_elegans/elegans_WS200/genebuild//signal_masks_withends.mat', 'signal_masks');
	%signal_masks = load_struct('/fml/ag-raetsch/share/databases/genomes/C_elegans/elegans_WS200/genebuild//signal_masks_withends_with_prot.mat', 'signal_masks');

	addpath('~/svn/projects/genefinding/data_processing_lsl/signal_masks')

	if blocks(1).chr==1
		save(sprintf('~/tmp/blocks_before%s', blocks(1).strand), 'blocks')
	end
	for j = 1:length(blocks)
		blocks(j) = apply_signal_mask_new_blocks(blocks(j), signal_masks);
	end
	if blocks(1).chr==1
		save(sprintf('~/tmp/blocks_after%s', blocks(1).strand), 'blocks')
	end
	rmpath('~/svn/projects/genefinding/data_processing_lsl/signal_masks')
end

% generate feature matrix states x pos x features 
%-------------------------------------------------------------
if generate_features,
  fprintf(1,'generate feature matrix for blocks \n')
  blocks = gen_features(blocks, model);
end ;

