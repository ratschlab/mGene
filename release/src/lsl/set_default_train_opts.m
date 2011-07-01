function opts = set_default_train_opts()

% following variables are used to store functions and 
% parameters of additional input tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.track_functions = {};
opts.track_params = {};
opts.track_files = {};
opts.track_names = {};
opts.track_monoton_functions = {};

opts.segment_feature_functions = {};
opts.segment_feature_params = {};
opts.segment_feature_files = {};
opts.segment_feature_names = {};
opts.segment_feature_monoton_functions = {};

% default values for training options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.C=1 ;
opts.maxNumIter  = 200 ;
opts.maxNumBlock = inf ;
opts.block_design = 'merge' ;
opts.block_design_merge_num = 3 ;
opts.block_design_merge_sep = 1000 ;
opts.block_design_regionauto_num = 3 ;
opts.block_design_regionlist = '' ;
opts.long_trans_thresh = 1000;
opts.use_loss_mask = 1;
opts.scale_loss_mask = 1;
opts.cleave_offset = 0;
opts.use_train_region = 0;
opts.use_rna_seq_for_label_gen = 0;
opts.subsample_reads = 0;
opts.use_rna_seq_for_label_gen = 1;
