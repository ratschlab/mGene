function print_train_opts(opts, fd_out)

for fd=unique([1 fd_out]),
	fprintf(fd, '\nThe following options are used:\n') ;
	fprintf(fd, '\tMaximal number of iterations:\t\t %i\n', opts.maxNumIter) ;
	fprintf(fd, '\tThreshold for long transitions:\t\t %i\n', opts.long_trans_thresh) ;
	fprintf(fd, '\tRegularization parameter C:\t\t %1.2f\n', opts.C) ;
	fprintf(fd, '\tCleave offset to tune SN SP tradeoff:\t %i\n', opts.cleave_offset) ;
	fprintf(fd, '\tMaximal number of training blocks:\t %1.2f\n', opts.maxNumBlock) ;
	fprintf(fd, '\tSuppress loss in alternative regions:\t %i\n', opts.use_loss_mask) ;
	if opts.use_rna_seq_for_label_gen
		fprintf(fd, '\tUse rna_seq data for label generation:\t track no. %i\n', opts.use_rna_seq_for_label_gen) ;
	end
	fprintf(fd, '\tSubsample reads:\t\t\t %i\n', opts.subsample_reads) ;
	fprintf(fd, '\tscale loss according to confirmation:\t %i\n', opts.scale_loss_mask) ;
	if opts.use_train_region
		fprintf(fd, '\ttake training examples only from regions:\n') ;
		fprintf(fd, '\t\t%s\n', opts.train_region_file) ;
	end
	fprintf(fd, '\tBlock design\t\t\t\t %s\n', opts.block_design) ;
	if isequal(opts.block_design, 'merge')
		fprintf(fd, '\t\tMerge %i genes into one block\n', opts.block_design_merge_num) ;
		fprintf(fd, '\t\tSeparation between genes\t %ibp\n', opts.block_design_merge_sep) ;
	elseif isequal(opts.block_design, 'regionauto')
		fprintf(fd, '\t\tMerge %i genes into one block\n', opts.block_design_regionauto_num) ;
	elseif isequal(opts.block_design, 'regionlist')
		fprintf(fd, '\t\tMerge %i genes into one block\n', opts.block_design_merge_num) ;
		fprintf(fd, '\t\tRegion list:\t\t\t%s\n', opts.block_design_regionlist) ;
	end
	assert(length(opts.track_functions)==length(opts.track_files))
	assert(length(opts.track_functions)==length(opts.track_params))
	assert(length(opts.track_functions)==length(opts.track_names))
	assert(length(opts.track_functions)==length(opts.track_monoton_functions))
	fprintf(fd, '\n');
	fprintf(fd, '\tNumber of tracks: \t%i\n', length(opts.track_functions)) ;
	for kk = 1:length(opts.track_functions)
		fprintf(fd, '\t\tUse %s as track_%i\n', opts.track_names{kk}, kk);
		fprintf(fd, '\t\tParser function: \t\t\t%s\n', opts.track_functions{kk});
		fprintf(fd, '\t\tFrom file: \t\t\t\t%s\n', opts.track_files{kk});
		fprintf(fd, '\t\tNumerical parameter: \t\t\t%s\n', opts.track_params{kk});
		fprintf(fd, '\t\tMonotonicity constraints defined by: \t%s\n', opts.track_monoton_functions{kk});
		fprintf(fd, '\n');
	end
	assert(length(opts.segment_feature_functions)==length(opts.segment_feature_files))
	assert(length(opts.segment_feature_functions)==length(opts.segment_feature_params))
	assert(length(opts.segment_feature_functions)==length(opts.segment_feature_names))
	assert(length(opts.segment_feature_functions)==length(opts.segment_feature_monoton_functions))
	fprintf(fd, '\n');
	fprintf(fd, '\tNumber of segment features: \t%i\n', length(opts.segment_feature_functions)) ;
	for kk = 1:length(opts.segment_feature_functions)
		fprintf(fd, '\t\tUse %s as segment_feature_%i\n', opts.segment_feature_names{kk}, kk);
		fprintf(fd, '\t\tReading features from file: \t\t%s\n', opts.segment_feature_files{kk});
		fprintf(fd, '\t\tWith function: \t\t\t\t%s\n', opts.segment_feature_functions{kk});
		fprintf(fd, '\t\tNumerical parameter: \t\t\t%s\n', opts.segment_feature_params{kk});
		fprintf(fd, '\t\tMonotonicity constraints defined by: \t%s\n', opts.segment_feature_monoton_functions{kk});
	end
	fprintf(fd, '\n') ;
end ;

