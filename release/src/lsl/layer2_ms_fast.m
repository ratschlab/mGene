function layer2_ms_fast(organism, experiment)

for C = [0.01]

	base_dir = get_pred_dir(organism, experiment, 0);
	unix(sprintf('mkdir -p %s', base_dir));
	
	output_dir 		= sprintf('%s/..', base_dir);
	
	load(sprintf('%s/lsl/train/init.mat', output_dir), 'PAR')
	blocks = load_struct(sprintf('%s/lsl/data/training_blocks.mat', output_dir), 'blocks')
	
	PAR.LSL.method.par_ms.C_regul.transitions_sq	=C;
	PAR.LSL.method.par_ms.C_regul.plif_ys_sq		=C;
	PAR.LSL.method.par_ms.C_regul.smoothness_sq		=C*10;
	
	QP_new = initialize_QP(blocks, PAR.LSL.method, PAR.model, 1, PAR.weight_names) ; 

	unix(sprintf('mkdir -p %s/lsl/train', output_dir));

	iteration = 1 ;
	while (1)
		fname = sprintf('%s/lsl/train/train_iteration%i.mat', output_dir,iteration) ;
		fname_tag = sprintf('%s/lsl/train/start_iteration%i_C%f', output_dir,iteration, C) ;
		if ~fexist(fname), 
			iteration=iteration-1 ; 
			fname2 = sprintf('%s/lsl/train/train_iteration%i.mat', output_dir,iteration) ;
			load(fname2, 'TRAIN_PAR', 'QP', 'num_examples', 'obj_solve');
			QP.Q = QP_new.Q;
			save(fname, 'TRAIN_PAR', 'QP', 'num_examples', 'obj_solve');
			unix(sprintf('touch %s', fname_tag));
			break ; 
		 end ;
		iteration=iteration+1 ;
	end ; 
	
	[fn_predictor, iter] = train_path_caller(PAR) ;
end

return

