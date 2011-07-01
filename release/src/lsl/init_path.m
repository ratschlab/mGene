function [QP,parameters]= init_path(PAR, blocks, lpenv, weight_names)
% [blocks,QP,parameters]= init_path(PAR)
%

%addpath /fml/ag-raetsch/share/software/matlab_tools/cplex10


%blocks = load_struct(PAR.FN.input_lsl.fn_training_blocks, 'blocks') ;


if fexist(PAR.FN.output_lsl.fn_init) & ~PAR.LSL.method.ignore_init_file
  disp('loading initialization data') ;
  load(PAR.FN.output_lsl.fn_init, 'QP') ;
  Old = load(PAR.FN.output_lsl.fn_init, 'PAR') ;
  assert(isequal(PAR.model,Old.PAR.model)) 
  assert(isequal(PAR.LSL.method,Old.PAR.LSL.method))
else
  % initialize matrices for quadratic program
  QP = initialize_QP(blocks, PAR.LSL.method, PAR.model, 1, weight_names) ;
  if ~fexist(fileparts(PAR.FN.output_lsl.fn_init))
    unix(sprintf('mkdir -p %s', fileparts(PAR.FN.output_lsl.fn_init)));
  end 
  save(PAR.FN.output_lsl.fn_init, 'QP', 'PAR', '-v7') ;
end ;


QP.res = zeros(size(QP.lb)) ;
QP.lpenv = lpenv;
QP = solve_qp(QP, PAR) ;

OBJ = sum(QP.res.*QP.f)+0.5*QP.res'*QP.Q*QP.res ;
parameters = init_weights(PAR.model) ;
parameters = vector2weights(QP.res(1:PAR.model.cnt_parameters)', parameters, PAR.weight_names) ;


