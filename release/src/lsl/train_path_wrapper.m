function train_path_wrapper(fname)
% train_path_wrapper(fname)

load(fname, 'inventory') ;
for i=1:length(inventory),
  load(fname, inventory{i}) ;
end ;

[fn_predictor, iter] = train_path(PAR) ;

%inventory={'fn_predictor', 'iter'} ;
%save(fname, 'inventory', 'fn_predictor', 'iter', '-V7') ;
