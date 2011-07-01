function lpenv = init_solver() ;
%  lpenv = init_solver() ;
%

[engine,env]=determine_engine() ;

if isequal(engine, 'octave'),

  lpenv = inf ;

else

global last_p_lp last_num_con
last_p_lp = -1 ;
%last_p_lp=[] ;
%last_num_con=0 ;

for j = 3:-1:1
	lpenv = cplex_license(0,j) ;
	if lpenv == 0
		pause(30)
	else 
		break;
	end
end
if lpenv == 0,
  pause(100) ;
  rproc_rerun('could not obtain cplex license') ; 
end ;

how = lp_set_param(lpenv,'CPX_PARAM_AGGIND', 1,1) ;
how = lp_set_param(lpenv,'CPX_PARAM_THREADS', 2,1) ;
how = lp_set_param(lpenv,'CPX_PARAM_BARTHREADS', 2,1) ;
disp('trying 4 threads') ;
how = lp_set_param(lpenv,'CPX_PARAM_THREADS', 4,1) ;
how = lp_set_param(lpenv,'CPX_PARAM_BARTHREADS', 4,1) ;
disp('trying 8 threads') ;
how = lp_set_param(lpenv,'CPX_PARAM_THREADS', 8,1) ;
how = lp_set_param(lpenv,'CPX_PARAM_BARTHREADS', 8,1) ;
how = lp_set_param(lpenv,'CPX_PARAM_PRECOMPRESS', 1,1) ;
% how = lp_set_param(lpenv,'CPX_PARAM_SIFTALG', 1,1) ; % dual
% how = lp_set_param(lpenv,'CPX_PARAM_SIFTALG', 4,1) ; % Barrier
how = lp_set_param(lpenv,'CPX_PARAM_SIFTALG', 3,1) ; 
how = lp_set_param(lpenv,'CPX_PARAM_SIFTDISPLAY', 2,1) ; % verbose

end
