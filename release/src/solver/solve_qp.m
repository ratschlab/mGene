function [QP frac_solved] = solve_qp(QP, PAR, percentile) ;
%   QP = solve_qp(QP, PAR) ;
%

[engine, env, tmp, tmp, licence_type]=determine_engine() ;

if isequal(engine, 'octave'),

  %[x, obj, info, lambda] = qp (@var{x0}, @var{H}, @var{q}, @var{A}, @var{b}, @var{lb}, @var{ub}, @var{A_lb}, @var{A_in}, @var{A_ub})

  if nnz(QP.Q)==0
    disp('Solving linear program') ;      

    if keyboard_allowed(),
      keyboard
    end ;

    %tic;[res] = mosek_qp(QP.res, QP.Q, QP.f, [], [], QP.lb, QP.ub, [], QP.A(idx,:), QP.b(idx)); toc;

    A_lb = -inf*ones(size(QP.b)) ;
    tic;[res] = mosek_qp(QP.res, QP.Q, QP.f, [], [], QP.lb, QP.ub, A_lb, QP.A(idx,:), QP.b(idx)); toc;

  else
    fprintf('Solving quadratic program. Licence: %s', licence_type) ;
	if nargin==3
		idx = find_relevant_constraints_prc(QP, percentile);
	else
		idx = find_relevant_constraints(QP);
	end
    QP.con_idx = idx;

    %tic;[res] = mosek_qp(QP.res, QP.Q, QP.f, [], [], QP.lb, QP.ub, [], QP.A(idx,:), QP.b(idx)); toc;

    if ~isequal(licence_type, 'commercial')
      A_lb = -inf*ones(size(QP.b(idx))) ;
      tic;[res] = mosek_qp(QP.res, QP.Q, QP.f, [], [], QP.lb, QP.ub, A_lb, QP.A(idx,:), QP.b(idx)); toc;
    else
      % use glpk QP solver
      ctype = char( 'U' * ones(1, length(idx))) ;
      tic;[res] = qpng(QP.Q, QP.f, QP.A(idx,:), QP.b(idx), ctype, QP.lb, QP.ub); toc;

      epsilon = 1e-2;
      if ~(all(QP.A(idx,:)*res-QP.b(idx)<=epsilon))
	num_violated = length(find(QP.A(idx,:)*res-QP.b(idx)>epsilon));
        fprintf('%i (%i) constraints are violated by the solution obtained with glpk (eps: %e)\n', num_violated, length(idx), epsilon);
      end
      %tic; [res, OBJ, INFO, LAMBDA] = qp (QP.res, QP.Q, QP.f, [], [], QP.lb, QP.ub, A_lb, QP.A(idx,:), QP.b(idx)); toc;
    end ;
  end ;
  last_num_con = size(QP.A,1) ;

  %if ~strcmp(how, 'OK')
  %  tic;[res,lambda,how]=lp_resolve(QP.lpenv, last_p_lp, 1, 'bar');toc;
  %end
  QP.res = res ;

  %% check if solution is feasible
  %assert(all(QP.A(idx,:)*res-QP.b(idx)<=1e-2)) ;
  diff = (QP.A*res-QP.b)<1e-2;
  frac_solved = mean(diff);
  fprintf('%f%% of all constraints are satisfied by the solution obtained on a subset of constraints\n', frac_solved*100);

else


  %how=lp_set_param(QP.lpenv,'CPX_PARAM_BARCROSSALG', 2,1) ; % create lp basis after barrier solving
  
  if nnz(QP.Q)==0
    disp('Solving linear program') ;      
    how=lp_set_param(QP.lpenv,'CPX_PARAM_PREDUAL', 1, 1) ; 
    tic;[res,lambda,how] = lp_solve(QP.lpenv, QP.f, QP.A, QP.b, QP.lb, QP.ub, 0, 1, 'dual');toc;
  else
    disp('Solving quadratic program') ;      
	if nargin==3
		idx = find_relevant_constraints_prc(QP, percentile);
	else
    	idx = find_relevant_constraints(QP);
	end
    QP.con_idx = idx;
    tic;[res, lambda, how] = qp_solve(QP.lpenv, QP.Q, QP.f, QP.A(idx,:), QP.b(idx), QP.lb, QP.ub, 0, 0, 'bar');toc;
    
    if ~strcmp(how, 'OK')
      how
    end
    QP.res = res ;
    
    %% check if solution is feasible
    assert(all(QP.A(idx,:)*res-QP.b(idx)<=1e-2)) ;
    diff = (QP.A*res-QP.b)<1e-2;
	frac_solved = mean(diff);
    fprintf('%1.3f%% of all constraints are satisfied by the solution obtained on a subset of constraints\n', frac_solved*100);
  end ;
  
  last_num_con = size(QP.A,1) ;
end ;

return

function idx = find_relevant_constraints_prc(QP, percentile)
	eps = 1e-2;
	Diff = (QP.b-QP.A*QP.res) ;
	% use all violated constraints and all active constraints 
	% plus a percetile of the inactive constraints
	cutoff = prctile(Diff(Diff>eps), percentile);
	idx = find(Diff<=cutoff);
    fprintf('number of violated constraints: %i (Diff<-%g)\n', sum(Diff<-eps), eps) ;
    fprintf('number of active constraints: %i (%g>Diff>-%g)\n', sum(Diff<eps&Diff>-eps), eps, eps) ;
    fprintf('cutoff: %1.4f \n', cutoff) ;
    fprintf('reduced problem has %i constraints (originally %i)\n', length(idx), size(QP.A,1)) ;
return

function idx = find_relevant_constraints(QP)

idx = 1:size(QP.A,1) ;
if 1,
  Diff = (QP.A*QP.res-QP.b) ;
  for p=95:-5:10
    ttt=my_prctile(Diff,p) ;
    if my_prctile(Diff,p)<-50, % some heuristic to make sure the gap is large enough
      idx = find(Diff>my_prctile(Diff,p)) ;
      fprintf('removing %2.0f%% of all constraints (d=%1.2f):\nreduced problem has %i constraints (originally %i)\n', p, my_prctile(Diff,p), length(idx), size(QP.A,1)) ;
      break ;
    end ;
  end ;
end ;


return
