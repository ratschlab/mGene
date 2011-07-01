function res=qpng(Q, f, A, b, ctype, lb, ub);
%Matlab MEX interface for the GLPK library
%
%[xopt, fmin, status, extra] = glpk (C, A, b, lb, ub, ctype, vartype, sense, param);
%
%Solve an LP/MILP problem using the GNU GLPK library. Given three
%arguments, glpk solves the following standard LP:
%
%min C'*x subject to A*x  <= b

error('qpng.m: NOT IMPLEMENTED')

%% assert that Q has entries only on its main diagonal
assert(isequal(Q.*eye(size(Q,1)), Q))
assert(all(all(Q>=0)))
C = max(Q);

% f states the linear regularizer term in the qp case
% so I simply add it here.
% I also add 0.1 here because the smoothness help variables 
% for the transitions were not regularized at all
C = full(C+f');

param.msglev = 3;

sense = 1;%minimization
vartype = repmat('C', 1, size(Q, 1));

[xopt, fmin, status, extra] = glpk (C, A, b, lb, ub, ctype, vartype, sense, param);

if ~(all((lb-xopt)<1e-2))
	fprintf('lower bound of variable violated by the solution obtained with glpk\n');
	fprintf('maximal violation: %i\n', max(xopt-ub));
end

if ~(all((xopt-ub)<1e-2))
	fprintf('upper bound of variable violated by the solution obtained with glpk\n');
	fprintf('maximal violation: %i\n', max(xopt-ub));
end

res = xopt;
