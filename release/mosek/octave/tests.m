x0 = [.1,.1]' ;

P = 2*[2,0.5; 0.5, 1] ;
q = [1.0, 1.0]' ;

A = [1.0,1.0];
b = [1.0]

lb = -inf * ones(2,1) ;
ub = inf * ones(2,1) ;

A_in = [-1.0,0.0;0.0,-1.0];

A_lb = -inf * ones(2,1) ;
A_ub =  zeros(2,1) ;

M = __mosek_qp__(x0,P,q,A,b,lb,ub,A_lb,A_in,A_ub) ;

N = qp(x0,P,q,A,b,lb,ub,A_lb,A_in,A_ub) ;

M
N

%   [ 2.50e-01]
%   [ 7.50e-01]
