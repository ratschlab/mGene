%
%

n=3;

x0 = [1.0,1.0,1.0]' ;
P = [ 2,0,-1; 0,0.2,0; -1,0,2 ]
q = [ 0,-1,0]';
A = [1,0,1;0,1,0] ;
b = [1,1]' ;

lb = -inf * ones(n,1) ;
ub =  inf * ones(n,1) ;

num_c = 4 ;

A_lb = -inf * ones(num_c,1) ;
A_in =  zeros(num_c,n)
A_ub =  inf * ones(num_c,1) ;

M = __mosek_qp__(x0,P,q,A,b) ;
M = __mosek_qp__(x0,P,q,A,b,lb,ub) ;
M = __mosek_qp__(x0,P,q,A,b,lb,ub,A_lb,A_in,A_ub) ;

M = mosek_qp(x0,P,q,A,b) ;
M = mosek_qp(x0,P,q,A,b,lb,ub) ;

##[x, obj, INFO, lambda] = mosek_qp(x0,P,q,A,b) ;
##[x, obj, INFO, lambda] = mosek_qp(x0,P,q,A,b,lb,ub) ;
##[x, obj, INFO, lambda] = mosek_qp(x0,P,q,A,b,lb,ub,A_lb,A_in,A_ub) ;


