function  [prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = learn_plifs_qp(Out,L,bins,min_inc);
% [prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = learn_plifs_qp(Out,L,bins,min_inc);

warning('learn_plifs_qp is an outdated function, use learn_plifs instead')  

if nargin<4, min_inc=1e-4 ; end ;
s = sort(Out) ;
limits = [s(round(linspace(1,length(s),bins-1))) inf] ;
 
prob = zeros(bins,1);
prob_cum = zeros(bins,1);
num = zeros(bins,1) ;
num_cum = zeros(bins,1) ;
for i=1:bins-1, 
  num_cum(i) = sum(Out>=limits(i)) ;
  prob_cum(i) = (sum(Out>=limits(i) & L==1 )+1e-10)/(sum(Out>=limits(i))+1e-10) ;
  num(i) = sum(Out>=limits(i) & Out<limits(i+1)) ;
  prob(i) = (sum(Out>=limits(i) & Out<limits(i+1) & L==1 )+1e-10)/(sum(Out>=limits(i) & Out<limits(i+1))+1e-10) ;
end 
prob_cum(bins) = 1;
prob(bins) = 1;

prob_unsort=prob ;
prob_cum_unsort=prob_cum ;

% prob
w1=num/sum(num) ;
Q1=spdiag(w1) ;
f1 = -prob.*w1 ;

% cum_prob
w2=num_cum/sum(num) ;
Q2=spdiag(w2) ;
f2 = -prob_cum.*w2 ;

Q=[Q1          zeros(bins)
   zeros(bins) Q2 ] ;
f=[f1; f2] ;

A=spzeros(1,2*bins) ;
b=[] ; q=1 ;

% monotonicity
for i=1:bins-1 ;
  A(q,i)=1 ; A(q,i+1)=-1 ; b(q)=-min_inc ; q=q+1 ;
  A(q,bins+i)=1 ; A(q,bins+i+1)=-1 ; b(q)=-min_inc ; q=q+1 ;
end ;

% prob<=prob_cum
for i=1:bins ;
  A(q,i)=1 ; A(q,bins+i)=-1 ; b(q)=0 ; q=q+1 ;
end ;

lb=zeros(2*bins,1) ;
ub=ones(2*bins,1) ;
lpenv=cplex_license(0) ;
[res,lambda,how]=qp_solve(lpenv,Q,f,A,b',lb,ub,0) ;

p1=res(1:bins) ;
p2=res(bins+1:2*bins) ;

prob=p1 ;
prob_cum=p2 ;

cplex_quit(lpenv) ;

return ;

