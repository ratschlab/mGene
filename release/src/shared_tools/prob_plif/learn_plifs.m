function  [prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = learn_plifs(Out,L,bins,do_balancing,min_inc);
% [prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = learn_plifs_qp(Output,Labels,numberofbins,do_balancing,min_inc);
% 
% learn piecewise linear functions (PLIF), with which outputs of an SVM
% can be transformed to probability and cumulative
% probability values
% 
% INPUT    Out                vector of SVM Outputs on an evaluation set
%
%          L                  vector of corresponding labels (-1,1)  
%
%          bins               number of bins used for interpolation (default 50)
%                             (limits are choosen such that the number of
%                             examples in every bin is similar)  
%
%          do_balancing       should be set to 1 if negative class is
%                             much bigger than positive class,
%                             (negatives are downsampled for
%                             determining the bins) 
%
%          min_inc            probabilities have to be monotonically
%                             increasing (higher SVM output value
%                             corresponds to higher probability of true positive) 
%                             (default 1e-5)
%
%
% OUTPUT   prob               vector of length 'bins'. prob_cum(i) is the
%                             estimated probability of examples with 
%                             Output = limits(i) to be true positive
% 
%          prob_cum           vector of length 'bins'. prob_cum(i) is the
%                             estimated probability of examples with 
%                             Output> limits(i) to be true positive
%
%          limits             vector of length 'bins' with supporting
%                             points (SVM Output values)
%
% see also: Out2Conf  
  
if nargin<5, min_inc = 1e-5 ; end ;
if nargin<4, do_balancing = 1; end ;
if nargin<3, bins = 50; end ;

if do_balancing
  s1 = [min(Out) max(Out)];
  s1 = [s1 Out(L==1)];
  idx = find(L==-1);
  num = min(sum(L==-1),sum(L==1));
  % ii = randperm(length(idx));
  idx = idx(randperm(length(idx)));
  s1 = sort([s1 Out(idx(1:num))]);
  limits = [s1(round(linspace(1,length(s1),bins-1))) inf] ;
else
  s1 = sort(Out) ;
  limits = [s1(round(linspace(1,length(s1),bins-1))) inf] ;
end
  
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

global NO_CPLEX 
if isequal(determine_engine(), 'octave') || isequal(NO_CPLEX,1),
  warning('Do not optimize plifs with QP')

  prob_cum = sort(prob_cum) ;
  prob = sort(prob) ;
  return ;
end ;

addpath /fml/ag-raetsch/share/software/matlab_tools/cplex9/
warning('using cplex')

which cplex_license


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
while ~lpenv
  pause(10)
  lpenv=cplex_license(0) ;
end
[res,lambda,how]=qp_solve(lpenv,Q,f,A,b',lb,ub,0) ;

p1=res(1:bins) ;
p2=res(bins+1:2*bins) ;

prob=p1 ;
prob_cum=p2 ;

cplex_quit(lpenv) ;

return ;

