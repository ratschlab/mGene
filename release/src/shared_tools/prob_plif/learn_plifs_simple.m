function  [prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = learn_plifs_simple(Out,L,bins);

% Input\
% 
% Functions called: -> pwl_add

s = sort(Out) ;
limits = [s(round(linspace(1,length(s),bins-1))) inf] ;
 
prob = zeros(bins,1);
prob_cum = zeros(bins,1);
for i=1:bins-1, 
  prob_cum(i) = (sum(Out>=limits(i) & L==1 )+1e-10)/(sum(Out>=limits(i))+1e-10) ;
  prob(i) = (sum(Out>=limits(i) & Out<limits(i+1) & L==1 )+1e-10)/(sum(Out>=limits(i) & Out<limits(i+1))+1e-10) ;
end 
prob_cum(bins) = 1;
prob(bins) = 1;
% if any(prob_cum<prob)
%   keyboard
% end

% prob = sum(XT.*repmat(L==1,bins,1),2)./sum(XT.*repmat(L==-1,bins,1)+XT.*repmat(L==1,bins,1),2);

% hack to make them monotonic increasing
prob_cum_unsort = prob_cum ;
prob_unsort = prob ;
prob_cum = sort(prob_cum) ;
prob = sort(prob) ;

