function  [Conf, Conf_cum] = Out2Conf(Out,prob, prob_cum,limits);
% [Conf, Conf_cum] = Out2Conf(Output ,prob, prob_cum,limits);
%  
% function transforms outputs of an SVM to probability and cumulative
% probability values by means of piecewise linear functions (PLIF). 
% Use learn_plifs to determine the function (it will return supporting
% points 'limits' and corresponding y-values, 'prob' and 'prob_cum')
%
% INPUT   Out       vector of SVM Output to be transformed
%         prob
%         prob_cum
%         limits    
%
% OUTPUT  Conf      vector (same length as Out) with estimated
%                   probability for true positive
%         Conf_cum  vector (same length as Out) with estimated 
%                   cumulative probability true positive
% 
% Functions called: -> pwl_add
%
% see also: learn_plifs
  
    
bins = length(prob); 

Conf_cum=zeros(1,length(Out)) ;
Conf=zeros(1,length(Out)) ;

for i=1:length(Out)
  x(:,1) = pwl_add(zeros(bins, 1), Out(i), limits) ;
  Conf_cum(i) = x'*prob_cum;
  Conf(i) = x'*prob;
end ;
