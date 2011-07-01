function [mask hist] = compute_mers(string, order, stepping, offset)
%function [mask hist] = compute_mers(string, order, stepping)

mask    = zeros(1,4^order) ;
hist    = zeros(1,4^order) ;

%if nargin<3,
%  stepping = 1;
%end;
%
%if nargin<4,
%  offset = 0;
%end;
%
%if isempty(string), 
%  return; 
%end;
string = lower(char(string)) ;

if ~isempty(find(~ismember(string,'acgt')))
  string(find(~ismember(string,'acgt'))) = 'a';
end;

assert(order<=8);

if length(string)<order,
  return ;
end
t = sg('translate_string', double(string), order, order-1) ;
if length(t)<order
  return ;
end ;

% transform t to matlab coordinates
t = t(1:end-order+1)+1 ;

% throw away kmers if stepping > 1
assert(stepping>0);
%t = t(1:stepping:end);
t = t(end-offset:-stepping:1);

% reduce t to the information whether or not a certain kmer occurs 
u = unique(t) ;

% get histogram counts
t = sort(t);
[dummy loc] = ismember(u, t);
occ = loc - [0 loc(1:end-1)];

% return histogram
hist(u) = occ/sum(occ);

% return "binary histogram" 
% telling whether or not a certain kmer is present
mask(u) = 1/sqrt(length(u)) ;
