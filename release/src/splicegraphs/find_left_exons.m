function exons_5prime = find_left_exons(vertices, edges, idx, num);
% exons_5prime = find_5prime_exons(vertices, edges, idx[, num=2]);

exons_5prime=[vertices(1,idx) vertices(2,idx)] ;

if nargin<4, num=2 ; end ;
if num<=0, return ; end ;
  
% find 5' exons
edges=triu(edges) ;
idx5=find(edges(:,idx)) ;
if isempty(idx5),
  return ;
end ;

% determine the 5' exons of each of the 5' exons 
% and compute the length
l=[] ; exons={} ;
for i=1:length(idx5),
  exons{i} = find_left_exons(vertices, edges, idx5(i), num-1);
  l(i) = sum(exons{i}(:,2)-exons{i}(:,1)) ;
  assert(l(i)>=0) ;
  assert(all(exons{i}(:,2)<vertices(1,idx))) ;
end ;

% select the largest transcript
[tmp,max_idx]=max(l) ;
exons_5prime = [exons{max_idx};
                exons_5prime] ;




