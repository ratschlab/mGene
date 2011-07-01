function exons_3prime = find_right_exons(vertices, edges, idx, num);
% exons_3prime = find_3prime_exons(vertices, edges, idx[, num=2]);
  
exons_3prime=[vertices(1,idx) vertices(2,idx)] ;

if nargin<4, num=2 ; end ;
if num<=0, return ; end ;

% find 3' exons
edges=triu(edges) ;
idx3=find(edges(idx,:)) ;
if isempty(idx3),
  return ;
end ;

% determine the 3' exons of each of the 3' exons 
% and compute the length
l=[] ;
for i=1:length(idx3),
  exons{i} = find_right_exons(vertices, edges, idx3(i), num-1);
  l(i) = sum(exons{i}(:,2)-exons{i}(:,1)) ;
  assert(l(i)>=0) ;
  assert(all(exons{i}(:,1)>vertices(2,idx))) ;
end ;

% select the largest transcript
[tmp,max_idx]=max(l) ;
exons_3prime = [exons_3prime;
                exons{max_idx}] ;




