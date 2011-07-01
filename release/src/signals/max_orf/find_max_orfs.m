function [max_start,max_end,max_atg]=find_max_orfs(str, atg)
% [max_start,max_end,max_atg]=find_max_orfs(str, atg)

str=lower(str) ;
if ~exist('atg'),
  atg = '' ;
end ;

idx = sort([strfind(str,'tag'), strfind(str,'taa'), strfind(str,'tga')]) ;
idx_atg = sort(strfind(str,atg)) ;

map=zeros(1,length(str)) ;
map(idx)=1 ;

max_len(1:3) = 0 ; max_start(1:3)=0; max_end(1:3)=0 ; max_atg(1:3) = 0 ;

for j=0:2
  %xif j==1, keyboard ; end ;
  last_start=-2 ; 
  for i=idx(mod(idx-1,3)==j)
    if ~isempty(atg)
      p = idx_atg(find(mod(idx_atg-1,3)==j & idx_atg>=last_start, 1, 'first')) ;
      if isempty(p),
        p=i ; % no atg -> no orf
      end ;
    else
      p=last_start+3; 
    end ;
    if i-p>max_len(j+1),
      max_atg(j+1)=p ;
      max_start(j+1)=last_start+3 ;
      max_end(j+1)=i-1 ;
      max_len(j+1)=i-p ;
    end ;
    last_start=i ; 
  end ;
  if isempty(idx(mod(idx-1,3)==j))
    if ~isempty(atg)
      p = idx_atg(find(mod(idx_atg-1,3)==j, 1, 'first')) ;
      if isempty(p), 
        p=length(p) ;
      end ;
    else
      p=-2+3; 
    end ;
    max_atg(j+1)=p ;
    max_start(j+1)=-2+3 ;
    max_end(j+1)=length(str) ;
    max_len(j+1)=length(str)-p ;
    continue ;
  end ;
  if length(str)-i>max_len(j+1),
    if ~isempty(atg)
      p = idx_atg(find(mod(idx_atg-1,3)==j & idx_atg>=i, 1, 'first')) ;
      if isempty(p),
        p=length(str) ;
      end ;
    else
      p=i ; 
    end ;
    max_atg(j+1)=p ;
    max_start(j+1)=i ;
    max_end(j+1)=length(str)-1 ;
    max_len(j+1)=length(str)-p ;
  end ;
end ;
