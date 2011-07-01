function [max_start,max_end]=find_max_orfs(str, tmp)

idx=sort([strfind(str,'tag'), strfind(str,'taa'), strfind(str,'tga')]) ;

map=zeros(1,length(str)) ;
map(idx)=1 ;

max_len(1:3) = 0 ; max_start(1:3)=0; max_end(1:3)=0 ;

for j=0:2
  last_start=-2 ; 
  for i=idx(mod(idx-1,3)==j)
    if i-last_start>max_len(j+1),
      max_start(j+1)=last_start ;
      max_end(j+1)=i-1 ;
      max_len(j+1)=i-last_start ;
    end ;
    last_start=i ; 
  end ;
  if isempty(idx(mod(idx-1,3)==j))
    max_start(j+1)=-2 ;
    max_end(j+1)=length(str) ;
    max_len(j+1)=length(str) ;
    continue ;
  end ;
  if length(str)-i>max_len(j+1),
    max_start(j+1)=i ;
    max_end(j+1)=length(str)-1 ;
    max_len(j+1)=length(str)-i ;
  end ;
end ;
max_start=max_start+3 ;
