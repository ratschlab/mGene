function [regions,level] = detectaltregions(gene)
% regions = detectaltregions(gene)
  
[max_exon_level,max_intron_level,exon_level,intron_level,level] = ...
    detectsplicegraph(gene) ;

if all(level(3,:)<=1), 
  regions=[] ;
  return ;
else
  regions = [] ;
  idx=find(level(3,:)>1) ;
  idx=unique([idx(1) find(level(3,2:end)>1 & level(3,1:end-1)<=1)+1]);

  
  for i=1:length(idx)
    %level(4,idx(i))
    %if level(4,idx(i))==0, keyboard; end ;

  %    if idx(i)==1
      start = level(1,idx(i)) ;
%    else
%      start = floor(mean([level(1,idx(i)-1) level(1,idx(i))])) ;
%    end ;
%    if idx(i)==size(level,2)
      stop = level(2,idx(i)) ;
%    else
%      stop = floor(mean([level(2,idx(i)) level(2,idx(i)+1)]))-1 ;
%    end ;
    depth = level(3,idx(i)) ;
    for j=idx(i):size(level,2),
      if level(3,j)>1 
%        if j==size(level,2)
          stop = level(2,j) ;
%        else
%          stop = floor(mean([level(2,j) level(2,j+1)]))-1 ;
%        end ;
        depth = max(depth, level(3,j)) ;
      else
        %if level(4,j)==0, keyboard; end ;
        break ;
      end ;
    end ;
    regions=[regions, [start; stop; depth]] ;
  end ;
end ;