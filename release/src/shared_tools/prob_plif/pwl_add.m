function vec = pwl_add(vec, value, limits)
% vec = pwl_add(vec, value, limits)
  
idx = sum(limits<=value) ;

if idx==0,
  vec(1)=vec(1)+1 ;
elseif idx==length(limits)-1,
  vec(end)=vec(end)+1 ;
else
  vec(idx+1) = vec(idx+1) + (value-limits(idx))/(limits(idx+1)-limits(idx)) ;  
  vec(idx)   = vec(idx)   + (limits(idx+1)-value)/(limits(idx+1)-limits(idx)) ;
end ;