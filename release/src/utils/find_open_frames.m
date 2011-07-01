function frames=find_open_frames(s, terminal)

s=lower(s) ;

if terminal,
  s=s(1:end-1) ;
end ;

s=lower(s) ;

frames = ones(1,3) ;
if length(s)>2,
  idx=[findstr(s,'taa') findstr(s,'tag') findstr(s,'tga')] ;
  frames(mod(idx-1,3)+1)=0 ;
end ;

