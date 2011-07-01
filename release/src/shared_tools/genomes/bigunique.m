function f=bigunique(s) ;
% f=bigunique(s) ;

f=[] ; step=ceil(length(s)/10) ;
for i=1:10,
  idx=1+(i-1)*step:min(length(s), i*step) ;
  f=unique([f unique(s(idx))]) ;
end ;

