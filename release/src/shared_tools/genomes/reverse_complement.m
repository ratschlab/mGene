function x=reverse_complement(y)

if length(y)>200*1024*1024
  x=uint8(y) ;
  for i=1:length(y),
    x(i)=y(end-i+1) ;
  end ;
  y=x ;
else
  y=y(end:-1:1) ;
  x=y ;
end ;

x(y=='a')='t' ;
x(y=='c')='g' ;
x(y=='g')='c' ;
x(y=='t')='a' ;
x(y=='n')='n' ;
x(y=='A')='T' ;
x(y=='C')='G' ;
x(y=='G')='C' ;
x(y=='T')='A' ;
x(y=='N')='N' ;

x=char(x) ;
