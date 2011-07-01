function [found,h] = find_hash(alph)
% handle = find_hash(alph)
% returns h(x) hash function or ID(x)
% searches a perfect hash function for the elements in (1xn) matrix alph to the interval (0:n-1).
  sub=min(alph);
  alph=alph-sub;
  found=false;
  asize=length(alph);
  fact=asize;
  h= @(x) x;
  for i=1:100
    for j=asize:i*max(alph)
        tmp=mod(mod(i*alph,j),asize);
        contained=zeros(1,asize);
        for k=1:asize
          contained(tmp(k)+1)=1;
        end
        if prod(contained)==1
          found=true;
          h = @(x) mod(mod(i*(x-sub),j),asize);
          return;
        end
    end
  end
  warning('no appropriate hash function found');
end