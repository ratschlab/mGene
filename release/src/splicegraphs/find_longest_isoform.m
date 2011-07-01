function find_longest_isoform

function FloydWarshall ()

  p = zeros(n,n,n) ;
  path = zeros(n,n,n) ;

  for k = 1:n
    for i=1:n,
      for j=1:n,
        path(i,j,k-1) = min( path(i,j,k-1), path(i,k,k-1)+path(k,j,k-1) );
      end
    end
  end 

  minval = inf ;
  mini=nan ;
  minj=nan ;
  for i=1:n,
    for j=1:n,
      if path(i,j,n)<minval, 
        mini=i;
        minj=j ;
        minvak=path(i,j,n) ;
      end ;
    end ;
  end 

  for k = n:-1:1
    if path(mini,minj,k-1) > path(mini,k,k-1)+path(k,minj,k-1) );
  end 
  