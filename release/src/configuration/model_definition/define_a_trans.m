function a_trans = define_a_trans(A, transition_pointers, seg_links)

a_trans = zeros(4,sum(~isinf(A(:)))) ;
% [from,to,val,segment]
k = 0 ;
for i=1:size(A,1)
  idx = find(~isinf(A(i,:))) ;
  val = A(i,idx) ;
  pen = transition_pointers(idx,i,1) ;
  a_trans(1,k+1:k+length(idx))=i-1 ;
  a_trans(2,k+1:k+length(idx))=idx-1 ;  
  a_trans(3,k+1:k+length(idx))=val ;    
  for q=1:length(pen) 
    idx_=find(seg_links(:,1)==pen(q)) ;
    if isempty(idx_), 
      assert(pen(q)==0) ; 
    else
      assert(length(idx_)==1) ;
      pen(q) = seg_links(idx_,2) ;
    end ;
  end ;
  a_trans(4,k+1:k+length(idx))= pen;
  k=k+length(idx) ;
end ;

return
