function dna_pos=map_rna_pos(exons, strand, rna_pos, offset) ;
% dna_pos=map_rna_pos(exons, strand, rna_pos[, offset]) ;

if nargin<4,
  offset=0 ;
end ;

if strand=='+',

  l=exons(:,2)-exons(:,1) ;
  L=[0; cumsum(l)] ; 
  L(end)=[] ;
  idx=find(rna_pos>L+offset, 1, 'last') ;
  rp=rna_pos-L(idx) ;
  assert(rp>=0) ;

  dna_pos = exons(idx,1)+rp-1 ;
else
  l=exons(end:-1:1,2)-exons(end:-1:1,1) ;
  L=[0; cumsum(l)] ;
  L(end)=[] ;
  idx = find(rna_pos>L+offset, 1, 'last') ;
  rp = rna_pos-L(idx) ;
  assert(rp>=0) ;
  
  dna_pos = exons(size(exons,1)+1-idx,2) - rp + 1;
end ;