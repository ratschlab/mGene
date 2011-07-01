
function seq = load_genomic_a( contig, start, stop );

assert( start < stop );
assert( (start<0) == (stop<0) );
strand = 2*(start>0) - 1;

if( strand == +1 )
  seq = load_genomic( contig, '+', start, stop );
end;

if( strand == -1 )
  seq = load_genomic( contig, '-', -stop, -start );
end;

