
function seq = mask_seq( seq, start, from, to );

assert( from >= start );
assert( to-start < length(seq) );
if( from <= to )
  seq( (from-start+1):(to-start+1) ) = '-';
end;


