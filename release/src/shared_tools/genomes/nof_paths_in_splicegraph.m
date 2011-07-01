
function [ nof_paths ] = nof_paths_in_splicegraph( splicegraph );

% --- get splice graph
exon_list = splicegraph{1};
A = tril( splicegraph{2} );
n = size( exon_list, 2 );
assert_all( size(A) == [n,n] );
assert_all( A==0 | A==1 );

% --- compute number of paths
is_first = ( sum(A,2) == 0 );
is_last = ( sum(A',2) == 0 );
t = double( is_first );
for( j = find(t'==0) )
  t(j) = sum( t( find(A(j,:)) ) );
end;
nof_paths = sum( t(is_last) );

