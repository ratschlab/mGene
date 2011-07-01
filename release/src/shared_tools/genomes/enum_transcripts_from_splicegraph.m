
function [ transcripts, paths ] = enum_transcripts_from_splicegraph( splicegraph, nof_paths );

% --- init
if( ~exist('nof_paths','var') )
  nof_paths = nof_paths_in_splicegraph( splicegraph );
end;
exon_list = splicegraph{1};

% --- trivial splice graph: no exon?
if( isempty(exon_list) )
  transcripts = {};
  paths = {};
  return;
end;

% --- work with lower triagonal adjacency matrix A
A = tril( splicegraph{2} );
n = size( exon_list, 2 );
if( n==1 & isempty(A) )
  A = [ 0 ];
end;
assert_all( size(A) == [n,n] );
assert_all( A==0 | A==1 );

% --- enable enough recursions
RECURSION_LIMIT = get( 0, 'RecursionLimit' );
if( RECURSION_LIMIT < n+12 )
  set( 0, 'RecursionLimit', n+12 );
  warning( sprintf('RecursionLimit set to %d',n+12 ) );
end;

% --- compute lists of outgoing edges
outedges = cell( 1, n );
for( i = 1:n )
  outedges{i} = find( A(:,i)' );
end;

% --- compute paths
paths = cell( 1, nof_paths );
p = 0;
is_first = ( sum(A,2) == 0 );
for( j = find(is_first)' )
  [ paths, p ] = add_transcripts( outedges, paths, p, j );
end;
assert( p == nof_paths );
set( 0, 'RecursionLimit', RECURSION_LIMIT );

% --- convert paths to transcripts
transcripts = cell( size(paths) );
for( p = 1:nof_paths )
  transcripts{p} = exon_list( :, paths{p} )';
  assert( size(transcripts{p},2) == 2 );
end;



function [ paths, p ] = add_transcripts( outedges, paths, p, path );

edges = outedges{ path(end) };
if( isempty(edges) )
  p = p + 1;
  assert( p <= length(paths) );
  paths{p} = path;
else
  for( k = edges )
    [ paths, p ] = add_transcripts( outedges, paths, p, [path,k] );
  end;
end;


