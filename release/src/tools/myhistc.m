
function [ n, bin ] = histc( x, edges );

assert( all( edges == edges ) );

last = repmat( edges(end), size(x) );
bin = zeros( size(x) );
for( i = 1:length(edges) )
  edge = repmat( edges(i), size(x) );
  bin = bin + 1 * ( edge <= x & x <= last );
end;

% --- the following code fails due to
%  /fml/ag-raetsch/share/software/matlab-7.6/toolbox/matlab/datafun/hist.m
% ---
%n = hist( bin, 0:length(edges) );
%n(1) = [];
% --- end of failing code ---

n = zeros( size(edges) );
for( i = 1:length(edges) )
  n(i) = sum( bin == i );
end;

