function found = issubset(A,B)
% function found = issubset(A,B)
%
% returns true if A is a subset of B, where 
% both A and B are vectors with 1 where elements exists and 0 otherwise.
found = isequal(find(A),intersect(find(A),find(B)));