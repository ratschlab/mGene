function [max_start,max_end]=find_max_orf(str)

[max_starts,max_ends] = find_max_orfs(str) ;
[tmp,idx]=max(max_ends-max_starts) ;
max_start=max_starts(idx) ;
max_end=max_ends(idx) ;