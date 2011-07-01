function [contig_pos,score] = interval_query_wrapper(fname, scorename, from, to) 

if nargin==2
	range = [0; 1e10];
elseif
	range = [from; to];
end

[contig_pos,score] = interval_query(fname,{scorename},range);
