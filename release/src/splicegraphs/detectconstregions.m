function [regions, level] = detectconstregions(gene)

% [regions, level] = detectconstregions(gene)
%
% Written by Georg Zeller, MPI Tuebingen, 2008
% Modified from detectaltregions.m

[max_exon_level,max_intron_level,exon_level,intron_level,level] = ...
    detectsplicegraph(gene);

const_blocks(1,:) = find(level(3,:)==1 & [1 level(3,1:end-1)~=1]);
const_blocks(2,:) = find(level(3,:)==1 & [level(3,2:end)~=1   1]);

regions = zeros(3,0);
for i=1:size(const_blocks,2),
  start = level(1,const_blocks(1,i));
  stop = level(2,const_blocks(2,i));
  depth = max(level(3,const_blocks(1,i):const_blocks(2,i)));
  regions = [regions, [start; stop; depth]];
  end
end