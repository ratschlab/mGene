function [regions,N]=get_connected_regions(bw);
  
 % similar to bw_label just for 1xm matrix 

reg_starts = find(bw(1:end)-[0 bw(1:end-1)]==1);
reg_ends = find(bw(1:end)-[bw(2:end) 0]==1);
assert(length(reg_starts)==length(reg_ends))

N = length(reg_starts);

regions = zeros(size(bw));
for i=1:N
  assert(reg_starts(i)<=reg_ends(i))
  regions(reg_starts(i):reg_ends(i))=i;
end