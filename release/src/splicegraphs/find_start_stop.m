function [tstart,tstop] = find_start_stop(vertices,edges,strand)
% function [tstart,tstop] = find_start_stop(vertices,edges,strand)
%
% find all the leftmost and rightmost exons (transcription starts and stops)
tstart = [];
tstop = [];
for ixv = 1:size(vertices,2)
  if sum(edges(ixv,1:ixv-1)) == 0
    tstart(end+1) = ixv;
  end
  if sum(edges(ixv+1:end,ixv)) == 0
    tstop(end+1) = ixv;
  end
end

if strand == '-'
  temp = tstart;
  tstart = tstop;
  tstop = temp;
end
