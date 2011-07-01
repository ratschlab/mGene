function blocks = count_mrna_level(vertices)
% function blocks = count_mrna_level(vertices)
%
% Counts the number of different exons in the splicegraph that covers a particular position.
% This estimates the number of alternative splicing and transcription events, but misses
% intron retention.
%
% Written by: Cheng Soon Ong, 31 May 2007

gstart =  vertices(1,1) - 1;
gstop =  vertices(2,end) + 1;
local = vertices - gstart;
mrnalevel = zeros(1,gstop-gstart);

for ixv = 1:size(local,2)
  mrnalevel(local(1,ixv):local(2,ixv)) = mrnalevel(local(1,ixv):local(2,ixv)) + 1;
end

% convert the array of levels into blocks, only reporting levels greater than 1
blocks = zeros(3,0);
level = 0;
bstart = -inf;
for ix = 1:length(mrnalevel)
  if mrnalevel(ix)~=level
    if level > 1
      blocks(:,end+1) = [bstart;ix-1;level];
    end
    level = mrnalevel(ix);
    bstart = ix;
  end
end
if ~isempty(blocks)
  blocks(1:2,:) = blocks(1:2,:) + gstart;
end
