function numsteps = count_steps(cover)
% function numsteps = count_steps(cover)
%
% Count the number of changes in cover, including the size of the step
%
% Written by: Cheng Soon Ong, 31 May 2007

numsteps = 0;
if isempty(cover)
  return
end

baselevel = cover(1);
curlevel = baselevel;
for ix = 1:length(cover)
  if curlevel ~= cover(ix)
    numsteps = numsteps + abs(cover(ix)-curlevel);
    curlevel = cover(ix);
  end
end
numsteps = numsteps + abs(curlevel-baselevel);
numsteps = numsteps/2;
assert(round(numsteps) == numsteps);

