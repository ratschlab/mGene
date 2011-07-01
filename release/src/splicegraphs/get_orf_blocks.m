function [blocks, rfpos, tstart] = get_orf_blocks(isostartpos,isostoppos,isoform,vertices)
% function [blocks, rfpos, tstart] = get_orf_blocks(isostartpos,isostoppos,isoform,vertices)
%
% Return the open reading frame in genomic coordinates.
% Assumes that there are no extra exons in isoform. This is trimmed by get_longest_orf.m.
% rfpos takes values of [0,1,2], NaN indicates non-mRNA.
%
% Written by: Cheng Soon Ong, 30 May 2007

% trim the exons at the ends of the isoforms by isostartpos and isostoppos
if length(isoform) > 2
  blocks = [[isostartpos;vertices(2,isoform(1))],...
	    vertices(:,isoform(2:end-1)),...
	    [vertices(1,isoform(end));isostoppos]];
elseif length(isoform) > 1
  blocks = [[isostartpos;vertices(2,isoform(1))],...
	    [vertices(1,isoform(end));isostoppos]];
else
  blocks = [isostartpos;isostoppos];
end

tstart =  blocks(1,1) - 1;
local = blocks - tstart;
mrna_length = sum(local(2,:)-local(1,:)+1);
if isostartpos>vertices(1,isoform(1)) && isostoppos<vertices(2,isoform(end))
  assert(mod(mrna_length,3)==0);
end
template = repmat([0,1,2],1,ceil(mrna_length/3)+1);

% compute the reading frame for each position
rfpos = repmat(nan,1,isostoppos-isostartpos+2);
if isostartpos>vertices(1,isoform(1))
  ixt = 1;
  for ix = 1:size(local,2)
    blen = local(2,ix)-local(1,ix);
    rfpos(local(1,ix):local(2,ix)) = template(ixt:ixt+blen);
    ixt = ixt + blen + 1;
  end
  if isostoppos<vertices(2,isoform(end))
    assert(template(ixt)==0);
  end
elseif isostoppos<vertices(2,isoform(end))
  ixt = length(template);
  for ix = size(local,2):-1:1
    blen = local(2,ix)-local(1,ix);
    rfpos(local(1,ix):local(2,ix)) = template(ixt-blen:ixt);
    ixt = ixt - blen - 1;
  end
  if isostartpos>vertices(1,isoform(1))
    assert(template(ixt)==2);
  end
end
