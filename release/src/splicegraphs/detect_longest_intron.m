function genes = detect_longest_intron(genes)
% function genes = detect_longest_introns(genes)
%
% set genes.len_longest_intron.
%
% Written by: Cheng Soon Ong, 14 May 2007.

for ix = 1:length(genes)
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  num_verts = size(vertices,2);
  max_len = 0;
  for ixi = 1:num_verts-1
    for ixj = ixi+1:num_verts
      if edges(ixi,ixj)
        intron_len = vertices(1,ixj) - vertices(2,ixi);
        assert(intron_len >= 0);
	max_len = max(max_len, intron_len);
      end
    end
  end
  genes(ix).len_longest_intron = max_len;
end
