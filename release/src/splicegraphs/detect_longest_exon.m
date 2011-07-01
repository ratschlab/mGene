function genes = detect_longest_exon(genes)
% function genes = detect_longest_exon(genes)
%
% sets genes.len_longest_exon
%
% Written by: Cheng Soon Ong, 14 May 2007

for ix = 1:length(genes)
  vertices = genes(ix).splicegraph{1};
  num_verts = size(vertices,2);
  max_len = 0;
  for ixv = 1:num_verts
    exon_len = vertices(2,ixv) - vertices(1,ixv);
    assert(exon_len >= 0);
    max_len = max(max_len, exon_len);
  end
  genes(ix).len_longest_exon = max_len;
end
