function genes = detect_coding_region(genes)
% function genes = detect_coding_region(genes)
%
% get the longest open reading frame over all isoforms
% Written by: Cheng Soon Ong

MAX_EXONS = 100;
for ix = 1:length(genes)
  if mod(ix,1000)==0
    fprintf('.');
  end

  % get all the isoforms
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  if size(vertices,2) > MAX_EXONS
    fprintf('skipping gene %d, with %d exons\n',ix,size(vertices,2));
    continue;
  end
  strand = genes(ix).strands(1);
  isoforms = find_all_isoforms(vertices,edges,strand);

  % annotate the gene with the longest one.
  isostartpos = zeros(1,length(isoforms));
  isostoppos = zeros(1,length(isoforms));
  isoORFlength = zeros(1,length(isoforms));
  for ixf = 1:length(isoforms)
    str = load_mrna(genes(ix).chr,strand,vertices(1,isoforms{ixf})',vertices(2,isoforms{ixf})');
    [strstart,strstop,isoORFlength(ixf)] = get_longest_orf(str,strand);
    [isostartpos(ixf),isostoppos(ixf),isoforms{ixf}] = strloc2genloc(isoforms{ixf}, vertices, strstart, strstop);
  end
  [len,idx_longest] = max(isoORFlength);
  genes(ix).coding_length = len;
  genes(ix).coding_region = [isostartpos(idx_longest),isostoppos(idx_longest)];
  genes(ix).coding_isoform = isoforms{idx_longest};
end

