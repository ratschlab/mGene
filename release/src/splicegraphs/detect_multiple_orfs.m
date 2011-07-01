function genes = detect_multiple_orfs(genes)
% function genes = detect_multiple_orfs(genes)
%
% Detect regions where there is more than one coding frame.
%
% Written by: Cheng Soon Ong, 14 May 2007

MAX_EXONS = 100;

for ix=1:length(genes)
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
  cutisoforms = {};

  % get all the coding regions.
  isostartpos = zeros(1,length(isoforms));
  isostoppos = zeros(1,length(isoforms));
  isoORFlength = zeros(1,length(isoforms));
  for ixf = 1:length(isoforms)
    str = load_mrna(genes(ix).chr,strand,vertices(1,isoforms{ixf})',vertices(2,isoforms{ixf})');
    [strstart,strstop,isoORFlength(ixf)] = get_longest_orf(str,strand);
    [isostartpos(ixf),isostoppos(ixf),cutisoforms{ixf}] = strloc2genloc(isoforms{ixf}, vertices, strstart, strstop);
  end

  % for each pair of isoforms, get the overlap length and the frame shift
  mult_orf = 0;
  for ix1 = 1:length(isoforms)
    for ix2 = ix1+1:length(isoforms)
      [blocks1,rfpos1,tstart1] = get_orf_blocks(isostartpos(ix1),isostoppos(ix1),cutisoforms{ix1},vertices);
      [blocks2,rfpos2,tstart2] = get_orf_blocks(isostartpos(ix2),isostoppos(ix2),cutisoforms{ix2},vertices);
      mult_orf(end+1) = numframeshifted(rfpos1,rfpos2,tstart1,tstart2);
    end
  end
  genes(ix).multiple_frame_lengths = mult_orf(find(mult_orf));
end
