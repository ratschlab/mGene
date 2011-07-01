function genes = count_alt_region_byintron(genes)
% function genes = count_alt_region_byintron(genes)
%
% Detect which regions are alternative by considering overlapping introns,
% and introns overlapping with exons.
%
% Written by: Cheng Soon Ong, 29 May 2007

MAX_EXONS = 100;

for ix=1:length(genes)
  if mod(ix,1000)==0
    fprintf('.');
  end
  if ~isfield(genes(ix),'coding_region')
    warning('Coding region unknown, run detect_coding_region');
    return
  end

  % convert everything to local coordinates
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  numcovered = zeros(1,genes(ix).stop-genes(ix).start+2);
  vertices = vertices - genes(ix).start + 1;
  coding_region = genes(ix).coding_region - genes(ix).start + 1;

  % for each pair of introns, if they overlap and are not the same, note region
  [ixstart,ixend] = find(triu(edges));
  for ix1 = 1:length(ixstart)
    intron1 = [vertices(2,ixstart(ix1))+1:vertices(1,ixend(ix1))-1];
    for ix2 = ix1+1:length(ixstart)
      intron2 = [vertices(2,ixstart(ix2))+1:vertices(1,ixend(ix2))-1];
      overlap = intersect(intron1,intron2);
      if ~isequal(intron1,intron2) && ~isempty(overlap)
	diffintron = setxor(intron1,intron2);
	numcovered(diffintron) = numcovered(diffintron) + 1;
	%numcovered(overlap) = numcovered(overlap) + 1;
      end
    end
  end

  % for each intron, if it is contained withing an  exon, 
  % or if the exon is contained in the intron,
  % then note region.
  for ixi = 1:length(ixstart)
    intron = [vertices(2,ixstart(ixi))+1:vertices(1,ixend(ixi))-1];
    for ixe = 1:size(vertices,2)
      exon = [vertices(1,ixe):vertices(2,ixe)];
      overlap = intersect(intron,exon);
      if isequal(overlap,intron) || isequal(overlap,exon)
	numcovered(overlap) = numcovered(overlap) + 1;
      end
    end
  end
  % for each of 5'UTR, coding and 3'UTR, get max of numcovered
  if genes(ix).strands(1) == '+'
    genes(ix).numalt_5utr = count_steps(numcovered(1:coding_region(1)-1));
    genes(ix).numalt_coding = count_steps(numcovered(coding_region(1):coding_region(2)));
    genes(ix).numalt_3utr = count_steps(numcovered(coding_region(2)+1:end));
  else
    genes(ix).numalt_3utr = count_steps(numcovered(1:coding_region(1)-1));
    genes(ix).numalt_coding = count_steps(numcovered(coding_region(1):coding_region(2)));
    genes(ix).numalt_5utr = count_steps(numcovered(coding_region(2)+1:end));
  end
end
