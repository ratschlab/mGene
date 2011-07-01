function genes = count_alt_region(genes)
% function genes = count_alt_region(genes)
%
% counts the number of alternative splicing events in
% 5' UTR, coding region and 3' UTR
%
% Note that alternative transcription starts and ends are also included in this count, but
% intron retention is missed.
%
% Written by: Cheng Soon Ong, 31 May 2007

MAX_EXONS = 100;

for ix = 1:length(genes)
  if mod(ix,1000)==0
    fprintf('.');
  end
  if ~isfield(genes(ix),'coding_region')
    warning('Coding region unknown, run detect_coding_region');
    return
  end

  genes(ix).numalt_5utr = 0;
  genes(ix).numalt_coding = 0;
  genes(ix).numalt_3utr = 0;

  % get all the isoforms
  vertices = genes(ix).splicegraph{1};
  %edges = genes(ix).splicegraph{2};
  %if size(vertices,2) > MAX_EXONS
  %  fprintf('skipping gene %d, with %d exons\n',ix,size(vertices,2));
  %  continue;
  %end
  %strand = genes(ix).strands(1);
  %isoforms = find_all_isoforms(vertices,edges,strand);

  region = count_mrna_level(vertices);
  if isempty(region), continue; end
  assert(all(region(3,:)>1));
  % the number of alternative events is one less than the number of levels.
  region(3,:) = region(3,:) - 1;

  % for each of 5'UTR, coding and 3'UTR, note the alternative regions.
  for ixr = 1:size(region,2)
    if genes(ix).strands(1) == '+'
      if segment_overlap(region(1,ixr),region(2,ixr),genes(ix).start,genes(ix).coding_region(1))
	genes(ix).numalt_5utr = max(genes(ix).numalt_5utr,region(3,ixr));
      end
      if segment_overlap(region(1,ixr),region(2,ixr),genes(ix).coding_region(1),genes(ix).coding_region(2))
	genes(ix).numalt_coding = max(genes(ix).numalt_coding,region(3,ixr));
      end
      if segment_overlap(region(1,ixr),region(2,ixr),genes(ix).coding_region(2),genes(ix).stop)
	genes(ix).numalt_3utr = max(genes(ix).numalt_3utr,region(3,ixr));
      end
    else
      if segment_overlap(region(1,ixr),region(2,ixr),genes(ix).coding_region(2),genes(ix).stop)
	genes(ix).numalt_5utr = max(genes(ix).numalt_5utr,region(3,ixr));
      end
      if segment_overlap(region(1,ixr),region(2,ixr),genes(ix).coding_region(1),genes(ix).coding_region(2))
	genes(ix).numalt_coding = max(genes(ix).numalt_coding,region(3,ixr));
      end
      if segment_overlap(region(1,ixr),region(2,ixr),genes(ix).start,genes(ix).coding_region(1))
	genes(ix).numalt_3utr = max(genes(ix).numalt_3utr,region(3,ixr));
      end
    end
  end
end
