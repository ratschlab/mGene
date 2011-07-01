function yes = match_first_exon(vertices,edges,exon_idx,test_exon_idx)

yes = 0;

%if (vertices(2,exon_idx) == vertices(2,test_exon_idx)) ||...
%    (vertices(2,exon_idx)+1 == vertices(2,test_exon_idx)) ||...
%    (vertices(2,exon_idx)-1 == vertices(2,test_exon_idx))
if (vertices(2,exon_idx) == vertices(2,test_exon_idx))
  yes = 1;

  cur_edge_left = 0;
  idx = 1;
  while (cur_edge_left == 0) && (idx<=size(vertices,2))
    if ((edges(exon_idx,idx)) &&...
	(vertices(2,idx) < vertices(1,exon_idx)))
      cur_edge_left = 1;
    end
  idx = idx + 1;
  end

  if cur_edge_left && ...
	(vertices(1,exon_idx)-5 > vertices(1,test_exon_idx))
    yes = 0;
    return
  end

  

  test_edge_left = 0;
  idx = 1;
  while (test_edge_left == 0) && (idx<=size(vertices,2))
    if ((edges(test_exon_idx,idx)) &&...
	(vertices(2,idx) < vertices(1,test_exon_idx)))
      test_edge_left = 1;
    end
  idx = idx + 1;
  end

  if test_edge_left && ...
	(vertices(1,exon_idx)+5 < vertices(1,test_exon_idx))
    yes = 0;
    return
  end

  if cur_edge_left && test_edge_left
    yes = 0;
  end

end


