function yes = match_exon_ends(vertices,edges,exon_idx,test_exon_idx)

yes = 0;

cur_edge_left = 0;
idx = 1;
while (cur_edge_left == 0) && (idx<=size(vertices,2))
  if ((edges(exon_idx,idx)) &&...
      (vertices(2,idx) < vertices(1,exon_idx)))
    cur_edge_left = 1;
  end
  idx = idx + 1;
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



cur_edge_right = 0;
idx = 1;
while (cur_edge_right == 0) && (idx<=size(vertices,2))
  if ((edges(exon_idx,idx)) &&...
      (vertices(1,idx) > vertices(2,exon_idx)))
    cur_edge_right = 1;
  end
  idx = idx + 1;
end


test_edge_right = 0;
idx = 1;
while (test_edge_right == 0) && (idx<=size(vertices,2))
  if ((edges(test_exon_idx,idx)) &&...
      (vertices(1,idx) > vertices(2,test_exon_idx)))
    test_edge_right = 1;
  end
  idx = idx + 1;
end


%disp([cur_edge_left,cur_edge_right,test_edge_left, ...
%      test_edge_right]);
%disp([vertices(:,exon_idx)',vertices(:,test_exon_idx)']);


if (cur_edge_left&&~cur_edge_right&&~test_edge_left&&test_edge_right)&&...
      (vertices(2,exon_idx) > vertices(1,test_exon_idx))&&...
      (vertices(1,exon_idx) < vertices(2,test_exon_idx))&&...
      (vertices(1,exon_idx) < vertices(1,test_exon_idx))&&...
      (vertices(2,exon_idx) < vertices(2,test_exon_idx))
  yes = 1;
end


if (~cur_edge_left&&cur_edge_right&&test_edge_left&& ~test_edge_right)&&...
      (vertices(1,exon_idx) < vertices(2,test_exon_idx))&&...
      (vertices(2,exon_idx) > vertices(1,test_exon_idx))&&...
      (vertices(1,exon_idx) > vertices(1,test_exon_idx))&&...
      (vertices(2,exon_idx) > vertices(2,test_exon_idx))
  yes = 1;
end

