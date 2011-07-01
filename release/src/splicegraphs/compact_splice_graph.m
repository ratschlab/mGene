


for gene_idx = 1:length(genes)
  vertices = genes(gene_idx).splicegraph{1};
  edges = genes(gene_idx).splicegraph{2};
  %fprintf(1,'gene index:\t%d\n',gene_idx);
  if (mod(gene_idx,100)==0)
    fprintf(1,'.');
  end
  
  if (size(edges,1) < 2)
    changed = 1;
    while changed
      changed = 0;
      exon_idx = 1;
      while exon_idx <= size(vertices,2)
	
	test_exon_idx = exon_idx+1;
	while test_exon_idx <= size(vertices,2)
	  if (vertices(1,exon_idx)<=vertices(1,test_exon_idx)) &&...
		(vertices(2,exon_idx)>=vertices(2,test_exon_idx))
	    vertices(1,exon_idx) = min(vertices(1,exon_idx),...
				       vertices(1,test_exon_idx));
	    vertices(2,exon_idx) = max(vertices(2,exon_idx),...
				       vertices(2,test_exon_idx));
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	  
	    %disp(vertices)
	    
	    %test_exon_idx = exon_idx+1;
	    changed = 1;
	    
	  end
	  test_exon_idx = test_exon_idx+1;
	end
	exon_idx = exon_idx + 1;
      end
    end
    
  else
    changed = 1;
    while changed
      changed = 0;
      exon_idx = 1;
      while exon_idx <= size(vertices,2)
	
	test_exon_idx = exon_idx + 1;
	while test_exon_idx <= size(vertices,2)
	  
	  
	  
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

	  

	  
	  % 1101 1110
	  if ((cur_edge_left&&cur_edge_right...
		   &&~test_edge_left&&test_edge_right)&&...
		  (vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
		  (vertices(1,exon_idx)-5<=vertices(1,test_exon_idx)))||...
		((cur_edge_left&&cur_edge_right...
		  &&test_edge_left&&~test_edge_right)&&...
		 (vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		 (vertices(2,exon_idx)+5>=vertices(2,test_exon_idx)))
	    
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	    changed = 1;

	  % 1011 0111  
	  elseif ((cur_edge_left&&~cur_edge_right...
		   &&test_edge_left&&test_edge_right)&&...
		  (vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		  (vertices(2,exon_idx)-5<=vertices(2,test_exon_idx)))||...
		((~cur_edge_left&&cur_edge_right...
		  &&test_edge_left&&test_edge_right)&&...
		 (vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
		 (vertices(1,exon_idx)+5>=vertices(1,test_exon_idx)))

	    vertices(:,exon_idx) = vertices(:,test_exon_idx);
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    changed = 1;
	    
	    
	  % 0110  
	  elseif ((~cur_edge_left&&cur_edge_right...
		  &&test_edge_left&&~test_edge_right)&&...
		 (vertices(1,exon_idx)+5>=vertices(1,test_exon_idx))&&...
		 (vertices(2,exon_idx)>=vertices(2,test_exon_idx)-5)&&...
		 (vertices(1,exon_idx)<=vertices(2,test_exon_idx)))


	    vertices(1,exon_idx) = vertices(1,test_exon_idx);
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    changed = 1;

	  % 1001
	  elseif ((cur_edge_left&&~cur_edge_right...
		  &&~test_edge_left&&test_edge_right)&&...
		 (vertices(1,exon_idx)<=vertices(1,test_exon_idx)+5)&&...
		 (vertices(2,exon_idx)-5<=vertices(2,test_exon_idx))&&...
		 (vertices(2,exon_idx)>=vertices(1,test_exon_idx)))


	    vertices(2,exon_idx) = vertices(2,test_exon_idx);
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    

	    changed = 1;

	  

	    
	  
	  % 0101 1010 0000
	  elseif ((vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		 (vertices(2,exon_idx)==vertices(2,test_exon_idx)))||...
		((~cur_edge_left&&cur_edge_right...
		  &&~test_edge_left&&test_edge_right)&&...
		 (vertices(2,exon_idx)==vertices(2,test_exon_idx)))||...
		((cur_edge_left&&~cur_edge_right...
		  &&test_edge_left&&~test_edge_right)&&...
		 (vertices(1,exon_idx)==vertices(1,test_exon_idx)))||...
		(~cur_edge_left&&~cur_edge_right&&...
		  ~test_edge_left&&~test_edge_right)

		
	    
	    vertices(1,exon_idx) = min(vertices(1,exon_idx),...
				       vertices(1,test_exon_idx));
	    vertices(2,exon_idx) = max(vertices(2,exon_idx),...
				       vertices(2,test_exon_idx));
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    
	    changed = 1;

	    
	  % 1100 0011
	  elseif ((cur_edge_left&&cur_edge_right...
		   &&~test_edge_left&&~test_edge_right)&&...
		  (vertices(1,exon_idx)<=vertices(1,test_exon_idx)+5)&&...
		  (vertices(2,exon_idx)>=vertices(2,test_exon_idx)-5))||...
		((~cur_edge_left&&~cur_edge_right...
		  &&test_edge_left&&test_edge_right)&&...
		 (vertices(1,exon_idx)+5>=vertices(1,test_exon_idx))&&...
		 (vertices(2,exon_idx)-5<=vertices(2,test_exon_idx)))

	    vertices(1,exon_idx) = max(vertices(1,exon_idx),...
				       vertices(1,test_exon_idx));
	    vertices(2,exon_idx) = min(vertices(2,exon_idx),...
				       vertices(2,test_exon_idx));
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    changed = 1;
	    
	    
	  % 0010
	  elseif (~cur_edge_left&&~cur_edge_right...
		   &&test_edge_left&&~test_edge_right)&&...
		  (vertices(1,exon_idx)+5>=vertices(1,test_exon_idx))&&...
		  (vertices(1,exon_idx)<=vertices(2,test_exon_idx))

	    vertices(1,exon_idx) = vertices(1,test_exon_idx);
	    vertices(2,exon_idx) = max(vertices(2,exon_idx),...
				       vertices(2,test_exon_idx));
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    changed = 1;
	    
	  % 0100
	  elseif (~cur_edge_left&&cur_edge_right...
		   &&~test_edge_left&&~test_edge_right)&&...
		  (vertices(2,exon_idx)>=vertices(2,test_exon_idx)-5)&&...
		  (vertices(1,exon_idx)<=vertices(2,test_exon_idx))

	    vertices(1,exon_idx) = min(vertices(1,exon_idx),...
				       vertices(1,test_exon_idx));
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    changed = 1;
	    

	  % 0001
	  elseif (~cur_edge_left&&~cur_edge_right...
		   &&~test_edge_left&&test_edge_right)&&...
		  (vertices(2,exon_idx)-5<=vertices(2,test_exon_idx))&&...
		  (vertices(2,exon_idx)>=vertices(1,test_exon_idx))

	    vertices(1,exon_idx) = min(vertices(1,exon_idx),...
				       vertices(1,test_exon_idx));
	    vertices(2,exon_idx) = vertices(2,test_exon_idx);
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    changed = 1;
	    
	    
	  % 1000
	  elseif (cur_edge_left&&~cur_edge_right...
		   &&~test_edge_left&&~test_edge_right)&&...
		  (vertices(1,exon_idx)<=vertices(1,test_exon_idx)+5)&&...
		  (vertices(2,exon_idx)>=vertices(1,test_exon_idx))

	    vertices(2,exon_idx) = max(vertices(2,exon_idx),...
				       vertices(2,test_exon_idx));
	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	    vertices = vertices(:,new_index);
	    edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	    edges = edges(new_index,new_index);
	    
	  
	    changed = 1;
	    
	    
	    
	    
	    
	  end
%	  disp([exon_idx,test_exon_idx]);
%	  disp(edges);
%	  disp(vertices);

	  test_exon_idx = test_exon_idx+1;

	end
	exon_idx = exon_idx + 1;
      end
    end

  
  
  
  end
    
  genes(gene_idx).splicegraph = {};
  genes(gene_idx).splicegraph = {vertices,edges};
end


fprintf(1,'\n');

return




