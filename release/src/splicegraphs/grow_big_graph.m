


for gene_idx = 1:length(genes)
  vertices = genes(gene_idx).splicegraph{1};
  edges = genes(gene_idx).splicegraph{2};
  %fprintf(1,'gene index:\t%d\n',gene_idx);
  if (mod(gene_idx,100)==0)
    fprintf(1,'.');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Some simple merging to reduce exons
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  start = min(vertices(2,:));
  stop = max(vertices(1,:));
  changed = 1;
  while changed
    changed = 0;
    exon_idx = 1;
    while exon_idx <= size(vertices,2)
      test_exon_idx = exon_idx+1;
      while test_exon_idx <= size(vertices,2)
	
	if (vertices(2,exon_idx)==start)&&(vertices(2,test_exon_idx)==start)
	  vertices(1,exon_idx) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  changed = 1;
	  
	elseif (vertices(1,exon_idx)==stop)&&(vertices(1,test_exon_idx)==stop)
	  vertices(2,exon_idx) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  changed = 1;
	  
	end
	test_exon_idx = test_exon_idx + 1;
      end
      exon_idx = exon_idx + 1;
    end
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Merge exons and keep the old ones too.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  changed = 1;
  while changed
    changed = 0;
    exon_idx = 1;
    while exon_idx <= size(vertices,2)
      
      test_exon_idx = exon_idx + 1;
      while test_exon_idx <= size(vertices,2)
	
	num_exons = size(vertices,2);
	
	
	cur_edge_left = 0;
	idx = 1;
	while (cur_edge_left == 0) && (idx<=num_exons)
	  if ((edges(exon_idx,idx)) &&...
	      (vertices(2,idx) < vertices(1,exon_idx)))
	    cur_edge_left = 1;
	  end
	  idx = idx + 1;
	end
	
  

	test_edge_left = 0;
	idx = 1;
	while (test_edge_left == 0) && (idx<=num_exons)
	  if ((edges(test_exon_idx,idx)) &&...
	      (vertices(2,idx) < vertices(1,test_exon_idx)))
	    test_edge_left = 1;
	  end
	  idx = idx + 1;
	end



	cur_edge_right = 0;
	idx = 1;
	while (cur_edge_right == 0) && (idx<=num_exons)
	  if ((edges(exon_idx,idx)) &&...
	      (vertices(1,idx) > vertices(2,exon_idx)))
	    cur_edge_right = 1;
	  end
	  idx = idx + 1;
	end
	

	test_edge_right = 0;
	idx = 1;
	while (test_edge_right == 0) && (idx<=num_exons)
	  if ((edges(test_exon_idx,idx)) &&...
	      (vertices(1,idx) > vertices(2,test_exon_idx)))
	    test_edge_right = 1;
	  end
	  idx = idx + 1;
	end

	  


	  

	new_vertex = zeros(2,1);
	  
	% 0000
	if ((~cur_edge_left&&~cur_edge_right&&...
		~test_edge_left&&~test_edge_right)&&...
	       ((vertices(1,exon_idx)~=vertices(1,test_exon_idx))||...
		(vertices(2,exon_idx)~=vertices(2,test_exon_idx))))
	  

	  new_vertex(1) = min(vertices(1,exon_idx),...
			      vertices(1,test_exon_idx));
	  new_vertex(2) = max(vertices(2,exon_idx),...
			      vertices(2,test_exon_idx));


	% 0101
	elseif ((~cur_edge_left&&cur_edge_right...
		 &&~test_edge_left&&test_edge_right)&&...
		(vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
		(vertices(1,exon_idx)~=vertices(1,test_exon_idx)))

	  new_vertex(1) = min(vertices(1,exon_idx),...
			      vertices(1,test_exon_idx));
	  new_vertex(2) = vertices(2,exon_idx);
	  

	    

	% 1010
	elseif ((cur_edge_left&&~cur_edge_right...
		 &&test_edge_left&&~test_edge_right)&&...
		(vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		(vertices(2,exon_idx)~=vertices(2,test_exon_idx)))
	  
	  new_vertex(1) = vertices(1,exon_idx);
	  new_vertex(2) = max(vertices(2,exon_idx),...
			      vertices(2,test_exon_idx));
	  
	  

	    

	% 0001
	elseif (~cur_edge_left&&~cur_edge_right...
		&&~test_edge_left&&test_edge_right)&&...
	      (vertices(2,exon_idx)-5<=vertices(2,test_exon_idx))&&...
	      (vertices(2,exon_idx)>=vertices(1,test_exon_idx))
	  
	  new_vertex(1) = min(vertices(1,exon_idx),...
				     vertices(1,test_exon_idx));
	  new_vertex(2) = vertices(2,test_exon_idx);
	  
	    
	% 0010
	elseif (~cur_edge_left&&~cur_edge_right...
		&&test_edge_left&&~test_edge_right)&&...
	      (vertices(1,exon_idx)+5>=vertices(1,test_exon_idx))&&...
	      (vertices(1,exon_idx)<=vertices(2,test_exon_idx))
	  
	  new_vertex(1) = vertices(1,test_exon_idx);
	  new_vertex(2) = max(vertices(2,exon_idx),...
			      vertices(2,test_exon_idx));
	  
	    
	% 0100
	elseif (~cur_edge_left&&cur_edge_right...
		&&~test_edge_left&&~test_edge_right)&&...
	      (vertices(2,exon_idx)>=vertices(2,test_exon_idx)-5)&&...
	      (vertices(1,exon_idx)<=vertices(2,test_exon_idx))
	  
	  new_vertex(1) = min(vertices(1,exon_idx),...
			      vertices(1,test_exon_idx));
	  new_vertex(2) = vertices(2,exon_idx);
	  
	    
	% 1000
	elseif (cur_edge_left&&~cur_edge_right...
		&&~test_edge_left&&~test_edge_right)&&...
	      (vertices(1,exon_idx)<=vertices(1,test_exon_idx)+5)&&...
	      (vertices(2,exon_idx)>=vertices(1,test_exon_idx))
	  
	  new_vertex(1) = vertices(1,exon_idx);
	  new_vertex(2) = max(vertices(2,exon_idx),...
			      vertices(2,test_exon_idx));
	    
  
	    
	end

	if (new_vertex(1)~=0)||(new_vertex(2)~=0)
	
	  % Include the new exon	  
	  search_idx = 1;
	  known_vertex = 0;
	  known_edges = 0;
	  while (search_idx<=num_exons) && ~known_vertex
	    if (new_vertex(1)==vertices(1,search_idx))&&...
		  (new_vertex(2)==vertices(2,search_idx))
	      known_vertex = 1;
	    end
	    new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
	    if (equal(new_edges,edges(exon_idx,:))||equal(new_edges,edges(exon_idx,:)))
	      known_edges = 1;
	    end
	    if known_vertex && known_edges, break; end;
	    search_idx = search_idx+1;
	  end
	  
		  
	  if ~known_vertex && ~known_edges
	    vertices(:,num_exons+1) = new_vertex;
	    
	  
	    edges(num_exons+1,:) = or(edges(exon_idx,:),...
				      edges(test_exon_idx,:));
	    edges(:,num_exons+1) = or(edges(:,exon_idx),...
				      edges(:,test_exon_idx));
	    
	    edges(exon_idx,:) = edges(num_exons+1,:);
	    edges(:,exon_idx) = edges(:,num_exons+1);
	    edges(test_exon_idx,:) = edges(num_exons+1,:);
	    edges(:,test_exon_idx) = edges(:,num_exons+1);
	    
	    changed = 1;
	  elseif known_vertex && ~known_edges
	    edges(exon_idx,:) = or(edges(exon_idx,:),...
				      edges(test_exon_idx,:));
	    edges(:,exon_idx) = or(edges(:,exon_idx),...
				      edges(:,test_exon_idx));
	    edges(test_exon_idx,:) = edges(exon_idx,:);
	    edges(:,test_exon_idx) = edges(:,exon_idx);

	    changed = 1;
	  end
	end

	
	
	test_exon_idx = test_exon_idx+1;
	
      end
      exon_idx = exon_idx + 1;
    end
  end
  
  
  
    
  genes(gene_idx).splicegraph = {};
  genes(gene_idx).splicegraph = {vertices,edges};
  
  
  
end


fprintf(1,'\n');

return




