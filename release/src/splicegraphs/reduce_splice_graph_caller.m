function genes=reduce_splice_graph_caller(genes) ;
% genes=reduce_splice_graph_caller(genes) ;

for gene_idx = 1:length(genes)
%gene_idx
%genes(gene_idx)
  vertices = genes(gene_idx).splicegraph{1};
  edges = genes(gene_idx).splicegraph{2};
  %fprintf(1,'gene index:\t%d\n',gene_idx);
  if (mod(gene_idx,1000)==0)
    fprintf(1,'%d\n',gene_idx);
  end
  if isempty(vertices), continue ; end ;

  % find all the intron locations
  [dummy,exon_order] = sort(vertices(1,:),2,'ascend');
  vertices = vertices(:,exon_order);
  edges = edges(exon_order,exon_order);
  intron_loc = [];
  for ix1 = 1:size(vertices,2)-1
    for ix2 = 2:size(vertices,2)
      if edges(ix1,ix2)
	intron_loc = union(intron_loc,[vertices(2,ix1)+1:vertices(1,ix2)-1]);
      end
    end
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

	[dummy,exon_order] = sort(vertices(1,:),2,'ascend');
	if ~isequal(vertices(1,:),dummy)
	  vertices = vertices(:,exon_order);
	  edges = edges(exon_order,exon_order);
	end
	
	test_exon_idx = exon_idx + 1;
	while test_exon_idx <= size(vertices,2)
	  reduce_now = 0;
	  
          if (test_exon_idx < exon_idx) && keyboard_allowed(), keyboard; end;
	  
	  cur_edge_left = sum(edges(1:exon_idx,exon_idx));
	  test_edge_left = sum(edges(1:test_exon_idx,test_exon_idx));
	  cur_edge_right = sum(edges(exon_idx:end,exon_idx));
	  test_edge_right = sum(edges(test_exon_idx:end,test_exon_idx));
	  
	  

	  
	  % 0000
	  if (~cur_edge_left&&~cur_edge_right&&~test_edge_left&&~test_edge_right)
	    
	    if ((vertices(2,exon_idx)>=vertices(1,test_exon_idx))&&...
		(vertices(1,exon_idx)<=vertices(1,test_exon_idx)))||...
		  ((vertices(2,test_exon_idx)>=vertices(1,exon_idx))&&...
		   (vertices(1,test_exon_idx)<=vertices(1,exon_idx)))&&...
		  (sum(ismember([min(vertices(1,exon_idx),vertices(1,test_exon_idx)):...
				 max(vertices(2,exon_idx),vertices(2,test_exon_idx))],intron_loc))==0)
	      
	      vertices(1,exon_idx) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
	      vertices(2,exon_idx) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
	      new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	      vertices = vertices(:,new_index);
	      edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	      edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	      edges = edges(new_index,new_index);
	    
	  
	    
	      reduce_now = 1;
	      changed = 1;
	    end
	    
	    
	  % 0101
	  elseif (~cur_edge_left&&cur_edge_right&&~test_edge_left&&test_edge_right)

	    if (vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
		  (sum(ismember([min(vertices(1,exon_idx),vertices(1,test_exon_idx)):...
				 vertices(2,exon_idx)],intron_loc))==0)
	      
	      vertices(1,exon_idx) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
	      new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	      vertices = vertices(:,new_index);
	      edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	      edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	      edges = edges(new_index,new_index);
	    
	  
	    
	      reduce_now = 1;
	      changed = 1;
	    end

	    
	  % 1010
	  elseif (cur_edge_left&&~cur_edge_right&&test_edge_left&&~test_edge_right)
	    if (vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		  (sum(ismember([vertices(1,exon_idx):...
				 max(vertices(2,exon_idx),vertices(2,test_exon_idx))],intron_loc))==0)
	      
	      vertices(2,exon_idx) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
	      new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	      vertices = vertices(:,new_index);
	      edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	      edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	      edges = edges(new_index,new_index);
	    
	  
	    
	      reduce_now = 1;
	      changed = 1;
	    end

	  
	  % 1111
	  elseif (cur_edge_left&&cur_edge_right&&test_edge_left&&test_edge_right)
	    if (vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		 (vertices(2,exon_idx)==vertices(2,test_exon_idx))
	    
	      new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	    
	    
	      vertices = vertices(:,new_index);
	      edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	      edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	      edges = edges(new_index,new_index);
	      
	      reduce_now = 1;
	      changed = 1;
	  
	    end
	    
	  end % cases
	  

          if ~reduce_now
	    test_exon_idx = test_exon_idx+1;
	  end
	  

	end
	exon_idx = exon_idx + 1;
      end %while exon_idx <= size(vertices,2)
    end %while changed

  
  
  
  end 
    
  genes(gene_idx).splicegraph = {};
  genes(gene_idx).splicegraph = {vertices,edges};
end


fprintf(1,'\n');

return




