function genes=merge_splice_graph(genes) ;

fprintf(1,'Number of genes:%d\n',length(genes));

for gene_idx = 1:length(genes)
  [dummy,exon_order] = sort(genes(gene_idx).splicegraph{1}(1,:),2,'ascend');
  vertices = genes(gene_idx).splicegraph{1}(:,exon_order);
  edges = genes(gene_idx).splicegraph{2}(exon_order,exon_order);
  %fprintf(1,'gene index:\t%d\n',gene_idx);
  if (mod(gene_idx,100)==0)
    fprintf(1,'.');
  end
  
%  if (size(edges,1) < 2)
%    changed = 1;
%    while changed
%      changed = 0;
%      exon_idx = 1;
%      while exon_idx <= size(vertices,2)
%	
%	test_exon_idx = exon_idx + 1;
%	while test_exon_idx <= size(vertices,2)
%	  
%	  
%	  if (vertices(1,exon_idx)<=vertices(1,test_exon_idx)) &&...
%		(vertices(2,exon_idx)>=vertices(2,test_exon_idx))
%	    vertices(1,exon_idx) = min(vertices(1,exon_idx),...
%				       vertices(1,test_exon_idx));
%	    vertices(2,exon_idx) = max(vertices(2,exon_idx),...
%				       vertices(2,test_exon_idx));
%	    new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
%	    
%	    
%	    vertices = vertices(:,new_index);
%	    %disp(vertices)
%	    
%	    %test_exon_idx = exon_idx+1;
%	    changed = 1;
%	    
%	  end
%	  test_exon_idx = test_exon_idx+1;
%	end
%	exon_idx = exon_idx + 1;
%      end
%    end
%    
%  else
%    intron_loc = [];
%    for idx1 = 1:size(vertices,2)-1
%      for idx2 = idx1:size(vertices,2)
%	if edges(idx1,idx2)
%	  intron_loc = union(intron_loc,[vertices(2,idx2):vertices(1,idx1)]);
%	end
%      end
%    end
    
    
    
    
  changed = 1;
  while changed
    changed = 0;
    exon_idx = 1;
    while exon_idx <= size(vertices,2)
      cur_exon = [vertices(1,exon_idx):vertices(2,exon_idx)];
      test_exon_idx = exon_idx + 1;
      while test_exon_idx <= size(vertices,2)
	test_exon = [vertices(1,test_exon_idx):vertices(2,test_exon_idx)];
	num_exons = size(vertices,2);
	[dummy,exon_order] = sort(vertices(1,:),2,'ascend');
	vertices = vertices(:,exon_order);
	edges = edges(exon_order,exon_order);
	
	cur_edge_left = sum(edges(1:exon_idx,exon_idx));
	test_edge_left = sum(edges(1:test_exon_idx,test_exon_idx));
	cur_edge_right = sum(edges(exon_idx:end,exon_idx));
	test_edge_right = sum(edges(test_exon_idx:end,test_exon_idx));
	  

	new_vertex = zeros(2,1);
	  
	% 0000 0001 0100 0010 1000 not considered since they have
        % exons without introns.
	  
	  
	  
	  
	  % 0011
	if ((~cur_edge_left&&~cur_edge_right...
		&&test_edge_left&&test_edge_right)&&...
	       (vertices(1,exon_idx)+5>=vertices(1,test_exon_idx))&&...
	       (vertices(2,exon_idx)-5<=vertices(2,test_exon_idx)))

	  vertices(1,exon_idx) = vertices(1,test_exon_idx);
	  vertices(2,exon_idx) = vertices(2,test_exon_idx);
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  
	  
	  
	  
	  
	  % 1100
	elseif ((cur_edge_left&&cur_edge_right...
		 &&~test_edge_left&&~test_edge_right)&&...
		(vertices(1,exon_idx)<=vertices(1,test_exon_idx)+5)&&...
		(vertices(2,exon_idx)>=vertices(2,test_exon_idx)-5))	    
	    
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  
	  
	    
	  % 0101
	elseif (vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
	      (~cur_edge_left&&cur_edge_right...
	       &&~test_edge_left&&test_edge_right)
	  
	  vertices(1,exon_idx) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));

	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  
	  

	    
	    
	  % 1010
	elseif (vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
	      (cur_edge_left&&~cur_edge_right...
	       &&test_edge_left&&~test_edge_right)
	  
	  vertices(2,exon_idx) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  
	  
	  
	    
	  % 0110 and 1001 not detected

	    
	  % 0111 
	elseif ((~cur_edge_left&&cur_edge_right...
		 &&test_edge_left&&test_edge_right)&&...
		(vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
		(vertices(1,exon_idx)+5>=vertices(1,test_exon_idx)))
	  
	  vertices(1,exon_idx) = vertices(1,test_exon_idx);
	  %vertices(2,exon_idx) = vertices(2,test_exon_idx);
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  
	  
	    
	  % 1101
	elseif ((cur_edge_left&&cur_edge_right...
		 &&~test_edge_left&&test_edge_right)&&...
		(vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
		(vertices(1,exon_idx)-5<=vertices(1,test_exon_idx)))
	  
	  %vertices(1,exon_idx) = vertices(1,test_exon_idx);
	  vertices(2,exon_idx) = vertices(2,test_exon_idx);
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  

	  % 1011  
	elseif ((cur_edge_left&&~cur_edge_right...
		 &&test_edge_left&&test_edge_right)&&...
		(vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		(vertices(2,exon_idx)-5<=vertices(2,test_exon_idx)))
		  
	  %vertices(1,exon_idx) = vertices(1,test_exon_idx);
	  vertices(2,exon_idx) = vertices(2,test_exon_idx);
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  
  
	    
	  % 1110
	elseif ((cur_edge_left&&cur_edge_right...
		 &&test_edge_left&&~test_edge_right)&&...
		(vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
		(vertices(2,exon_idx)+5>=vertices(2,test_exon_idx)))
	  
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  
	  
	  changed = 1;
	  
  
	    
	  % 1111 not detected
	end
	test_exon_idx = test_exon_idx+1;

      end
      exon_idx = exon_idx + 1;
    end
  end % while changed

  genes(gene_idx).splicegraph = {};
  genes(gene_idx).splicegraph = {vertices,edges};


end % for each gene_idx




fprintf(1,'\n');

return




