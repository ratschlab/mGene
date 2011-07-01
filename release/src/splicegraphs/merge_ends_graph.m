


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
  
  is_intron = [];
  for exon_idx = 1:size(vertices,2)-1
    is_intron = union(is_intron,[vertices(2,exon_idx):1:vertices(1,exon_idx+1)]);
  end
  
  
  changed = 1;
  while changed
    changed = 0;
    exon_idx = 1;
    while exon_idx <= size(vertices,2)
      test_exon_idx = exon_idx+1;
      while test_exon_idx <= size(vertices,2)

	% merge the exons which are at the 5 prime end
	if (vertices(2,exon_idx)==start)&&(vertices(2,test_exon_idx)==start)
	  vertices(1,exon_idx) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  changed = 1;

	% merge the exons which are at the 3 prime end
	elseif (vertices(1,exon_idx)==stop)&&(vertices(1,test_exon_idx)==stop)
	  vertices(2,exon_idx) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
	  new_index = [1:test_exon_idx-1,test_exon_idx+1:size(vertices,2)];
	  
	  
	  vertices = vertices(:,new_index);
	  edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
	  edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
	  edges = edges(new_index,new_index);
	  changed = 1;
	  
	% merge exons 
	
	
	end
	test_exon_idx = test_exon_idx + 1;
      end
      exon_idx = exon_idx + 1;
    end
  end

  
  
  
  
    
  genes(gene_idx).splicegraph = {};
  genes(gene_idx).splicegraph = {vertices,edges};
  
  
  
end


fprintf(1,'\n');

return




