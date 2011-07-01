
disp('Building Compact Splice Graphs...')

%for gene_idx = 1: length(genes)
for gene_idx = 18

  % construct a new empty splice graph, based on
  % http://proline.bic.nus.edu.sg/dedb/
  vertices =  [];
  edges = [];
  est_ends = [];

  for transcript_idx = 1:length(genes(gene_idx).transcripts)
    exon_start_end = genes(gene_idx).exons{transcript_idx};

    for exon_idx = 1:(size(exon_start_end,1)-1)
      if (genes(gene_idx).strands(transcript_idx) == '+')
	exon1_start = exon_start_end(exon_idx,1);
	exon1_end = exon_start_end(exon_idx,2);
	exon2_start = exon_start_end(exon_idx+1,1);
	exon2_end = exon_start_end(exon_idx+1,2);
      else
	exon1_start = exon_start_end(exon_idx,2);
	exon1_end = exon_start_end(exon_idx,1);
	exon2_start = exon_start_end(exon_idx+1,2);
	exon2_end = exon_start_end(exon_idx+1,1);
      end      
      
      match1 = 0;
      match2 = 0;
      if isempty(vertices)
	vertices(1,1) = exon1_start;
	vertices(2,1) = exon1_end;
	vertices(1,2) = exon2_start;
	vertices(2,2) = exon2_end;
	edges = zeros(2);
	edges(1,2) = 1;
	edges(2,1) = 1;
	num_exons = 2;
	est_ends = ['s'];
	if (exon_idx == size(exon_start_end,1)-1)
	  est_ends(2) = 'e';
	else
	  est_ends(2) = 'm';
	end
      else
	exon1_idx = 0;
	exon2_idx = 0;
	idx = 1;
	while (idx <= num_exons)
	  match1 = exon_match([exon1_start;exon1_end],vertices(:,idx),...
			      exon_start_end,est_ends(idx));
	  if (match1~=0)
	    exon1_idx = idx;
	    break;
	  end
	  
	  idx = idx+1;
	end

	idx = 1;
	while (idx <= num_exons)
	  match2 = exon_match([exon2_start;exon2_end],vertices(:,idx),...
			      exon_start_end,est_ends(idx));
	  if (match2~=0)
	    exon2_idx = idx;
	    break;
	  end
	  
	  idx = idx+1;
	end

	
	
	if (match1~=0) && (match2~=0)
	  edges(exon1_idx,exon2_idx) = 1;
	  edges(exon2_idx,exon1_idx) = 1;
	  if (match1==2)
	    vertices(1,exon1_idx) = min([vertices(1,exon1_idx), ...
		    exon1_start]);
	    if (exon_idx~=1) && (est_ends(exon1_idx)=='s')
	      est_ends(exon1_idx)='m';
	    end
	  end
	  if (match2==3)
	    vertices(2,exon2_idx) = max([vertices(2,exon2_idx), ...
		    exon2_end]);
	    if (exon_idx~=size(exon_start_end,1)-1)...
		  && (est_ends(exon1_idx)=='e')
	      est_ends(exon1_idx)='m';
	    end
	  end
	  
	else
	  if (match1==0) && (match2~=0)
	    vertices(1,num_exons+1) = exon1_start;
	    vertices(2,num_exons+1) = exon1_end;
	    if (exon_idx==1)
	      est_ends(num_exons+1)='s';
	    else
	      est_ends(num_exons+1)='m';
	    end
	    edges(exon2_idx,num_exons+1) = 1;
	    edges(num_exons+1,exon2_idx) = 1;
	    num_exons = num_exons + 1;
	  elseif (match1~=0) && (match2==0)
	    vertices(1,num_exons+1) = exon2_start;
	    vertices(2,num_exons+1) = exon2_end;
	    if (exon_idx==size(exon_start_end,1)-1)
	      est_ends(num_exons+1)='e';
	    else
	      est_ends(num_exons+1)='m';
	    end
	    edges(exon1_idx,num_exons+1) = 1;
	    edges(num_exons+1,exon1_idx) = 1;
	    num_exons = num_exons + 1;
	  else
	    assert((exon1_idx==0)&&(exon2_idx==0));
	    vertices(1,num_exons+1) = exon1_start;
	    vertices(2,num_exons+1) = exon1_end;
	    if (exon_idx==1)
	      est_ends(num_exons+1)='s';
	    else
	      est_ends(num_exons+1)='m';
	    end
	    num_exons = num_exons + 1;

	    vertices(1,num_exons+1) = exon2_start;
	    vertices(2,num_exons+1) = exon2_end;
	    if (exon_idx==size(exon_start_end,1)-1)
	      est_ends(num_exons+1)='e';
	    else
	      est_ends(num_exons+1)='m';
	    end
	    num_exons = num_exons + 1;	  
	    
	    edges(num_exons-1,num_exons) = 1;
	    edges(num_exons,num_exons-1) = 1;
	  end
	end
	
	
	
      end
      genes(gene_idx).exons{transcript_idx}
      disp([match1,match2])
      est_ends
    end
  end
  genes(gene_idx).splicegraph = {vertices,edges};
%  viewsplicegraph(genes(gene_idx).exons,genes(gene_idx).name, ...
%		  genes(gene_idx).transcripts, [], ...
%		  genes(gene_idx).strands, ...
%		  genes(gene_idx).splicegraph{1}, ...
%		  genes(gene_idx).splicegraph{2})
%  keyboard

end



return




