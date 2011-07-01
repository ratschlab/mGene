function genes = build_splice_graph_caller(genes) ;
% genes = build_splice_graph_caller(genes) ;
  
for gene_idx = 1:length(genes)
  % construct a new empty splice graph, based on
  % http://proline.bic.nus.edu.sg/dedb/
  vertices =  [];
  edges = [];
%  fprintf(1,'gene index:\t%d\n',gene_idx);
  if (mod(gene_idx,100)==0)
    fprintf(1,'.');
  end

  for transcript_idx = 1:length(genes(gene_idx).transcripts)
    exon_start_end = genes(gene_idx).exons{transcript_idx};
    
    if (size(exon_start_end,1) == 1)
      exon1_start = exon_start_end(1,1);
      exon1_end = exon_start_end(1,2);

      
      if isempty(vertices)
	vertices(1,1) = exon1_start;
	vertices(2,1) = exon1_end;
	edges = 0;
	num_exons = 1;
      else
	vertices(1,num_exons+1) = exon1_start;
	vertices(2,num_exons+1) = exon1_end;
	edges(1,num_exons+1) = 0;
	edges(num_exons+1,1) = 0;
	num_exons = num_exons + 1;
      end

      
    else
      

      for exon_idx = 1:(size(exon_start_end,1)-1)
	exon1_start = exon_start_end(exon_idx,1);
	exon1_end = exon_start_end(exon_idx,2);
	exon2_start = exon_start_end(exon_idx+1,1);
	exon2_end = exon_start_end(exon_idx+1,2);


	
	if isempty(vertices)
	  vertices(1,1) = exon1_start;
	  vertices(2,1) = exon1_end;
	  vertices(1,2) = exon2_start;
	  vertices(2,2) = exon2_end;
	  edges = zeros(2);
	  edges(1,2) = 1;
	  edges(2,1) = 1;
	  num_exons = 2;
	else
	  exon1_idx = 0;
	  exon2_idx = 0;
	  idx = 1;
	  while idx <= num_exons
%	    if ((vertices(1,idx)==exon1_start)||...
%		(vertices(1,idx)+1==exon1_start)||...
%		(vertices(1,idx)-1==exon1_start))...
%		  && ((vertices(2,idx)==exon1_end)||...
%		      (vertices(2,idx)+1==exon1_end)||...
%		      (vertices(2,idx)-1==exon1_end))
	    if ((vertices(1,idx)==exon1_start) &&...
		(vertices(2,idx)==exon1_end))
	      exon1_idx = idx;
	    end
	    idx = idx+1;
	  end
	  idx = 1;
	  while idx <= num_exons
%	    if ((vertices(1,idx)==exon2_start)||...
%		(vertices(1,idx)+1==exon2_start)||...
%		(vertices(1,idx)-1==exon2_start))...
%		  && ((vertices(2,idx)==exon2_end)||...
%		      (vertices(2,idx)+1==exon2_end)||...
%		      (vertices(2,idx)-1==exon2_end))
	    if ((vertices(1,idx)==exon2_start) &&... 
		(vertices(2,idx)==exon2_end))
	      exon2_idx = idx;
	    end
	    idx = idx+1;
	  end
	  if (exon1_idx~=0) && (exon2_idx~=0)
	    edges(exon1_idx,exon2_idx) = 1;
	    edges(exon2_idx,exon1_idx) = 1;
	  else
	    if ((exon1_idx==0) && (exon2_idx~=0))
	      vertices(1,num_exons+1) = exon1_start;
	      vertices(2,num_exons+1) = exon1_end;
	      edges(exon2_idx,num_exons+1) = 1;
	      edges(num_exons+1,exon2_idx) = 1;
	      num_exons = num_exons + 1;
	    elseif ((exon2_idx==0) && (exon1_idx~=0))
	      vertices(1,num_exons+1) = exon2_start;
	      vertices(2,num_exons+1) = exon2_end;
	      edges(exon1_idx,num_exons+1) = 1;
	      edges(num_exons+1,exon1_idx) = 1;
	      num_exons = num_exons + 1;
	    else
	      assert((exon1_idx==0)&&(exon2_idx==0));
	      vertices(1,num_exons+1) = exon1_start;
	      vertices(2,num_exons+1) = exon1_end;
	      num_exons = num_exons + 1;
	      vertices(1,num_exons+1) = exon2_start;
	      vertices(2,num_exons+1) = exon2_end;
	      num_exons = num_exons + 1;	  
	      
	      edges(num_exons-1,num_exons) = 1;
	      edges(num_exons,num_exons-1) = 1;
	    end
	  end
	end
      end
    end
  end




  genes(gene_idx).splicegraph = {vertices,edges};
end


fprintf(1,'\n');


return





