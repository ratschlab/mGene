


fsock = fopen('xor_exons.txt','w');
%fsock=1;



fprintf(fsock,['The following exons were alternatively spliced in' ...
	       ' an exclusive fashion.\n\n']);

for ix=unique(idx_xor_exons)
  fprintf(fsock,'Gene index (for debugging only):%d\n',ix);
  fprintf(fsock,'Gene name:%s\n',genes(ix).name);
  fprintf(fsock,'Chromosome:%s\n',genes(ix).chr);
  fprintf(fsock,'Strand:%c\n',genes(ix).strands(1));
  
  
  fprintf(fsock,'Exons exclusively spliced\n');
  num_exons = size(genes(ix).splicegraph{1},2);
  for exon_idx = 1:num_exons-2
    for exon_idx1 = (exon_idx+1):(num_exons-1)
      for exon_idx2 = (exon_idx1+1):num_exons
	for exon_idx3 = (exon_idx2+1):num_exons
	  edges = genes(ix).splicegraph{2};
	  vertices = genes(ix).splicegraph{1};
	  if (edges(exon_idx1,exon_idx3))&&(edges(exon_idx,exon_idx2))&&...
		(edges(exon_idx,exon_idx1))&&(~edges(exon_idx1,exon_idx2))&&...
		(edges(exon_idx2,exon_idx3))&&...
		(vertices(1,exon_idx2)>vertices(2,exon_idx1))
	    %fprintf(1,'%d [%d xor %d] %d\n',...
	    %	    exon_idx,exon_idx1,exon_idx2,exon_idx3);
	    fprintf(fsock,'(%d,%d) [(%d,%d) xor (%d,%d)] (%d,%d)\n',...
		    vertices(1,exon_idx),vertices(2,exon_idx),...
		    vertices(1,exon_idx1),vertices(2,exon_idx1),...
		    vertices(1,exon_idx2),vertices(2,exon_idx2),...
		    vertices(1,exon_idx3),vertices(2,exon_idx3));
	  end
	end
      end
    end
  end

  
%  for ix2 =1:length(genes(ix).transcripts)
%    fprintf(fsock,'Transcript %d, strand %c :%s\n',ix2, ...
%	    genes(ix).strands(ix2),genes(ix).transcripts{ix2});
%  end
  
  
  fprintf(fsock,'\n\n');
  
  
end

fclose(fsock);
