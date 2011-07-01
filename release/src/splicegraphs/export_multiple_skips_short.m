
fsock = fopen('multiple_exons_skipped.txt','w');
%fsock=1;




fprintf(fsock,['The following exons were alternatively spliced in' ...
	       ' with multiple skips.\n\n']);


for ix=unique(idx_multiple_skips)
  
  num_exons = size(genes(ix).splicegraph{1},2);
  edges = genes(ix).splicegraph{2};
  labels = repmat([1:num_exons]',1,num_exons);
  for exon_idx_first = 1:num_exons-3
    for exon_idx_last = (exon_idx_first+3):num_exons
      if (edges(exon_idx_first,exon_idx_last))
	
	
	
	% find all pairs shortest path
	exist_path = triu(edges);
	exist_path(exon_idx_first,exon_idx_last)=0;
	exist_path(exist_path==0)=inf ;
	for i=1:num_exons; exist_path(i,i)=0 ; end ;

	
	path = ~isinf(exist_path).*labels;
	
	for k=1:num_exons,
	  for i=1:num_exons,
	    idx=find(exist_path(i,k)+exist_path(k,:)<exist_path(i,:));
	    exist_path(i,idx)=exist_path(i,k)+exist_path(k,idx) ;
	    path(i,idx) = path(k,idx) ;
	  end 
	end 
	
	if exist_path(exon_idx_first,exon_idx_last)>2 & ~isinf(exist_path(exon_idx_first,exon_idx_last))
	  backtrace = path(exon_idx_first,exon_idx_last);
	  while (backtrace(end) > exon_idx_first)
	    backtrace = [backtrace,path(exon_idx_first,backtrace(end))];
	  end
	  backtrace = backtrace(1:end-1);
	  backtrace = backtrace(length(backtrace):-1:1);
	  if length(backtrace)<2,keyboard;end


	  fprintf(fsock,'Gene index (for debugging only):%d\n',ix);
	  fprintf(fsock,'Gene name:%s\n',genes(ix).name);
	  fprintf(fsock,'Chromosome:%s\n',genes(ix).chr);
	  fprintf(fsock,'Strand:%c\n',genes(ix).strands(1));
	  

  
	  
	  
	  fprintf(fsock,['Exons that are consecutively skipped, first' ...
			 ' and last are the ends\n']);
	  vertices = genes(ix).splicegraph{1};
	  fprintf(fsock,'(%d,%d) -- ',vertices(:,exon_idx_first));
	  
	  for skipcount = 1:length(backtrace)
	    fprintf(fsock,'(%d,%d) ',vertices(:,backtrace(skipcount))');
	  end
	  fprintf(fsock,'-- (%d,%d)\n',vertices(:,exon_idx_last));
	  
	  idx_multiple_skips = [idx_multiple_skips,repmat(ix,1,length(backtrace))];
	  exon_multiple_skips = [exon_multiple_skips,backtrace];
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
