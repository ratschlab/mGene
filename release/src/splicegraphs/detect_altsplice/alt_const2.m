function [idx_alt, idx_con, genes] = alt_const2(genes) ;
%function [idx_alt, idx_con, genes] = alt_const2(genes) ;

idx_alt=[] ;
idx_con=[] ;
same_5prime=0 ;
same_3prime=0 ;

for ix=1:length(genes)
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  vertices = genes(ix).splicegraph{1} ;
  edges = genes(ix).splicegraph{2} ;
  genes(ix).is_alt=0 ;
  for idx1=1:num_exons

    if idx1<num_exons
      if (sum(edges(idx1,idx1+1:end))>1)
        where=find(edges(idx1,idx1+1:end)==1) ;
        num_edges1=sum(edges(idx1,idx1:end)) ;
        for i=1:sum(edges(idx1,idx1+1:end))-1
          for j=i+1:sum(edges(idx1,idx1+1:end))
            if vertices(1,idx1+where(i))==vertices(1,idx1+where(j))
              same_3prime=same_3prime+1 ;
              num_edges1=num_edges1-1 ;
              break
            end
          end
        end
        if num_edges1>1
          genes(ix).is_alt=1 ;
          idx_alt=[idx_alt, ix] ;
          break
        end
      end 
    end
    
    if 1<idx1
      if (sum(edges(1:idx1-1,idx1))>1)
        where=find(edges(1:idx1-1,idx1)) ;
        num_edges2=sum(edges(1:idx1-1,idx1)) ; 
        for i=1:sum(edges(1:idx1-1,idx1))-1
          for j=i+1:sum(edges(1:idx1-1,idx1))
            if vertices(2,where(i))==vertices(2,where(j))
              same_5prime=same_5prime+1 ;
              num_edges2=num_edges2-1 ;
              break
            end
          end
        end
        if num_edges2>1
          genes(ix).is_alt=1 ;
          idx_alt=[idx_alt, ix] ;
          break ;
        end
      end
    end
    
  end
  if genes(ix).is_alt==0
    idx_con=[idx_con, ix] ;
  end
end

same_3prime
same_5prime

fprintf(1,'\n\nTotal alternatively spliced:\t\t\t\t\t%d\n',...
	length(idx_alt));
fprintf(1,'Total constitutively spliced:\t\t\t\t\t%d\n',...
	length(idx_con));

