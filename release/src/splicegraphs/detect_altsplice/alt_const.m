function genes= alt_const(genes) ;
%function genes = alt_const(genes) ;

tot_exons = 0;
for ix=1:length(genes)
  if (mod(ix,1000)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  tot_exons = tot_exons + num_exons ;
  vertices = genes(ix).splicegraph{1} ;
  edges = genes(ix).splicegraph{2} ;

  if isempty(edges)
    genes(ix).is_alt_spliced=0 ;
    genes(ix).is_alt=0 ;
    continue ;
  end 
  
  init_alt_idx = [] ;
  term_alt_idx = [] ;
  for i=1:num_exons,
    if ~any(edges(1:i-1,i)) % initial
      idx=setdiff(find(vertices(2,i)==vertices(2,:)),i) ;
      if ~isempty(idx) % other exon with same end
        is_simple=1 ; 
        for j=idx,
          if ~isequal(find(edges(j,j+1:end))+j,find(edges(i,i+1:end))+i)
            is_simple = 0 ;
            break ;
          end ;
          idx_prev = find(edges(1:j-1,j)) ;
          for k=idx_prev,
            if vertices(2,k)>vertices(1,i)
              is_simple = 0 ;
              break ;
            end ;
          end ;
        end ;
        if ~is_simple, continue ; end ;
        init_alt_idx(end+1) = i ;
      end
    end ;
    if ~any(edges(i,i+1:num_exons)) % terminal
      idx=setdiff(find(vertices(1,i)==vertices(1,:)),i) ;
      if ~isempty(idx) % other exon with same start
        is_simple=1 ; 
        for j=idx,
          if ~isequal(find(edges(1:j-1,j)),find(edges(1:i-1,i)))
            is_simple = 0 ;
            break ;
          end ;
          idx_next = find(edges(j,j+1:end))+j ;
          for k=idx_next,
            if vertices(1,k)<vertices(2,i)
              is_simple = 0 ;
              break ;
            end ;
          end ;
        end ;
        if ~is_simple, continue ; end ;
        term_alt_idx(end+1) = i ;
      end ;
    end ;
  end ;
  take_idx = setdiff(1:num_exons, [init_alt_idx,term_alt_idx]) ;
  vertices = genes(ix).splicegraph{1}(:,take_idx) ;
  edges = genes(ix).splicegraph{2}(take_idx,take_idx) ;

  start=min(vertices(1,:)-50) ;
  stop=max(vertices(2,:)+50) ;
  
  exon_loc = zeros(1,stop-start);

  for i=1:size(edges,1),
    for j=i+1:size(edges,1),
      if edges(i,j)==1,
        cur_edge = [vertices(2,i)+1-start, vertices(1,j)-1-start] ;
        exon_loc(cur_edge(1):cur_edge(2)) = ...
            exon_loc(cur_edge(1):cur_edge(2)) + 1;
      end ;
    end ;
  end ;

  for i=1:size(vertices,2)
    cur_vertex = vertices(:,i) - [start;start];
    exon_loc(cur_vertex(1):cur_vertex(2)) = ...
        exon_loc(cur_vertex(1):cur_vertex(2)) + 1;
  end ;

  if max(exon_loc)>1
    genes(ix).is_alt_spliced=1 ;
    genes(ix).is_alt=1 ;
  else
    genes(ix).is_alt_spliced=0 ;
    if isempty([init_alt_idx,term_alt_idx])
      genes(ix).is_alt=0 ;
    else
      genes(ix).is_alt=1 ;
    end ;
  end
end

fprintf(1,'\n\nTotal genes:\t\t\t\t\t\t\t%d\n',...
        length(genes));
fprintf(1,'Total exons:\t\t\t\t\t\t\t%d\n',...
        tot_exons);
fprintf(1,'Total genes with alternative isoforms:\t\t\t\t%d\n',...
        sum([genes(:).is_alt]));
fprintf(1,'Total genes alternatively spliced:\t\t\t\t%d\n',...
        sum([genes(:).is_alt_spliced]));
fprintf(1,'Total constitutively spliced:\t\t\t\t\t%d\n',...
        sum(~[genes(:).is_alt_spliced]));

