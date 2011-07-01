function [transcripts,map,pos]=gentranscripts(splicegraph, id) ;
  
  edges = triu(splicegraph{2}) ;
  if nargin<2 | isempty(id)
    id = find(sum(edges)==0) ;
    if length(id)>=1,
      transcripts={} ;
      map=[] ;
      for i=1:length(id),
        [transcripts_,map_,pos]=gentranscripts(splicegraph, id(i)) ;
        transcripts(end+1:end+length(transcripts_)) = transcripts_ ;
        map=[map;map_] ;
      end ;
      return ;
    end ;
  end ;

  transcripts={} ;
  idx=find(edges(id(end),:)) ;
  if isempty(idx),
    term_id = find(sum(edges')==0) ;
    %assert(length(term_id)==1) ;
    assert(any(id(end)==term_id)) ;
    transcripts{1} = id ;
  else
    for i=idx,
      transcripts_=gentranscripts(splicegraph, [id i]) ;
      transcripts(end+1:end+length(transcripts_)) = transcripts_ ;
    end ;
  end ;
  
  if nargout>1,
    exons = splicegraph{1} ;
    pos = unique(exons(:)) ;
    start=inf ;
    stop=-inf ;
    for i=1:length(transcripts)
      start=min(start,exons(1,transcripts{i}(1)))-1 ;
      stop=max(stop,exons(2,transcripts{i}(end))) ;
    end ;
    
    exon_loc = zeros(length(transcripts),stop-start);
    for i=1:length(transcripts)
      for j=1:length(transcripts{i})
        cur_vertex = exons(:,transcripts{i}(j)) - [start;start];
        exon_loc(i,cur_vertex(1):cur_vertex(2)-1) = ...
            exon_loc(i,cur_vertex(1):cur_vertex(2)-1) + 1;
      end
    end ;
    map = exon_loc(:,pos-start) ;
  end ;
  