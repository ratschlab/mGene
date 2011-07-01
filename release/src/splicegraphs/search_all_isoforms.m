function isoform = search_all_isoforms(vertices,edges,strand)
% isoform = search_all_isoforms(vertices,edges,strand)

% find all the leftmost and rightmost exons (transcription starts and stops)
tstart = [];
for ixv = 1:size(vertices,2)
  if sum(edges(ixv,1:ixv-1)) == 0
    tstart(end+1) = ixv;
  end
end

% find all transcription ends
tstop = [];
for ixv = 1:size(vertices,2)
  if sum(edges(ixv+1:end,ixv)) == 0
    tstop(end+1) = ixv;
  end
end
  
  
% do a depth first search through the graph, generating a list of isoforms
% - for each transcription start, do depth first search.
if strand == '+'
  best_isoform = [] ;
  max_len = 0 ;
  for ixt = tstart
    isoform = dfs_search_max_len(vertices, triu(edges),ixt,[]);
    len = sum(vertices(2,isoform)-vertices(1,isoform)) ;
    if len>max_len,
      max_len=len ;
      best_isoform = isoform ;
    end ;
  end
else
  best_isoform = [] ;
  max_len = 0 ;
  for ixt = tstop
    isoform = dfs_search_max_len(vertices, tril(edges),ixt,[]);
    len = sum(vertices(2,isoform)-vertices(1,isoform)) ;
    if len>max_len,
      max_len=len ;
      best_isoform = isoform ;
    end ;
  end
  isoform = best_isoform(end:-1:1);
end

