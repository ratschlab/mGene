function isoforms = find_all_isoforms(vertices,edges,strand)
% function isoforms = find_all_isoforms(vertices,edges,strand)
%
% For a DAG determined by (vertices,edges), find all the isoforms.
% Consider all transcription starts and get all paths from each start.
%
% Written by Cheng Soon Ong, 20 April 2006

global g_isoforms

% find all the leftmost and rightmost exons (transcription starts and stops)
[tstart,tstop] = find_start_stop(vertices,edges,strand);

% do a depth first search through the graph, generating a list of isoforms
% - for each transcription start, do depth first search.
g_isoforms = {};
for ixt = tstart
  dfs_visit(triu(edges),ixt,[]);
end
if strand == '-'
  for ixf = 1:length(g_isoforms)
    g_isoforms{ixf} = g_isoforms{ixf}(end:-1:1);
  end
end

isoforms = g_isoforms;
