function sitegraph = convert_exon2ss(splicegraph,strand)
% function sitegraph = convert_exon2ss(splicegraph,strand)
%
% Convert the representation of splice graphs from exon based to splice site based.
%
% Exon based representation:
% - nodes are exons
% - edges are introns
%
% Splice site based representation (bipartite graph):
% - nodes are:
%   * transcription start sites (1)
%   * donor splice sites (3)
%   * acceptor splice sites (6)
%   * polyA sites (8), transcription stop
% - edges are:
%   * exons (3)
%   * introns (7)
%
% The reason for the strange numbering is to try to be similar to the mgene model.
%
% The bipartite graph is represented 2 matrices:
% - a 2 by m matrix, where the top row denotes the node type, 
%   and the bottom row the position in the genome.
% - an adjacency matrix with two types of edges (3 and 7)


num_exons = size(splicegraph{1},2);
num_nodes = num_exons*2;
sitegraph{1} = zeros(2,num_nodes);
sitegraph{1}(2,:) = reshape(splicegraph{1},1,num_nodes);
if strand == '+'
  sitegraph{1}(1,:) = repmat([6 3],1,num_exons);
else
  sitegraph{1}(1,:) = repmat([3 6],1,num_exons);
end

% find all the leftmost and rightmost exons (transcription starts and stops)
[tstart,tstop] = find_start_stop(splicegraph{1},splicegraph{2},strand);
if strand == '+'
  sitegraph{1}(1,tstart*2-1) = 1;
  sitegraph{1}(1,tstop*2) = 8;
else
  sitegraph{1}(1,tstart*2) = 1;
  sitegraph{1}(1,tstop*2-1) = 8;
end


sitegraph{2} = zeros(num_nodes);
% build all exon edges
for ix = 1:2:num_nodes
  sitegraph{2}(ix,ix+1) = 3;
  sitegraph{2}(ix+1,ix) = 3;
end


% build all intron edges
for ix = 1:num_exons
  for iy = ix+1:num_exons
    if splicegraph{2}(ix,iy) == 1
      sitegraph{2}(2*ix,2*iy-1) = 7;
      sitegraph{2}(2*iy-1,2*ix) = 7;
    end
  end
end



sitegraph = unique_sitegraph(sitegraph);
