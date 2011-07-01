function sitegraph = unique_sitegraph(sitegraph)
% function sitegraph = unique_sitegraph(sitegraph)
%
% Merge the nodes of a sitegraph which are identical


num_nodes = size(sitegraph{2},1);
[unique_nodes,idx_orig,idx_unique] = unique(sitegraph{1}(2,:));
if length(unique_nodes) == num_nodes
  return
end

% need to check whether the nodes are of the same type (acceptor or donor)
to_merge = {};
for ix=1:length(unique_nodes)
  eq_idx = find(idx_unique==ix);
  if length(eq_idx)>1 && all(sitegraph{1}(1,eq_idx) == sitegraph{1}(1,eq_idx(1)))
    to_merge{end+1} = eq_idx;
  end
end

% merge exons
del_list = [];
for ix = 1:length(to_merge)
  del_list = [del_list,to_merge{ix}(1:end-1)];
end
new_idx = [1:num_nodes];
new_idx(del_list) = [];

sitegraph{1} = sitegraph{1}(:,new_idx);

% merge introns
for ix = 1:length(to_merge)
  cur_merge = to_merge{ix};
  rowsub = sitegraph{2}(cur_merge,:);
  idx = zeros(1,num_nodes);
  for ixr = 1:length(cur_merge)
    idx = or(idx,rowsub(ixr,:));
  end
  colidx = find(idx);
  exon_idx = colidx(find(mod(sum(sitegraph{2}(cur_merge,colidx)),3)==0));
  intron_idx = colidx(find(mod(sum(sitegraph{2}(cur_merge,colidx)),7)==0));
  assert(all(sort(union(exon_idx,intron_idx))==sort(colidx)));
  rowsub(:,exon_idx) = 3;
  rowsub(:,intron_idx) = 7;
  sitegraph{2}(cur_merge,:) = rowsub;
  sitegraph{2}(:,cur_merge) = rowsub';
end
sitegraph{2} = sitegraph{2}(new_idx,new_idx);


% sort the nodes
[dummy, node_order] = sort(sitegraph{1}(2,:),'ascend');
sitegraph{1} = sitegraph{1}(:,node_order);
sitegraph{2} = sitegraph{2}(node_order,node_order);

