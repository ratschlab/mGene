function blocks = split_regions(REGS,genes,regions_par, do_assertion)

% blocks=split_regions(REGS,genes,regions_par, do_assertion)
% 
% REGS are splitted into short blocks containing only region_par.min_genes_per_region genes if 
% split criterion holds
% split criterion: splits must have minimal distance of region_par.minCutDist (default 250)
% to both genes
% If region_par.merge_confirmed is set, only 1 gene is put into region,
% regions will be merged later with only complete and confirmed genes

  
if nargin<4
  do_assertion=0;
end
fprintf(1,'Start split_regions...\n')

if ~isfield(regions_par,'minCutDist')
  regions_par.minCutDist = 250;
end

if ~isempty(genes)
  %   load(fn_genes,'genes'); 
  map = zeros(1,length(genes));

  if regions_par.cluster_paralogs
    
    % filter out genes from regions other than the given trivial regions
    % this is necessary for correct detection of paralog blocks
    assert(isequal([genes.id],1:length(genes)))
    
    all_genes_in_REGS_idx = [];
    fprintf(1,'sorting paralogs\n')
    for r=1:length(REGS)
      fprintf(1,'check region %i\r',r)
      contig_genes_idx = find(REGS(r).chr_num==[genes.chr_num]&[genes.strand]==REGS(r).strand);
      all_genes_in_REGS_idx = [all_genes_in_REGS_idx contig_genes_idx];
    end
    for j=1:length(genes)
      fprintf(1,'check genes %i\r',j)
      other_genes_idx = find(~ismember(genes(j).paralogs,all_genes_in_REGS_idx));
      genes(j).paralogs(other_genes_idx) = [];
    end
    all_genes_in_REGS_idx = unique(all_genes_in_REGS_idx);
    genes = genes(all_genes_in_REGS_idx);
  end
else
  warning('no genes: blocks are produced without looking at gene boundaries!!!')
  % keyboard
  genes=[];
  map  =[];
end

blocks=[];
for r=1:size(REGS,2)
  REGS(r).TR_idx = r;
  [blocks,map] = split_regions_helper(blocks,REGS(r),genes,regions_par,map);
end % loop over REGS

if ~isempty(blocks)&&~isempty([blocks.gene_ids])
  assert(isequal(find(map),sort([blocks.gene_ids])))
else
  assert(~any(map))
end
cluster = {};
covered = [];
if isfield(genes,'paralogs')
  idx1 = find(~cellfun('isempty',{genes.paralogs}));
  for i = idx1
    paralogs = unique([map(genes(i).paralogs)]);
    if isempty(intersect(paralogs,covered))
      covered = [covered paralogs];
      cluster{end+1} = paralogs;
      for j= paralogs
        assert(isempty(blocks(j).paralogs))
        blocks(j).paralogs = paralogs;
      end
    else
      idx_c = [];
      covered = unique([covered paralogs]);
      for j=1:length(cluster)
        if any(ismember(paralogs,cluster{j}))
          idx_c = [idx_c,j];
        end
      end
      cluster{idx_c(1)} = [cluster{idx_c(1)},paralogs];
      for j=2:length(idx_c)
        cluster{idx_c(1)} = ([cluster{idx_c(1)},cluster{idx_c(j)}]);
        cluster{idx_c(j)} = [];
      end
      cluster{idx_c(1)} = unique(cluster{idx_c(1)});
      for j= cluster{idx_c(1)}
        blocks(j).paralogs = cluster{idx_c(1)};
      end
      
    end
  end
end
 

if do_assertion
  %check split
  fprintf(1,'check pairwise overlap of blocks\n')
  sum1 = 0;
  %ids = zeros(1,length(blocks));
  for r1 = 1:length(blocks)
    fprintf(1,'%i of %i blocks checked \r',r1,length(blocks)) 
    sum1 = sum1 + blocks(r1).stop-blocks(r1).start+1;
    ids(r1) =  blocks(r1).id;
    num_genes(r1) = length(blocks(r1).gene_ids);
   % rgenes = [];
   % for g=1:num_genes(r1)
   %   gidx = find([genes.id]==blocks(r1).gene_ids(g));
   %   assert(length(gidx)==1)
   %   rgenes = [rgenes genes(gidx)];
   % end
   % assert(max([rgenes.stop])<=blocks(r1).stop-PAR.blocks.minCutDist)
   % assert(min([rgenes.start])>=blocks(r1).start+PAR.blocks.minCutDist)
   % for r2 = (r1+1):size(blocks,2)
   %   check_overlap(blocks(r1),blocks(r2));
   % end
  end
  sum2 = 0;
  for r = 1:size(REGS,2)
    sum2 = sum2 + REGS(r).stop-REGS(r).start+1;
  end
  assert(sum1==sum2)
  assert(all(sort(ids)==1:length(blocks)))
  assert(min(num_genes)>= regions_par.min_genes_per_region)
  assert(isequal(sort([blocks.ids]),1:length(genes)))%all genes are assigned to exactly one region
  fprintf(1,'created %i blocks with min %i and max %i genes\n',length(blocks),min(num_genes), max(num_genes))
  fprintf(1,'num of genes: %ii\n',sum(num_genes))
end



if 0% ~isempty(genes)
  assert(isempty(setdiff([genes.id],[blocks.gene_ids])))
end
