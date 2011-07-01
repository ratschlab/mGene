function [regions,map] = split_regions_helper(regions,REG,GENES,regions_par,map)
% [regions,map] = split_regions_helper(regions,REG,GENES,regions_par,map)
  
new_tr=1;
g=0;
contigstart = REG.start;
contigstop = REG.stop;

if regions_par.merge_confirmed
  min_genes_per_region = 1;
  fpritnf('putting only 1 gene in region; regions will be merged later with only complete genes')
else 
  min_genes_per_region = regions_par.min_genes_per_region;
end


if 0%~isempty(GENES)
  GENES_chr_num=[GENES.chr_num] ;
  GENES_strand=[GENES.strand] ;
  GENES_start=[GENES.start]+5 ; % hack for some minor problems
  GENES_stop=[GENES.stop]-5 ;   % hack for some minor problems

  idx =  find(REG.chr_num==GENES_chr_num& GENES_strand==REG.strand & ...
              REG.start <= GENES_start & REG.stop >= GENES_stop);
  idx_overlap = find(REG.chr_num==GENES_chr_num & GENES_strand==REG.strand & ...
                     REG.start > GENES_start & REG.start < GENES_stop);
  assert(isempty(idx_overlap))
  idx_overlap = find(REG.chr_num==GENES_chr_num & GENES_strand==REG.strand & ...
                     REG.stop > GENES_start & REG.stop < GENES_stop);
  assert(isempty(idx_overlap))
else 
  idx=[];
end

fprintf(1,'split region %s%s, num of genes found: %i  \r',REG.chr,REG.strand, length(idx)) 
if ~isempty(idx)
  genes = GENES(idx);
  [tmp pos] = sort([genes.start]);
  genes = genes(pos);
  clear pos idx start_pos
else
  block_length = REG.stop-REG.start;
  if ~isfield(regions_par,'max_length')|isempty(regions_par.max_length)
    num_blocks = 1;
  else
    num_blocks = ceil(block_length/regions_par.max_length);
    if num_blocks>0,
      block_length = ceil(block_length/num_blocks);
    else
      block_length = 0 ;
    end ;
  end
  start = REG.start+2*regions_par.minCutDist;
  genes = [];
  randn('seed', 1234);
  %if keyboard_allowed(), keyboard ; end ;
  for j=1:num_blocks
    dist = round(abs(randn(1,2)*regions_par.minCutDist/2+2*regions_par.minCutDist));
    stop = min(start+block_length, REG.stop);
    genes(j).id = [];
    genes(j).transcript_status =[];
    genes(j).is_complete = [];
    genes(j).start = start+dist(1);
    genes(j).stop = stop-dist(2);
    start = stop;
  end
end

%genes
cnt=0 ;
while g<size(genes,2)
  cnt = cnt + 1 ;
  if cnt>g+20, fprintf('break') ; break ; end ;

  region = [];
  region.organism = REG.organism;
  region.chr = REG.chr;
  region.chr_num  = REG.chr_num;
  region.strand  = REG.strand;
  region.id = length(regions)+1;
  region.TR_idx = REG.TR_idx ;
  region.num_fragments = 1;
  region.config = REG.config;
  region.paralogs = [];
  r_genes = [];
  region.gene_ids=[];
  region.gene_status=[];
  region.gene_complete=[];
  genedist_too_short=0;
  if new_tr
    region.start = REG.start;
    region.stop = REG.start;
    new_tr=0;
  else
    region.start = region_end_stop;
  end
  %add genes to region as long as min gene num is not achieved, genedist is too short or more if you are at 
  %the end of the contig and last region would be to short
  while  g<size(genes,2)&&...
        (length(r_genes)<min_genes_per_region||...
         g+min_genes_per_region>size(genes,2)||...
         genedist_too_short ||...
         (length(r_genes)>0 && isfield(region, 'stop') && ...
          (region.stop+regions_par.offset>contigstop||region.stop-regions_par.offset<1))) 
    % is this the last gene of the contig?
    % enough genes for this region collected?
    % still enough genes to fill another region?
    % next gene closer than 2*minCutDist?
    g=g+1;
    r_genes = [r_genes genes(g)];
    region.gene_ids = [region.gene_ids genes(g).id];
	if isfield(genes(g), 'transcript_status')
	    region.gene_status = [region.gene_status max(genes(g).transcript_status)];
	end
    region.gene_complete = [region.gene_complete max(genes(g).is_complete)];
    map(genes(g).id) = region.id;
    gene_stop = max([r_genes.stop]);
    if g==size(genes,2)
      region.stop = REG.stop;
    else
      genedist = genes(g+1).start-gene_stop-1;
      rand_num = ceil(rand*(genedist-2*regions_par.minCutDist));
      region.stop = gene_stop+regions_par.minCutDist+rand_num;
      genedist_too_short = (genedist<=2*regions_par.minCutDist+1);
    end
  end%while add genes
  region.left_dist = min([r_genes.start])- region.start;
  region.right_dist = region.stop - gene_stop;
  %   assert(length(region.gene_ids)>=min_genes_per_region)
  region_end_stop = region.stop+1;
  
  %if ~isempty(GENES)
   % check_block_gene_ids(region, GENES);
  %end

  region.reg_no_overlap(1) = max([contigstart+regions_par.offset region.start]);
  region.reg_no_overlap(2) = min([contigstop-regions_par.offset region.stop]);
  region.start = max([contigstart region.start-regions_par.offset]);
  region.stop = min([contigstop region.stop+regions_par.offset]);
  region.reg_no_overlap = region.reg_no_overlap-region.start+1;
 
  switch region.strand
   case '+',
    region.offset = region.start - 1;
   case '-',
    region.offset = region.stop + 1;
  end
  % keyboard
  regions = [regions region];
end % loop over genes in contig

