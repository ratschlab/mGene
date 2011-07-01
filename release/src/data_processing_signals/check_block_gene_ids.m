function check_block_gene_ids(blocks, genes)
% find all genes that overlap with a block and check if the gene_ids list of that block is correct and complete
% 
% this is also a check for the blocks being produced correctly without cutting any gene.


chrs = unique([blocks.chr_num]);
strands = unique([blocks.strand]);
wrong_idx = [];
cut_idx = [];

block_chr_num = [blocks.chr_num] ;
genes_chr_num = [genes.chr_num] ;
block_strand  = [blocks.strand] ;
genes_strand  = [genes.strand] ;

for c = chrs
  %cblocks = blocks(block_chr_num==c);
  %cgenes  = genes(genes_chr_num==c); 
  for s=strands 
    %sblocks = cblocks([cblocks.strand]==s);
    %sgenes  = cgenes([cgenes.strand]==s); 
    sblocks = blocks(block_chr_num==c & block_strand==s) ;
    sgenes = genes(genes_chr_num==c & genes_strand==s) ;
    for j=1:length(sblocks),
      if mod(j,100)==0, fprintf('block: %i (%i)\r',j,length(sblocks)), end 
      idx = find([sgenes.stop]>sblocks(j).start & [sgenes.start]<sblocks(j).stop);
      %assert(isequal([sgenes(idx).id],sblocks(j).gene_ids))
      if ~isequal(sort([sgenes(idx).id]),sort(sblocks(j).gene_ids)),
        if keyboard_allowed(),
          keyboard
        end ;
        wrong_idx = [wrong_idx sblocks(j).id]; 
        g = sgenes(idx);
        for k = 1:length(g)
          if g(k).stop>sblocks(j).stop|g(k).start<sblocks(j).start
            cut_idx = [cut_idx g(k).id];
          end
        end
        %keyboard
      end, 
    end%sblocks
  end%strand
end%chr
if length(wrong_idx)>0
  fprintf('%i blocks have incorrect lists of associated genes\n',length(wrong_idx))
end
if length(cut_idx)>0
  fprintf('%i genes have been cut by a block boundary\n',length(cut_idx))
end