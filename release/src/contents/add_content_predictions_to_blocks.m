function blocks = add_content_predictions_to_blocks(blocks, cont_pred_dirs)
% blocks = add_content_predictions_to_blocks(blocks, cont_pred_dirs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if binary files exist, use the faster function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fexist(sprintf('%s/pred/contig_1+.pos',cont_pred_dirs{1})),
  fprintf('using add_bin_content_prediction()...\n');
  blocks = add_bin_content_predictions_to_blocks(blocks, cont_pred_dirs);
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% else we end up here using the slower function for mat files:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contigs = unique({blocks.chr}) ;

for i=1:length(cont_pred_dirs),
  fnames{i}=  sprintf('%s/output_spf.mat', cont_pred_dirs{i}) ;
  L = load(fnames{i}, 'SPF_info') ;
  SPF_info{i} = L.SPF_info ;
end ;

fprintf('Processing %i contigs: \n', length(contigs)) ;
for j=1:length(contigs)
  if mod(j, 5)==0 || length(contigs)<10, 
    fprintf('%i ... ', j) ; 
  end ;

  contig_block_ids = find(ismember({blocks.chr}, contigs{j}));
  strands = unique({blocks.strand});
  for k=1:length(strands)
    strand_block_ids = intersect(strmatch(strands(k), {blocks.strand}, 'exact'), contig_block_ids);
    for m=1:length(strand_block_ids)
      len = blocks(strand_block_ids(m)).stop - blocks(strand_block_ids(m)).start +1;
      blocks(strand_block_ids(m)).content_pred = zeros(length(cont_pred_dirs),len);
    end
    for l = 1:length(cont_pred_dirs)
      SPF_contig_idx = strmatch(contigs{j}, {SPF_info{l}.contig_name}, 'exact') ;
      SPF_strand_idx = strmatch(strands{k}, {SPF_info{l}.strand}, 'exact') ;
      SPF_contig_strand_idx = intersect(SPF_contig_idx, SPF_strand_idx) ;
      assert(length(SPF_contig_strand_idx)==1) ;

      field_str = sprintf('SPF_%i', SPF_contig_strand_idx-1) ;
      LL=load(fnames{l}, field_str) ;
      LL=LL.(field_str) ;
      assert(isequal(LL.contig_name, contigs{j})) ;
      assert(isequal(LL.strand, strands{k})) ;
      %keyboard
      for m=1:length(strand_block_ids)
        start_pos = max(1, bsearch(LL.pos,  blocks(strand_block_ids(m)).start)-10) ;
        stop_pos = min(length(LL.pos), bsearch(LL.pos,  blocks(strand_block_ids(m)).stop)+10) ;
        idx11 = start_pos:stop_pos ;
        %[tmp,idx1,idx2] = intersect(LL.pos, blocks(strand_block_ids(m)).start:blocks(strand_block_ids(m)).stop) ;
        %idx11 = find(LL.pos>=blocks(strand_block_ids(m)).start & LL.pos<=blocks(strand_block_ids(m)).stop) ;
        [tmp, idx1, idx2] = intersect(LL.pos(idx11), blocks(strand_block_ids(m)).start:blocks(strand_block_ids(m)).stop) ;
        idx1=idx11(idx1) ;
        start_idx = idx11(find(LL.pos(idx11)==blocks(strand_block_ids(m)).start)) ;
        if isempty(start_idx),
          [tmp,start_idx] = min(abs(LL.pos(idx11)-blocks(strand_block_ids(m)).start)) ;
          start_idx = idx11(start_idx) ;
        end ;
        start_score = LL.score(start_idx) ;
        blocks(strand_block_ids(m)).content_pred(l,idx2) = LL.score(idx1) - start_score ;
      end 
    end
  end
end 
fprintf('Done.\n') ;

for i=1:length(blocks), 
  idx=find(isnan(blocks(i).content_pred)) ; 
  if ~isempty(idx), 
    fprintf('found %i nan/inf values in content predictions of block %i (%s, %c) => setting to 0\n', length(idx), i, blocks(i).chr, blocks(i).strand) ; 
    blocks(i).content_pred(idx) = 0 ; 
  end ;
end ;
