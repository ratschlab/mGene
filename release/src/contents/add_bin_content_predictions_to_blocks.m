function blocks = add_bin_content_predictions_to_blocks(blocks, cont_pred_dirs)
% blocks = add_content_predictions_to_blocks(blocks, cont_pred_dirs)

fields = {'output'};
contigs = unique({blocks.chr}) ;

%for i=1:length(cont_pred_dirs),
%  fnames{i}=  sprintf('%s/output_spf.mat', cont_pred_dirs{i}) ;
%  L = load(fnames{i}, 'SPF_info') ;
%  SPF_info{i} = L.SPF_info ;
%end ;

fprintf('Processing %i contigs: \n', length(contigs)) ;
for j=1:length(contigs)
  if mod(j, 5)==0 || length(contigs)<60, 
    fprintf('%i ... ', j) ; 
  end ;

  contig_block_ids = find(ismember({blocks.chr}, contigs{j}));
  strands = unique({blocks.strand});
  for k=1:length(strands)
    strand_block_ids = intersect(strmatch(strands(k), {blocks.strand}, 'exact'), contig_block_ids);
    for m=strand_block_ids
      len = blocks(m).stop - blocks(m).start + 1;
      blocks(m).content_pred = zeros(length(cont_pred_dirs),len);
    end
    for l = 1:length(cont_pred_dirs)
      %SPF_contig_idx = strmatch(contigs{j}, {SPF_info{l}.contig_name}, 'exact') ;
      %SPF_strand_idx = strmatch(strands{k}, {SPF_info{l}.strand}, 'exact') ;
      %SPF_contig_strand_idx = intersect(SPF_contig_idx, SPF_strand_idx) ;
      %assert(length(SPF_contig_strand_idx)==1) ;

      %field_str = sprintf('SPF_%i', SPF_contig_strand_idx-1) ;
      %LL=load(fnames{l}, field_str) ;
      %LL=LL.(field_str) ;
      %assert(isequal(LL.contig_name, contigs{j})) ;
      %assert(isequal(LL.strand, strands{k})) ;
      %keyboard
      for m = strand_block_ids
        %start_pos = max(1, bsearch(LL.pos,  blocks(m).start)-10) ;
        %stop_pos = min(length(LL.pos), bsearch(LL.pos,  blocks(m).stop)+10) ;
        %idx11 = start_pos:stop_pos ;
        %[tmp,idx1,idx2] = intersect(LL.pos, blocks(m).start:blocks(strand_block_ids(m)).stop) ;
        %idx11 = find(LL.pos>=blocks(m).start & LL.pos<=blocks(strand_block_ids(m)).stop) ;
        %[tmp, idx1, idx2] = intersect(LL.pos(idx11), blocks(m).start:blocks(strand_block_ids(m)).stop) ;
        %idx1=idx11(idx1) ;
        %start_idx = idx11(find(LL.pos(idx11)==blocks(m).start)) ;
        %if isempty(start_idx),
        %  [tmp,start_idx] = min(abs(LL.pos(idx11)-blocks(m).start)) ;
        %  start_idx = idx11(start_idx) ;
        %end ;
        %start_score = LL.score(start_idx) ;
        %blocks(m).content_pred(l,idx2) = LL.score(idx1) - start_score;
        f_name=sprintf('%s/pred/contig_%i%s', cont_pred_dirs{l}, blocks(m).chr_num, strands{k});
        if ~fexist(sprintf('%s.pos', f_name))
          error('file %s.pos not found.', f_name);
        end
        [contig_pos,score] = interval_query(f_name,fields,[blocks(m).start;blocks(m).stop]) ;
        [tmp,idx1,idx2] = intersect(contig_pos, blocks(m).start:blocks(m).stop) ;
		score = score(idx1) ;
        blocks(m).content_pred(l,idx2) = score(:)  - score(1) ;
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
