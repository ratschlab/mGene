function [segments,genes,ok] = path2segmentation(my_path,my_pos_idx,block, model)
% [segmentation,ok] = path2segmentation(path,pos_idx,block, model)
%pos = [0     length(block.all_pos)-1] ;
%path = [model.state_ids.seq_start    model.state_ids.seq_end]-1 ;
  
segment_names = fieldnames(model.segments) ;
segments = [] ;
ok = 1 ;
genes = {} ;

if length(my_pos_idx)>length(unique(my_pos_idx))
  ok=0;
  return
end
if my_path(1)~=model.state_ids.seq_start
  ok=0;
  return
end

%my_pos_idx(my_path==getfield(model.state_ids,'intergenic_long'))=[];
%my_path(my_path==getfield(model.state_ids,'intergenic_long'))=[];

%if isfield(model.state_ids, 'rna_seq_polya')
%	idxp = find(my_path==getfield(model.state_ids,'rna_seq_polya'));
%	my_pos_idx(idxp)=[];
%	my_path(idxp)=[];
%end


% same for intron long
if isfield(model.state_ids,'intron_long')
  for i=getfield(model.state_ids,'intron_long')
    del_idx=my_path==i;
    my_pos_idx(del_idx)=[];
    my_path(del_idx)=[];
  end
end
for i=1:length(my_path)-1
  lengths_idx = model.transition_pointers(my_path(i+1), my_path(i)) ;
  if my_path(i+1) == model.state_ids.tis(2),
    if ~(my_path(i) == model.state_ids.don(1)) ;
      error('invalid path') ;
      %warning('invalid path') ;
      %ok=0 ;
      %return ;
    end ;
    % handle the agatg-case properly 
    segment_type = model.seg_links(find(model.seg_links(:,1)==lengths_idx),2) ;
    segments = [segments; [block.all_pos(my_pos_idx(i)) block.all_pos(my_pos_idx(i+1)), model.segments.intron ]] ;
    continue ;
  end ;
  
  if lengths_idx~=0
    segment_type = model.seg_links(find(model.seg_links(:,1)==lengths_idx),2) ;
    % hack
    if isfield(model.segments,'intercistronic')&&segment_type==model.segments.intercistronic & ...
          block.all_pos(my_pos_idx(i+1)) - block.all_pos(my_pos_idx(i)) > model.lengths_range.intercistronic(2),
      segment_type = model.segments.intergenic ;
	elseif isfield(model.segments,'rna_seq_polya')&&segment_type==model.segments.rna_seq_polya
		assert(segments(end, 3)==model.segments.utr3exon)
		segments(end, 2) = block.all_pos(my_pos_idx(i+1));
		continue;
    end 
  else
    segment_type = 0 ;
    warning('unable to decode path') ;
    ok = 0 ;
  end ;
  segments = [segments; [block.all_pos(my_pos_idx(i)) block.all_pos(my_pos_idx(i+1)),segment_type ]] ;
end

segments(1,1) = 1;
segments(end,2) = length(block.seq);
idx =  find(segments(:,3) ==getfield(model.segments,'intergenic'));
if isfield(model.segments,'intercistronic')
  idx = [idx find(segments(:,3) ==getfield(model.segments,'intercistronic'))];
end
n_genes = length(idx)-1 ;
genes{1} = [];
for i=1:n_genes
  genes{i} = segments(idx(i)+1:idx(i+1)-1,:);
end
