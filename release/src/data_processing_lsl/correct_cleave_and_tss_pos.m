function [SEG]=correct_cleave_and_tss_pos(SEG, block, model)
% 
% if there are no predictions for the annotated cleave and tss positions then 
% set cleave and tss to the nearest predicted  position 

col4 = SEG.segments(:,4);
num_genes = length(unique(col4))-1;
assert(all(unique(col4)'==0:num_genes))


for i=1:num_genes
  gene_lines = find(col4==i);
  segments = SEG.segments(gene_lines,1:3);

  upos = unique(segments(:,1:2));
  compos = intersect(block.all_pos,upos);
  if 1
    diffpos = setdiff(upos,compos);

    % find cleave position
    %------------------------------------------------------------
    cleave_idx = find(segments(:,3)== model.segments.polya_tail,1,'last');
    if isempty(cleave_idx)
      cleave_idx = find(segments(:,3)== model.segments.utr3exon,1,'last');
    end
	if isfield(model.use, 'non_coding')&&isempty(cleave_idx)
		cleave_idx = find(segments(:,3)== model.segments.nc_exon,1,'last');
	end
    cleave_pos = segments(cleave_idx,2);

    lower_bound = segments(cleave_idx,1)+1;
    upper_bound = SEG.segments(min(gene_lines)-1+cleave_idx+1,2)-1;

    new_pos = find_nearest_pos(cleave_pos, block.Signals.cleave.pos, ...
                               lower_bound, upper_bound);
    if isempty(new_pos), SEG = []; return, end
    SEG.segments(min(gene_lines)-1+cleave_idx,2) = new_pos;
    SEG.segments(min(gene_lines)-1+cleave_idx+1,1) = new_pos;

    % find tss position
    %------------------------------------------------------------
    tss_idx = find(segments(:,3)== model.segments.utr5exon,1,'first');
	if isfield(model.use, 'non_coding')&&isempty(tss_idx)
		tss_idx = find(segments(:,3)== model.segments.nc_exon,1,'first');
	end
    tss_pos = segments(tss_idx,1);

    lower_bound = SEG.segments(min(gene_lines)-1+tss_idx-1,1)+1;
    upper_bound = segments(tss_idx,2)-1;

    new_pos = find_nearest_pos(tss_pos, block.Signals.tss.pos, lower_bound, upper_bound);
    if isempty(new_pos), SEG = []; return, end
    SEG.segments(min(gene_lines)-1+tss_idx,1) = new_pos;
    SEG.segments(min(gene_lines)-1+tss_idx-1,2) = new_pos;

    if ~(all(ismember(diffpos,cleave_pos)|ismember(diffpos, tss_pos)))
      SEG = []; 
      return
    end
    upos = unique(SEG.segments(gene_lines,1:2));
    assert(length(intersect(block.all_pos,upos))==length(upos));
  end
end

return


