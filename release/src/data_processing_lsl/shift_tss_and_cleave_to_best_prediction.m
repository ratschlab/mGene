function [SEG]=shift_tss_to_best_prediction(SEG, block, model, win)

if nargin<4
	win = 200;
end
col4 = SEG.segments(:,4);
num_genes = length(unique(col4))-1;
assert(all(unique(col4)'==0:num_genes))


for i=1:num_genes
  gene_lines = find(col4==i);
  segments = SEG.segments(gene_lines,1:3);

  if 1
    % find tss position
    %------------------------------------------------------------
    tss_idx = find(segments(:,3) == model.segments.utr5exon);
	if isfield(model.use, 'non_coding')&&isempty(tss_idx)
		tss_idx = find(segments(:,3)== model.segments.nc_exon,1,'first');
	end
    tss_pos = segments(tss_idx(1),1);

	% find tss predictions in upstream intergenic segment
    %------------------------------------------------------------
    lower_bound = SEG.segments(min(gene_lines)-1+tss_idx(1)-1,1)+1; % end of last gene
    upper_bound = segments(tss_idx(1),2)-1;% end of 5'UTR
    lower_bound = max(lower_bound+50, tss_pos-win);

	if lower_bound>=upper_bound % intergenic region to short
		continue
	end

    all_tss_idx = find(block.Signals.tss.pos>lower_bound & block.Signals.tss.pos<upper_bound);

	% favour tss predictions that are close to the upper bound	
	all_tss_pos = block.Signals.tss.pos(all_tss_idx);
	linear_penalty= ((upper_bound-all_tss_pos)*0.8+all_tss_pos-lower_bound)/(upper_bound-lower_bound);

	[max_score, max_idx] = max(block.Signals.tss.Conf_cum(all_tss_idx).*linear_penalty);
    if isempty(all_tss_idx)||max_score<0.6
		continue
	end

	% shift tss position
    %------------------------------------------------------------
	new_pos = block.Signals.tss.pos(all_tss_idx(max_idx));
    SEG.segments(min(gene_lines)-1+tss_idx(1),1) = new_pos;
    SEG.segments(min(gene_lines)-1+tss_idx(1)-1,2) = new_pos;

	%fprintf('moved positon of tss from %i to %i\n',tss_pos, new_pos)
  end
    % find cleave position
    %------------------------------------------------------------
	cleave_idx = [];
	if isfield(model.segments, 'polya_tail')
    	cleave_idx = find(segments(:,3)== model.segments.polya_tail,1,'last');
	end
    if isempty(cleave_idx)
      cleave_idx = find(segments(:,3)== model.segments.utr3exon,1,'last');
    end
	if isfield(model.use, 'non_coding')&&isempty(cleave_idx)
		cleave_idx = find(segments(:,3)== model.segments.nc_exon,1,'last');
	end

    cleave_pos = segments(cleave_idx,2);

    lower_bound = segments(cleave_idx,1)+1;% start of the last 3'UTR exon
    upper_bound = SEG.segments(min(gene_lines)-1+cleave_idx+1,2)-1; % start of next gene
	upper_bound = min(cleave_pos+win, upper_bound-50);

	all_cleave_idx = find(block.Signals.cleave.pos>lower_bound & block.Signals.cleave.pos<upper_bound);
	% favour cleave predictions that are close to the upper bound	
	all_cleave_pos = block.Signals.cleave.pos(all_cleave_idx);
	linear_penalty= ((upper_bound-all_cleave_pos)*0.8+all_cleave_pos-lower_bound)/(upper_bound-lower_bound);

	[max_score, max_idx] = max(block.Signals.cleave.Conf_cum(all_cleave_idx).*linear_penalty);
    if isempty(all_cleave_idx)||max_score<0.6
		continue
	end

	% shift cleave position
	%------------------------------------------------------------
    new_pos = block.Signals.cleave.pos(all_cleave_idx(max_idx)); 
    SEG.segments(min(gene_lines)-1+cleave_idx,2) = new_pos;
    SEG.segments(min(gene_lines)-1+cleave_idx+1,1) = new_pos;

end

return


