function [SEG]=move_cleave_pos(SEG, block, model)
% 
% if there is no feature position between last segment end 
% and cleavage site the the cleavage site is moved 

if isempty(SEG)
	return
end

col4 = SEG.segments(:,4);
num_genes = length(unique(col4))-1;
assert(all(unique(col4)'==0:num_genes))


for i=1:num_genes
	gene_lines = find(col4==i);
	segments = SEG.segments(gene_lines,1:3);


	% find cleave position
	%------------------------------------------------------------
	cleave_idx = find(segments(:,3)== model.segments.polya_tail,1,'last');
	if isempty(cleave_idx)
		cleave_idx = find(segments(:,3)== model.segments.utr3exon,1,'last');
	end
	cleave_pos = segments(cleave_idx,2);
	last_pos = segments(cleave_idx,1);
	
	pp = find_nearest_pos(cleave_pos-1, block.all_pos, last_pos+1, cleave_pos-1);

	if isempty(pp)
		%% next segment end
		ub = SEG.segments(min(gene_lines)-1+cleave_idx+1,2);

		%% move cleave to the next cleavage prediciton
		lb = cleave_pos+1; 	
		canditates = block.Signals.cleave.pos;
		new_pos = find_nearest_pos(cleave_pos+2, canditates, lb, ub);

		if isempty(new_pos)||new_pos-cleave_pos>50%too far away
			SEG = []; 
			return, 
		end
		SEG.segments(min(gene_lines)-1+cleave_idx,2) = new_pos;
		SEG.segments(min(gene_lines)-1+cleave_idx+1,1) = new_pos;

		upos = unique(SEG.segments(gene_lines,1:2));
		assert(length(intersect(block.all_pos,upos))==length(upos));
	end
end

return


