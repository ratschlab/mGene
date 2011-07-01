function [path,pos_idx,pos]=segmentation2path(SEG, block, model)
% [path,pos_idx,pos]=segmentation2path(SEG, block, model)
	

%pos = [0		 length(block.all_pos)-1] ;
%path = [model.state_ids.seq_start		model.state_ids.seq_end]-1 ;
	

genestr=lower(block.seq) ;

col4 = SEG.segments(:,4);
num_genes = length(unique(col4))-1;
assert(all(unique(col4)'==0:num_genes))

path = [] ;
pos = [] ;
path(end+1) = model.state_ids.seq_start ;
pos(end+1) = block.all_pos(1) ;

any_invalid = 0 ;
for i=1:num_genes
	segments = SEG.segments(col4==i,1:3);

	% path generation for non coding genes
	idx = find(segments(:,3)==model.segments.cds_exon) ;
	if isempty(idx)
		[path pos] = define_nc_path(path, pos, segments, model);
		continue;
	end

	% dealing with invalid genes
	if (size(segments,1)==1 & isnan(segments(1,3)))
		path(end+1:end+2) = [model.state_ids.tss(1) model.state_ids.cdsStop(1)] ;
		pos(end+1:end+2) = segments(1:2)';
		warning('gene is invalid') ;
		any_invalid = 1 ;
		continue	 
	end

	% consider 5' utr
	idx = find(segments(:,3)==model.segments.utr5exon) ;
	if ~isempty(idx),
		utr5_segments = segments(idx,:) ;
		path(end+1) = model.state_ids.tss(1) ;
		pos(end+1) = utr5_segments(1,1) ;
		
		for i=1:size(utr5_segments,1)-1
			path(end+1)= model.state_ids.don(1) ;
			pos(end+1) = utr5_segments(i, 2) ;
			path(end+1)= model.state_ids.acc(1) ;
			pos(end+1) = utr5_segments(i+1, 1) ;
		end ;
	end ;

	if isfield(model.segments,'transexon')
		idx = find(segments(:,3)==model.segments.transexon) ;
	else 
		idx = [];
	end
	if ~isempty(idx),
		if length(idx)==1 && segments(idx(1),1)==segments(idx(1),2),
			% this is the trans-splice agatg case
			path(end+1) = model.state_ids.trans(2) ;
			pos(end+1) = segments(idx,1) ;
		else
			utr5_segments = segments(idx,:) ;
			path(end+1) = model.state_ids.trans(1) ;
			pos(end+1) = utr5_segments(1,1) ;
			
			for i=1:size(utr5_segments,1)-1
				path(end+1)= model.state_ids.don(1) ;
				pos(end+1) = utr5_segments(i, 2) ;
				path(end+1)= model.state_ids.acc(1) ;
				pos(end+1) = utr5_segments(i+1, 1) ;
			end ;
		end ;
	end ;

	% consider coding region
	idx = find(segments(:,3)==model.segments.cds_exon) ;
	coding_segments = segments(idx,:) ;
	if idx(1)>1 && segments(idx(1)-1,3)==model.segments.intron,
		% special case agatg
		path(end+1) = model.state_ids.don(1) ;
		pos(end+1) = segments(idx(1)-1,1) ;
		path(end+1) = model.state_ids.tis(2) ;
		pos(end+1) = segments(idx(1)-1,2) ;
	else
		if ~isfield(model.state_ids,'trans')||path(end)~=model.state_ids.trans(2),
			path(end+1) = model.state_ids.tis(1) ;
			pos(end+1) = coding_segments(1,1) ;
		end ;
	end ;
	%spliced = '' ;
	%for i=1:size(coding_segments,1)
	%	spliced = [spliced genestr(coding_segments(i,1):coding_segments(i,2)-1)] ;
	%end ;
	%assert(mod(length(spliced),3)==0) ;
	%f=find_open_frames(spliced,1)

	frame = mod(coding_segments(1,2)-coding_segments(1,1),3) ;
	for i=1:size(coding_segments,1)-1
		pos(end+1) = coding_segments(i,2) ;
		pos(end+1) = coding_segments(i+1,1) ;
		if frame==0,
			path(end+1)= model.state_ids.don(1+1) ;
			path(end+1)= model.state_ids.acc(1+1) ;
		elseif frame==1,
			if genestr(coding_segments(i,2)-1)~='t'
				path(end+1)=model.state_ids.don(1+2) ;
				path(end+1)=model.state_ids.acc(1+2) ;
			else
				path(end+1)=model.state_ids.don(1+3) ;
				path(end+1)=model.state_ids.acc(1+3) ;
			end ;
		elseif frame==2,
			if genestr(coding_segments(i,2)-2)=='t' & genestr(coding_segments(i,2)-1)=='a'
				path(end+1)=model.state_ids.don(1+5) ;
				path(end+1)=model.state_ids.acc(1+5) ;
			elseif genestr(coding_segments(i,2)-2)=='t' & genestr(coding_segments(i,2)-1)=='g'
				path(end+1)=model.state_ids.don(1+6) ;
				path(end+1)=model.state_ids.acc(1+6) ;
			else
				path(end+1)=model.state_ids.don(1+4) ;
				path(end+1)=model.state_ids.acc(1+4) ;
			end ;
		end ;
		frame = mod(frame+coding_segments(i+1,2)-coding_segments(i+1,1),3) ;
	end ;
	
	% consider 3' utr
	utr3_segments = [] ;

	idx = find(segments(:,3)==model.segments.cds_exon) ;
	if segments(idx(end)+1,3)==model.segments.intron 
		% special case: agstop
		path(end+1)= model.state_ids.don(2) ;
		pos(end+1) = segments(idx(end)+1,1) ;
		path(end+1)= model.state_ids.cdsStop(2) ;
		pos(end+1) = segments(idx(end)+1,2) ;
	else
		% normal stop
		pos(end+1)=coding_segments(end,2) ;
		path(end+1)=model.state_ids.cdsStop(1) ;
	end ;

	idx = find(segments(:,3)==model.segments.utr3exon) ;
	if ~isempty(idx),
		utr3_segments = segments(idx,:) ;
		for i = 1:size(utr3_segments,1)-1
			path(end+1)= model.state_ids.don(end) ;
			pos(end+1) = utr3_segments(i, 2) ;
			path(end+1)= model.state_ids.acc(end) ;
			pos(end+1) = utr3_segments(i+1, 1) ;
		end ;
	end ;

	% consider 3' utr
	idx = find(segments(:,3)==model.segments.polya_tail) ;
	if ~isempty(idx),
		assert(length(idx)==1) ;
		utr3_segments = segments(idx,:) ;
		path(end+1) = model.state_ids.polya(1) ;
		pos(end+1) = utr3_segments(1,1) ;
	end ;

	if isfield(model.state_ids, 'rna_seq_polya')
		if ~isempty(utr3_segments)
			len = utr3_segments(end,2)-pos(end);
	
			if len>=10
				lb = pos(end)+1;
				ub = utr3_segments(end,2)-1;
				pp = find_nearest_pos(utr3_segments(end,2)-len+1, block.all_pos, lb, ub);
			elseif len>=2
				lb = pos(end)+1;
				ub = utr3_segments(end,2)-1;
				pp = find_nearest_pos(pos(end)+1, block.all_pos, lb, ub);
			else
				pp = [];
			end
			if isempty(pp)
				error('segmentation2path: could not assign rna_seq_polya state\n');
			end
			path(end+1) = model.state_ids.rna_seq_polya;
			pos(end+1) = pp;
		else
			error('segmentation2path: utr3_segments empty');
		end
	end

	if ~isempty(utr3_segments),
		path(end+1) = model.state_ids.cleave ;
		pos(end+1) = utr3_segments(end,2) ;
	end ;
end ;

path(end+1) = model.state_ids.seq_end ;
pos(end+1)	= block.all_pos(end) ;


if ~any_invalid,
	for i=1:length(pos)-1
		assert(model.A(path(i),path(i+1))>=-1000) ;
	end ;
end ;


[tmp, pos_idx, idx] = intersect(block.all_pos, pos) ;

assert(isempty(setdiff(pos, block.all_pos(pos_idx))))
%if ~isempty(setdiff(pos, block.all_pos(pos_idx)))
%	[tmp,idx1]=setdiff(pos, block.all_pos(pos_idx)) ;
%	warning('%i (%i) position(s) not found', length(idx1),length(pos)) ;
%	keyboard
%end ;
path = path(idx) ;
pos = block.all_pos(pos_idx) ;

%%%%

return 


function [path pos] = define_nc_path(path, pos, segments, model)

path(end+1) = model.state_ids.tss_nc;
pos(end+1) = segments(1, 1);

idx = find(segments(:, 3)==model.segments.intron); 

for j = idx' 
	path(end+1) = model.state_ids.don_nc;
	pos(end+1) = segments(j, 1);
	path(end+1) = model.state_ids.acc_nc;
	pos(end+1) = segments(j, 2);
end	
path(end+1) = model.state_ids.cleave_nc;
pos(end+1) = segments(end, 2);

return




