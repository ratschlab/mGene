function blocks = add_reads_from_bam(blocks, base_dir, which_data, tmp, maxval, filter)

if nargin<5
	maxval = 50;
end
if nargin<6
	%filter.intron = 100000; 
	filter.intron = 20000;
	filter.exon_len = 8;
	%filter.exon_len = 20;
	filter.mismatch = 1;
end

if isstruct(which_data)
	%% assume this is the genome info struct
	which_data = tmp;
end

types = separate(which_data, ',');
if isempty(types{1})
	fprintf('add_reads_from_bam: nothing to do\n');
	return
end

%reads from both strands are in one sam file
filenames = separate(base_dir, ',');

for b=1:length(blocks)
    clear introns

	if mod(b, 10)==0
		fprintf('\radd_exon_track_from_bam: %i(%i)', b, length(blocks));
	end
	block_len = blocks(b).stop-blocks(b).start+1;

	pair = 0;
	if any(strcmp(types, 'pair_coverage'))
		pair = 1;
	end

	%% get data from bam
	if any(strcmp(types, 'exon_track'))
		mapped = 1;
		spliced = 1;
		[introns, coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair);
	end
	if any(strcmp(types, 'mapped_exon_track'))
		mapped = 1;
		spliced = 0;
		[introns, mapped_coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair);
	end
	if any(strcmp(types, 'spliced_exon_track'))
		mapped = 0;
		spliced = 1;
		[introns, spliced_coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair);
	end	
	if any(strcmp(types, 'exon_intron_cov_opp'))
		% coverage from opposite strand
		mapped = 0;
		spliced = 1;
		bb = blocks(b);
		bb.strand = flip_strand(bb.strand);
		[introns_opp, spliced_coverage_opp, pair_cov] = get_all_data(bb, mapped, spliced, filenames, filter, pair);
		if blocks(b).strand=='+'
			introns_opp = introns_opp-blocks(b).start+1;
			introns_opp(:, 2) = introns_opp(:, 2)+1;
		elseif ~isempty(introns_opp)
			introns_opp = blocks(b).stop-introns_opp(:,2:-1:1)+1;
			introns_opp(:, 2) = introns_opp(:, 2)+1;
		end
	end	
	if ~exist('introns', 'var')
		% no exon coverage needed at all
		mapped = 0;
		spliced = 1;
		[introns, spliced_coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair);
	end

	%% process introns
	if blocks(b).strand=='+'
		introns = introns-blocks(b).start+1;
		introns(:, 2) = introns(:, 2)+1;
	elseif ~isempty(introns)
		introns = blocks(b).stop-introns(:,2:-1:1)+1;
		introns(:, 2) = introns(:, 2)+1;
	end
 
	% add requested data to block
	for j = 1:length(types)
		if strcmp(types{j}, 'pair_coverage')
			if ~isfield(blocks, 'tracks')
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end
			pair_cov(pair_cov>maxval) = maxval;
			if blocks(b).strand=='+' 
				blocks(b).tracks(end, :) = pair_cov;
			else
				blocks(b).tracks(end, :) = pair_cov(end:-1:1);
			end 	
		elseif strcmp(types{j}, 'exon_track')
			%% add exon track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			if ~isfield(blocks, 'tracks')
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end
			coverage(coverage>maxval) = maxval;
			if blocks(b).strand=='+' 
				blocks(b).tracks(end, :) = coverage;
			else
				blocks(b).tracks(end, :) = coverage(end:-1:1);
			end 
		elseif strcmp(types{j}, 'mapped_exon_track')
			%% add mapped exon track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			if ~isfield(blocks, 'tracks')
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end
			mapped_coverage(mapped_coverage>maxval) = maxval;
			if blocks(b).strand=='+' 
				blocks(b).tracks(end, :) = mapped_coverage;
			else
				blocks(b).tracks(end, :) = mapped_coverage(end:-1:1);
			end 
		elseif strcmp(types{j}, 'spliced_exon_track')
			%% add spliced exon track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			if ~isfield(blocks, 'tracks')
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end
			spliced_coverage(spliced_coverage>maxval) = maxval;
			if blocks(b).strand=='+' 
				blocks(b).tracks(end, :) = spliced_coverage;
			else
				blocks(b).tracks(end, :) = spliced_coverage(end:-1:1);
			end 
		elseif strcmp(types{j}, 'exon_intron_cov_opp')
			%% add spliced exon and intron track from opposite strand
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			if ~isfield(blocks, 'tracks')
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end
			intron_list = introns_opp;
			intron_coverage = zeros(1, block_len);
			if ~isempty(intron_list)
				for k = 1:size(intron_list, 1)
					from_pos = max(1, intron_list(k,1));
					to_pos = min(block_len, intron_list(k,2));
					intron_coverage(from_pos:to_pos) = intron_coverage(from_pos:to_pos)+1;
				end
			end
			if blocks(b).strand=='+'
				spliced_coverage_opp = spliced_coverage_opp+intron_coverage;
			else
				spliced_coverage_opp = spliced_coverage_opp(end:-1:1)+intron_coverage;
			end
			spliced_coverage_opp(spliced_coverage_opp>maxval) = maxval;
			blocks(b).tracks(end, :) = spliced_coverage_opp;
		elseif strcmp(types{j}, 'intron_track')
			%% add intron track to file
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			intron_list = introns;
			if ~isfield(blocks, 'tracks')
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end
		
			intron_coverage = zeros(1, block_len);
			if ~isempty(intron_list)
				for k = 1:size(intron_list, 1)
					from_pos = max(1, intron_list(k,1));
					to_pos = min(block_len, intron_list(k,2));
					intron_coverage(from_pos:to_pos) = intron_coverage(from_pos:to_pos)+1;
				end
			end
			intron_coverage(intron_coverage>maxval) = maxval;
			blocks(b).tracks(end, :) = intron_coverage;
		elseif strcmp(types{j}, 'intron_list')
			%% add intron list to file
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			intron_list = introns;
		
			if ~isempty(intron_list)
				rm_idx = find(intron_list(:,1)<1);
				rm_idx = [rm_idx; find(intron_list(:,2)>block_len)];
				intron_list(rm_idx, :) = [];
			end

			% calc number of occurences
			[tmp  fidx] = unique(intron_list, 'rows', 'first');
			[intron_list  lidx] = unique(intron_list, 'rows', 'last');
			intron_quality = lidx-fidx+1;
		
			if ~isempty(intron_list)
				[tmp sortidx] = sort(intron_list(:,2));
				intron_list = intron_list(sortidx, :);
				intron_quality = double(intron_quality(sortidx));
			end
            
            if isfield(filter, 'mincount') ;
                take_map=zeros(1,length(intron_list)) ;
                for kk=1:length(intron_quality)
                    take_map(kk)=intron_quality(kk)>=filter.mincount ;
                end ;
                %fprintf('dropped %i introns, keeping %i introns\n', sum(take_map==0), sum(take_map~=0)) ;
                intron_list=intron_list(find(take_map~=0),:) ;
                intron_quality=intron_quality(find(take_map~=0),:) ;
            end ;
		
			if ~isfield(blocks, 'segment_lists')
				blocks(b).segment_lists = {intron_list};
				blocks(b).segment_scores = {intron_quality};
			else
				blocks(b).segment_lists{end+1} = intron_list;
				blocks(b).segment_scores{end+1} = intron_quality;
			end
		else 
			fprintf('unknown type of data: %s', types{j})
		end
	end
end

return

function [introns, coverage, pair_cov] = get_all_data(block, mapped, spliced, filenames, filter, pair) 

	block_len = block.stop-block.start+1;
	% get all data from bam file
	coverage = zeros(1, block_len);
	pair_cov = zeros(1, block_len);
	introns = [];
	for j = 1:length(filenames)
		fname = filenames{j};
		if ~(exist(fname)==2)
			fprintf('add_exon_track_from_bam: did not find file %s\n', fname);
			continue;
		end

		subsample = 1000;%% no subsampling
		if isfield(block, 'subsample_reads')
			subsample = floor(block.subsample_reads*1000);
		end

		contig_name = block.chr;
		strand = block.strand;
		VERBOSE = 0;
		collapse = 1;
		maxminlen = 0;
		%warning('constants changed')
		if pair
			[coverage_tmp, introns_cell, pair_cov_tmp] = get_reads(fname, contig_name, block.start, block.stop, strand, collapse, subsample, filter.intron, filter.exon_len, filter.mismatch, mapped, spliced, maxminlen, pair);
			pair_cov = pair_cov+double(pair_cov_tmp);
		else
			[coverage_tmp, introns_cell] = get_reads(fname, contig_name, block.start, block.stop, strand, collapse, subsample, filter.intron, filter.exon_len, filter.mismatch, mapped, spliced);
		end
		coverage = coverage+double(coverage_tmp);
		if iscell(introns_cell)
			introns = [introns; double([introns_cell{:}]')];
		else
			introns = [introns; double(introns_cell)'];
		end
	end
return

