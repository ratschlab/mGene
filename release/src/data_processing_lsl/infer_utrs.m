function genes = infer_utrs(blocks, genes, PAR, track_function, track_file)

genome_info = init_genome(PAR.FN.genome.fn_genome_config);

if nargin==5
	blocks = feval(track_function, blocks, track_file, 'exon_track,intron_track');
end
if size(blocks(1).tracks, 1)~=1
	for j = 1:length(blocks)
		blocks(j).tracks=sum(blocks(j).tracks);
	end
end
for j=1:length(blocks)
	gene_id = blocks(j).gene_ids;
	assert(length(gene_id)==1)
	assert(genes(gene_id).id==gene_id)
	base_dir = PAR.FN.output_sig.tss.fn_pred;
	f_name_tss = sprintf('%scontig_%i%s', base_dir, genes(gene_id).chr_num, genes(gene_id).strand );
	base_dir = PAR.FN.output_sig.cleave.fn_pred;
	f_name_cleave = sprintf('%scontig_%i%s', base_dir, genes(gene_id).chr_num, genes(gene_id).strand );
	genes(gene_id) = find_tss_and_cleave(blocks(j), genes(gene_id), f_name_tss, f_name_cleave);
end

which_checks.exons_sorted = 1;
which_checks.intron_length = 1;
which_checks.exon_length = 1;
which_checks.splicesites = 1;
which_checks.orf = 1;
which_checks.gene_length = 0;
which_checks.graph = 0;
which_checks.transacc = 0;
which_checks.complete = 1;
gg = genes([blocks.gene_ids]);
num_complete = sum([gg.is_complete]);
gg = check_genes(gg, genome_info, which_checks); 
fprintf('now %i (%i) genes are complete (%i before)\n', sum([gg.is_complete]), length(gg), num_complete)
genes([blocks.gene_ids]) = gg;

return


function gene = find_tss_and_cleave(block, gene, f_name_tss, f_name_cleave)

	if ~fexist([f_name_tss, '.pos'])||~fexist([f_name_cleave, '.pos'])
		error('file not found');
	end
	for k = 1:length(gene.transcripts)
		if length(gene.tis)<k||isempty(gene.tis{k})
			continue;
		end
		if length(gene.tss)>=k&&~isempty(gene.tss{k})
			continue
		end
		%% find tss position
		if gene.strand=='+'
			if isfield(block, 'tracks')
				[lower_bound, upper_bound] = find_upstream_boundaries(block, gene.tis{k}-1, block.tracks);
			else
				lower_bound = max(block.start+1, gene.tis{k}-500);
				upper_bound = gene.tis{k}-1;
			end
			[all_tss_pos, all_tss_score] = interval_query(f_name_tss, {'Conf_cum'} ,[lower_bound;upper_bound]);
			[max_score, max_idx] = max(all_tss_score);
			gene.tss{k} = all_tss_pos(max_idx);
			gene.utr5_exons{k} = [gene.tss{k} gene.tis{k}];
			gene.start = min(gene.start, gene.tss{k});
		else
			if isfield(block, 'tracks')
				% revert track because all other coordinates are on watson strand
				[lower_bound, upper_bound] = find_downstream_boundaries(block, gene.tis{k}+1, block.tracks(end:-1:1));
			else
				upper_bound = min(block.stop-1, gene.tis{k}+500);
				lower_bound = gene.tis{k}+1;
			end
			[all_tss_pos, all_tss_score] = interval_query(f_name_tss, {'Conf_cum'} ,[lower_bound;upper_bound]);	
			[max_score, max_idx] = max(all_tss_score);
			gene.tss{k} = all_tss_pos(max_idx);
			gene.utr5_exons{k} = [gene.tis{k} gene.tss{k}];
			gene.stop = max(gene.stop, gene.tss{k});
		end
	end

	for k = 1:length(gene.transcripts)
		if length(gene.cdsStop)<k||isempty(gene.cdsStop{k})
			continue;
		end
		if length(gene.cleave)>=k&&~isempty(gene.cleave{k})
			continue
		end		
		%% find cleave position
		if gene.strand=='+'
			if isfield(block, 'tracks')
				[lower_bound, upper_bound] = find_downstream_boundaries(block, gene.cds_exons{k}(end, 2)+1, block.tracks);
			else
				upper_bound = min(block.stop-1, gene.cds_exons{k}(end, 2)+500);
				lower_bound = gene.cds_exons{k}(end, 2)+1;
			end
			[all_cleave_pos, all_cleave_score] = interval_query(f_name_cleave, {'Conf_cum'} ,[lower_bound;upper_bound]);
			[max_score, max_idx] = max(all_cleave_score);
			gene.cleave{k} = all_cleave_pos(max_idx);
			gene.utr3_exons{k} = [gene.cds_exons{k}(end, 2) gene.cleave{k}];
			gene.stop = max(gene.stop, gene.cleave{k});
		else
			if isfield(block, 'tracks')
				% revert track because all other coordinates are on watson strand
				[lower_bound, upper_bound] = find_upstream_boundaries(block, gene.cds_exons{k}(1,1)-1, block.tracks(end:-1:1));
			else
				upper_bound = gene.cds_exons{k}(1,1)-1; 
				lower_bound = max(block.start+1, gene.cds_exons{k}(1,1)-500);
			end
			[all_cleave_pos, all_cleave_score] = interval_query(f_name_cleave, {'Conf_cum'} ,[lower_bound;upper_bound]);	
			[max_score, max_idx] = max(all_cleave_score);
			gene.cleave{k} = all_cleave_pos(max_idx);
			gene.utr3_exons{k} = [gene.cleave{k} gene.cds_exons{k}(1,1)];
			gene.start = min(gene.start, gene.cleave{k});
		end

	end 
return

function[lower_bound, upper_bound] = find_upstream_boundaries(block, gene_start, track);

	assert(size(track, 1)==1)
	threshold = 3;
	win = 10;
	for j = gene_start:-1:block.start+win+1
		track_pos = (j-win:j)-block.start+1;% local coordinates
		win_cov = mean(track(track_pos));
		if j == gene_start
			start_cov = win_cov;
		end
		if win_cov<threshold
			break
		end
		if win_cov<start_cov/2 %drop in coverage
			break
		end
	end
	upper_bound = j;
	lower_bound = max(block.start, upper_bound-200);
	if upper_bound-lower_bound<50%some freedom to choose good signal
		upper_bound = min(gene_start, lower_bound+50);
	end
return
function[lower_bound, upper_bound] = find_downstream_boundaries(block, gene_end, track);

	assert(size(track, 1)==1)
	threshold = 3;
	win = 10;
	for j = gene_end:block.stop-win-1
		track_pos = (j:j+win)-block.start+1;% local coordinates
		win_cov = mean(track(track_pos));
		if j == gene_end
			start_cov = win_cov;
		end
		if win_cov<threshold
			break
		end
		if win_cov<start_cov/2 %drop in coverage
			break
		end
	end
	lower_bound = j;
	upper_bound = min(block.stop, lower_bound+200);
	if upper_bound-lower_bound<50%some freedom to choose good signal
		lower_bound = max(gene_end, upper_bound-50);
	end
return



