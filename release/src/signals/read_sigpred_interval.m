function [chr, strand, pos, score]=read_sigpred_interval(sigpred_fname, genome_info)

chr=[] ; 
strand = '' ;
pos = [] ;
score = [] ;

for j = 1:length(genome_info.contig_names)
	len = get_chr_len(genome_info, j);
	for s = '+-'
		f_name = sprintf('%s/contig_%i%s', sigpred_fname, j, s);

		if fexist([f_name '.pos'])&&fexist([f_name '.label'])
			[contig_pos, contig_score] = interval_query(f_name, {'label'},[1;len]);
		else
			error('file not found: %s ', f_name);
		end
		chr = [chr j*ones(1, length(contig_pos))];
		pos = [pos contig_pos'];
		score = [score contig_score'];
		strand = [strand repmat(s, 1, length(contig_pos))];
	end	
end
