function compare_aa_fasta(fasta1, fasta2, fn_out)

if nargin==3
	[tmp names1]=read_fasta_file(fasta1, 1);
	range = 1:1000:length(names1)
	if length(range)==1
		range(end+1) = inf;
	end
	cnt = 0;
	jobinfo = rproc_empty();
	for j = 1:length(range)-1
		PAR.fasta1 = fasta1;
		PAR.fasta2 = fasta2; 
		PAR.fn_out = fn_out; 
		PAR.range = [range(j), range(j+1)-1];
		if j == length(range)-1
			PAR.range(2) = inf;
		end
		if 0%fexist(sprintf('%s%i.mat', PAR.fn_out, PAR.range(1)))
			continue
		end
		opts.addpaths = {fileparts(which('compare_aa_fasta')), fileparts(which('get_strain_names'))};
		opts.waitonfull = 1;
		opts.priority = 57;
		opts.maxjobs = 2500;
		cnt = cnt+1;
		jobinfo(cnt) = rproc('compare_aa_fasta', PAR, 5000, opts, 100000);
		%compare_aa_fasta(PAR);
	end
	jobinfo = rproc_wait(jobinfo, 20, 1, -1);

	% collect results
	load(sprintf('%s%i.mat', fn_out, range(1)), 'mindist', 'd01_cnt', 'match_names', 'd01_match_names', 'diff_match_names', 'mindist_diff')
	for j = 2:length(range)-1
		try
			l = load(sprintf('%s%i.mat', fn_out, range(j)), 'mindist', 'd01_cnt', 'match_names', 'd01_match_names', 'diff_match_names', 'mindist_diff');
		catch
			sprintf('%s%i.mat', fn_out, range(j))
			warning('file not found')
		end
		mindist(range(j):range(j+1)-1) = l.mindist(range(j):range(j+1)-1);	
		mindist_diff(range(j):range(j+1)-1) = l.mindist_diff(range(j):range(j+1)-1);	
		d01_cnt(range(j):range(j+1)-1) = l.d01_cnt(range(j):range(j+1)-1);
		match_names(range(j):range(j+1)-1) = l.match_names(range(j):range(j+1)-1);
		d01_match_names(range(j):range(j+1)-1) = l.d01_match_names(range(j):range(j+1)-1);
		diff_match_names(range(j):range(j+1)-1) = l.diff_match_names(range(j):range(j+1)-1);
	end
	save(fn_out, 'mindist', 'd01_cnt', 'match_names', 'd01_match_names', 'diff_match_names', 'mindist_diff')
	for j = 2:length(range)-1
		unix(sprintf('rm -f %s%i.mat', fn_out, range(j)));
	end
	return
elseif nargin==1
	fasta2 = fasta1.fasta2;
	range = fasta1.range;
	fn_out = fasta1.fn_out;
	fasta1 = fasta1.fasta1;
	paths
else
	error('wrong number of args')
end

addpath /fml/ag-raetsch/home/raetsch/seqan/demos/

[seqs1,names1]=read_fasta_file(fasta1);
[seqs2,names2]=read_fasta_file(fasta2);

mindist = 2*ones(1, length(seqs1));
mindist_diff = 2*ones(1, length(seqs1));
match_names = cell(1, length(seqs1));
d01_match_names = repmat({''}, 1, length(seqs1));
diff_match_names = repmat({''}, 1, length(seqs1));
d01_cnt = zeros(1, length(seqs1));
t = cputime;
for j = range(1):min(length(seqs1), range(2))
	s1 = seqs1{j};
	for k = 1:length(seqs2)
		if mod(k, 1000)==0
			fprintf('seqs1 %i (%i) vs seqs2 %i (%i)\r', j, length(seqs1), k, length(seqs2));
		end
		s2 = seqs2{k};
		minimal_possible_distance = abs(length(s1)-length(s2))/max(length(s2), length(s2));
		if minimal_possible_distance>max(mindist(j), 0.05)
			continue
		end
		[matches, mismatches] = run_mex_align(s1, s2);
		distance = aa_dist(matches, length(s1), length(s2));
		if distance<mindist(j)
			mindist(j) = distance;
			match_names{j} = sprintf('%s,%s', names1{j}, names2{k});
		end	
		if distance<0.01
			d01_cnt(j) = d01_cnt(j)+1;
			d01_match_names{j} = sprintf('%s;%s,%s',d01_match_names{j}, names1{j}, names2{k});
		end
		if ~compare_names(names1{j}, names2{k}) && distance < mindist_diff(j)
			diff_match_names{j} = sprintf('%s;%s,%s',diff_match_names{j}, names1{j}, names2{k});
			mindist_diff(j) = distance;
		end
	end
end
runtime = cputime-t
save(sprintf('%s%i.mat',fn_out, range(1)), 'mindist', 'd01_cnt', 'match_names', 'd01_match_names', 'diff_match_names', 'mindist_diff')
return 

function ret = compare_names(name1, name2)
	items1 = separate(name1,'.');
	items2 = separate(name2,'.');
	items1 = separate(items1{1}, '_'); % remove Can_0_...
	items2 = separate(items2{1}, '_');
	ret = strcmp(items1{end}, items2{end});
return

function [matches, mismatches] = run_mex_align(aa1, aa2)
	score = mexalign(aa1, aa2);
	score_gap = -1;
	score_mismatch = -1;
	mismatches = (score-abs(length(aa1)-length(aa2))*score_gap)/score_mismatch;
	matches = min(length(aa1), length(aa2))-mismatches;
return




