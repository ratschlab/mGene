function [names, matches, mismatches, q_gap_bases, t_gap_bases, all_items] = blat_wrapper(fasta1_seq, fasta2_seq, seq_type) 
% [names, matches, mismatches, q_gap_bases, t_gap_bases] =blat_wrapper(fasta1_seq, fasta2_seq) 

if nargin<3,
    seq_type='prot' ;
end ;

if nargin==0
	% test
	fasta1 = '~/tmp/blat_exp1.fasta'
	fasta2 = '~/tmp/blat_exp2.fasta'
	outfile = '~/tmp/align.psl';
	%outfile = '-';
end

fasta1=tempname ;
fasta2=tempname ;
outfile=tempname ;

fd=fopen(fasta1, 'w+') ;
write_fasta(fd, 'fasta1', fasta1_seq) ;
fclose(fd);

fd=fopen(fasta2, 'w+') ;
write_fasta(fd, 'fasta1', fasta2_seq) ;
fclose(fd);

% align all sequences in fasta file 1 to all sequences in fasta file 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blat_exe = '/fml/ag-raetsch/share/software/blat/x86_64/blat';
fprintf('call: %s -t=%s -q=%s %s %s %s\n', blat_exe, seq_type, seq_type, fasta1, fasta2, outfile);
[a b] = unix(sprintf('%s -t=prot -q=prot %s %s %s', blat_exe, fasta1, fasta2, outfile));

% parse output file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fd = fopen(outfile, 'r');

names=[] ;
matches=[] ;
mismatches=[] ;
q_gap_bases=[] ;
t_gap_bases=[] ;

header_skip = 1;
cnt = 0; 
all_items={} ;
while ~feof(fd)
	line = fgetl(fd);

	if length(line)>0 && line(1)=='-'
		header_skip=0; 
		continue
	end
	if header_skip
		continue
	end
	items = separate(line) ;
    all_items{end+1}=items ;
	cnt = cnt+1;
	matches(cnt) = str2num(items{1});
	mismatches(cnt) = str2num(items{2});
	q_gap_bases(cnt) = str2num(items{6});% gap bases in query
	t_gap_bases(cnt) = str2num(items{8});% gap bases in target
	names{cnt} = items{10};
end

fclose(fd);

unix(sprintf('rm %s %s %s', fasta1, fasta2, outfile)) ;
