function str = load_mrna(contig,strand,start,stop, genome_info)
% function str=load_mrna(contig,strand,start,stop, genome_info)
%
% Similar to load genomic, except check for consensus when doing multi blocks.
% Also assumes a global genome_info.

if nargin<5,
  global g_genome_info ;
  genome_info = g_genome_info ;
end ;

if length(start)>1 | length(stop)>1,
  assert(length(start)==length(stop)) ;
  assert(all([start(2:end)-stop(1:end-1)]>0));
  if strand == '+', 
    idx = 1:length(start) ; 
  else 
    idx = length(start):-1:1 ; 
  end ;
  str = '' ;
  for i=idx
    seq = load_mrna(contig,strand,start(i)-5,stop(i)+5) ;
    str = [str seq(6:end-5)] ;

    % check internal splicing sites
    if (i~=1 && strand == '+') || (i~=length(start) && strand == '-')
      assert(isequal(seq(4:5),'ag'))
    end ;
    if (i~=1 && strand == '-') || (i~=length(start) && strand == '+')
      assert(isequal(seq(end-4:end-3),'gt')||isequal(seq(end-4:end-3),'gc'))
    end ; 
  end
  return ;
end ;


if ischar(contig),
  contig_idx = find(ismember(genome_info.contig_names,contig)) ;
  assert(all(size(contig_idx)==1) & ~isempty(contig_idx)) ;
else
  contig_idx = contig ;
end ;
  
fname=genome_info.flat_fnames{contig_idx} ;
if genome_info.alphabet(1)=='a',
  NN='n' ;
else
  NN='N' ;
end ;

d=dir(fname) ; left_n = 0 ; right_n=0 ;
if isinf(stop), stop=d.bytes; end ;

if start<1 | stop>d.bytes,
  warning('boundary of contig reached, padding with "n"') ;
  if start<1,
    left_n = -(start-1) ;
    start = 1 ;
  end ;
  if stop>d.bytes,
    right_n = stop-d.bytes ;
  end ;
  if start>d.bytes && stop>d.bytes,
    str = char(NN*ones(1,stop-start)) ;
    return ;
  end ;
end ;

fd=fopen(fname,'r') ;
fseek(fd, start-1, -1) ;
str=char(fread(fd,stop-start+1, 'char=>char'))' ;
fclose(fd) ;
if left_n>0,
  str=[char(NN*ones(1,left_n)) str] ;
end ;
if right_n>0,
  str=[str char(NN*ones(1,right_n))] ;
end ;

assert(length(str)==stop-start+left_n+1)
let = unique(str) ;
if genome_info.alphabet(1)=='a' && ~isempty(setdiff(let,'acgtn')),
  warning(sprintf('invalid letter in genome! ("%s")', setdiff(let,'acgtn')));
  let = setdiff(let, 'acgtn') ;
  for c=let,
    str(str==c) = 'n';
  end ;
elseif genome_info.alphabet(1)=='A' && ~isempty(setdiff(let,'ACGTN')),
  warning(sprintf('invalid letter in genome! ("%s")', setdiff(let,'ACGTN')));
  let = setdiff(let, 'ACGTN') ;
  for c=let,
    str(str==c) = 'N';
  end ;
end

if strand=='-',
  str=reverse_complement(str) ;
end ;

return ;
