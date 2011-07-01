function colmap=compute_colocation_map(genes1, genes2, genome_info, prep_genes)
% colmap=compute_colocation_map(genes1, genes2, genome_info)
  

if nargin<4
  prep_genes =1;
end
    
%% prepare

CHR = genome_info.contig_names;
STR = '+-';
chr_maps = {} ;
for i=1:length(CHR)
  d=dir(genome_info.flat_fnames{i}) ;
  chr_maps{1,i}=uint32(zeros(1,d(1).bytes)) ;
  chr_maps{2,i}=uint32(zeros(1,d(1).bytes)) ;
end ;

%% prepare genes

if prep_genes
  if ~isfield(genes1,'strand')
    for g=1:length(genes1)
      genes1(g).strand = genes1(g).strands(1);
    end
  end
  

  if ~isfield(genes2,'strand')
    for g=1:length(genes2)
    genes2(g).strand=genes2(g).strands(1);
    end
  end
  
  for i=1:length(genes1) ;
    % start=inf; stop=0 ;
    start = genes1(i).start ;
    stop = genes1(i).stop ;
    for j=1:length(genes1(i).exons)
      if start>min(genes1(i).exons{j}(:)), tmp=genes1(i).exons{j}(:,1:2) ; start=min(tmp(:)) ; end ;
      if stop<max(genes1(i).exons{j}(:)), tmp=genes1(i).exons{j}(:,1:2) ; stop=max(tmp(:)) ; end ;
    end ;
    if ~isinf(start) && stop~=0 ;
      genes1(i).start=start ;
      genes1(i).stop=stop ;
    end ;
    genes1(i).chr_num =  strmatch(genes1(i).chr,CHR,'exact') ;
  end
  
  for i=1:length(genes2) ;
    % start=inf; stop=0 ;
    start = genes2(i).start ;
    stop = genes2(i).stop ;
    for j=1:length(genes2(i).exons)
      if start>min(genes2(i).exons{j}(:)), tmp=genes2(i).exons{j}(:,1:2) ; start=min(tmp(:)) ; end ;
      if stop<max(genes2(i).exons{j}(:)), tmp=genes2(i).exons{j}(:,1:2) ; stop=max(tmp(:)) ; end ;
    end ;
    if ~isinf(start) & stop~=0,
      genes2(i).start=max(1, start) ;
      genes2(i).stop=stop ;
    end ;
    genes2(i).chr_num =  strmatch(genes2(i).chr,CHR,'exact') ;
  end
end

%%% start colocation

chr_maps1 = fill_chr_map(genes1,chr_maps);

[chr_maps2, colmap] = fill_chr_map(genes2,chr_maps,chr_maps1);



function [chr_maps, colmap] = fill_chr_map(genes, chr_maps ,chr_maps1)

if isempty(genes)
  colmap={} ;
  return ;
end ;

colmap{length(genes)}={}; ; 

CHR = unique([genes.chr_num]);
STR = '+-';
num = 0;
chr_nums = [genes.chr_num];
strands = [genes.strand];
for c=CHR
  num=num+1;
  if (mod(num,100)==0), 
    fprintf('%i/%i\r',num,length(CHR))
  end
  for s=1:length(STR)
    idx = find(chr_nums==c&strands==STR(s));
    [starts,ii1] = sort([genes(idx).start]);
    stops = [genes(idx).stop];
    stops = stops(ii1);
    ids = idx(ii1);
    for i=1:length(ids)
      start = max(starts(i),1);
      stop = min(stops(i),length(chr_maps{s,c}));
      chr_maps{s,c}(start:stop) = ids(i) ;
      if exist('chr_maps1', 'var')
        colmap{ids(i)}= unique(chr_maps1{s,c}(start:stop)); 
        colmap{ids(i)}(colmap{ids(i)}==0)=[] ;
      end
    end
  end
end

fprintf('\n') ;