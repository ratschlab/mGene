function genes = find_paralogs(genes,genome_info,thresh,min_len,max_exons)
% function genes = find_paralogs(genes,genome_info,thresh,min_len,max_exons)
%
% find all pairwise paralogs of genes, and label the genes.
% uses longest mRNA sequence greater than min_len,
% and BLAT at thresh% sequence identity.
% Ignores genes with more than max_exons exons.
% thresh defaults to 50
% min_len defaults to 100
% max_exons defaults to 100

if nargin<3,
  thresh = 0.5
end ;
if nargin<4,
  min_len = 100
end ;
if nargin<5,
  max_exons = 300
end ;


set(0,'RecursionLimit',300);
fname = genes2fasta(genes,min_len,max_exons,genome_info);

fprintf('BLAT with -minIdentity=%2.0f\n', thresh*100) ;
unix(sprintf('/fml/ag-raetsch/home/ong/bin/blat -minIdentity=%2.0f %s.fa %s.fa %s.psl -dots=1000 -noHead', ...
	     thresh*100, fname, fname, fname)) ;

similarity = psl2similarity(fname, length(genes), length(genes));

unix(sprintf('rm %s.fa %s.psl', fname, fname)) ;

fprintf('Cluster via single linkage above a certain threshold\n');
% cluster
sets = {} ;
for i=1:length(genes),
  sets{i} = i ;
end ;

changed = 1 ;
while(changed)
  changed = 0 ;
  for i=1:length(genes),
    % find homologuos sets
    idx = find(any(similarity(:,sets{i}),2)) ;
    % join sets with current
    %if length(idx)>1, keyboard ; end ;
    for j=1:length(idx)
      if ~any(idx(j)==sets{i}), 
        if ~isempty(sets{idx(j)}),
          changed = changed + 1 ; 
        end ;
        sets{i} = union(sets{idx(j)}, sets{i}) ;
        sets{idx(j)}=[] ; 
      end ;
    end ;
  end ;
  num_elem = 0 ; num_used = 0 ;
  for i=1:length(genes),
    num_elem = num_elem + length(sets{i}) ;
    num_used = num_used + (length(sets{i})>0) ;
  end ;
  assert(num_elem==length(genes)) ;
  [changed num_used]
end ;







% label the genes
for ix = 1:length(genes)
  genes(ix).has_paralog = 0;
end
for ix = 1:length(genes)
  if length(sets{ix}) > 1
    for ixp = sets{ix}
      genes(ixp).has_paralog = 1;
      genes(ixp).paralogs = sets{ix};
    end
  end    
end

