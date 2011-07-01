function [frac_exon, frac_intron, frac_5prime, frac_3prime,...
          nums, totals] = exon_counter(genes,fid)

genes = alt_const(genes);
idx_alt=[] ;
for i=1:length(genes)
  if genes(i).is_alt_spliced
    idx_alt=[idx_alt, i] ;
  end
end


for ix = 1:length(genes)
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};

  [dummy,exon_order] = sort(vertices(1,:),2,'ascend');
  if ~isequal(vertices(1,:),dummy)
    vertices = vertices(:,exon_order);
    edges = edges(exon_order,exon_order);
  end ;
  genes(ix).splicegraph{1} = vertices;
  genes(ix).splicegraph{2} = edges;
  
end




tot_skip = 0;
[idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt);
if ~isempty(idx_exon_skips)
  num_exon = size(unique([idx_exon_skips;exon_exon_skips(2,:)]','rows'),1);
else
  num_exon = 0;
end;
tot_skip = tot_skip + num_exon;
fprintf(1,'Number of unique single exon skips:\t\t\t\t%d\n',num_exon);
fprintf(fid,'Number of unique single exon skips:\t\t\t\t%d\n',num_exon);



[idx_xor_exons, exon_xor_exons] = detect_xorexons(genes, idx_alt) ;
if ~isempty(idx_xor_exons)
  num_xor = size(unique([idx_xor_exons;exon_xor_exons(2:3,:)]','rows'),1);
else
  num_xor = 0;
end
num_xor = num_xor*2;
tot_skip = tot_skip + num_xor;
fprintf(1,'Number of unique XOR exons :\t\t\t\t\t%d\n',num_xor);
fprintf(fid,'Number of unique XOR exons :\t\t\t\t\t%d\n',num_xor);





[idx_multiple_skips, exon_multiple_skips] = detect_multipleskips(genes, idx_alt) ;
num_unique = size(unique([idx_multiple_skips;exon_multiple_skips]','rows'),1);
%tot_skip = tot_skip + num_unique;
fprintf(1,'Number of unique multiple exon skips (not in total):\t\t\t%d\n',num_unique);
fprintf(fid,'Number of unique multiple exon skips (not in total):\t\t\t%d\n',num_unique);


fprintf(1,'Number of unique exon skips:\t\t\t\t\t%d\n',tot_skip);
fprintf(fid,'Number of unique exon skips:\t\t\t\t\t%d\n',tot_skip);








[idx_intron_reten,intron_intron_reten] = detect_intronreten(genes, idx_alt) ;
if ~isempty(idx_intron_reten)
  num_intron = size(unique([idx_intron_reten;intron_intron_reten(1:2,:)]','rows'),1);
else
  num_intron = 0;
end
fprintf(1,'Number of unique intron retentions:\t\t\t\t%d\n',num_intron);
fprintf(fid,'Number of unique intron retentions:\t\t\t\t%d\n',num_intron);




[idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] = detect_altprime(genes, idx_alt);
%num_5prime = size(unique([idx_alt_5prime;exon_alt_5prime(:).threeprimesite]','rows'),1);
%if ~isempty(idx_alt_5prime)
%  num_5prime = length([exon_alt_5prime(:).fiveprimesites]);
%else
%  num_5prime = 0;
%end
num_5prime = length(idx_alt_5prime);
fprintf(1,'Number of unique alternative 5 prime sites:\t\t\t%d\n',num_5prime);
fprintf(fid,'Number of unique alternative 5 prime sites:\t\t\t%d\n',num_5prime);
%num_3prime = size(unique([idx_alt_3prime;exon_alt_3prime(:).fiveprimesite]','rows'),1);
%if ~isempty(idx_alt_3prime)
%  num_3prime = length([exon_alt_3prime(:).threeprimesites]);
%else
%  num_3prime =0;
%end
num_3prime = length(idx_alt_3prime);
fprintf(1,'Number of unique alternative 3 prime sites:\t\t\t%d\n',num_3prime);
fprintf(fid,'Number of unique alternative 3 prime sites:\t\t\t%d\n',num_3prime);




cand_exons = 0;
cand_introns = 0;
cand_5prime = 0;
cand_3prime = 0;
for ix = 1:length(genes)
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  cand_intron_str = {};
  for exon_idx = 1:size(vertices,2)
    cur_edge_left = sum(edges(1:exon_idx,exon_idx));
    cur_edge_right = sum(edges(exon_idx:end,exon_idx));
    if cur_edge_left && cur_edge_right
      cand_exons = cand_exons + 1;
    end
    intron_map = zeros(1,size(edges,1));
    intron_map(exon_idx+1:end) = edges(exon_idx,exon_idx+1:end);
    for exon_idx2 = find(intron_map)
      cand_intron_str{end+1} = sprintf('%i\t%i',vertices(2,exon_idx),vertices(1,exon_idx2));
    end
    %cand_introns = cand_introns + sum(edges(exon_idx,exon_idx+1:end));
    if cur_edge_left
      cand_3prime = cand_3prime + 1;
    end
    if cur_edge_right
      cand_5prime = cand_5prime + 1;
    end
  end
  cand_introns = cand_introns + length(unique(cand_intron_str));
end
fprintf(1,'Number of candidates for exon skips:\t\t\t\t%d\n',cand_exons);
fprintf(fid,'Number of candidates for exon skips:\t\t\t\t%d\n',cand_exons);
fprintf(1,'Number of candidates for intron retentions:\t\t\t%d\n',cand_introns);
fprintf(fid,'Number of candidates for intron retentions:\t\t\t%d\n',cand_introns);
fprintf(1,'Number of candidates for 5 prime sites:\t\t\t\t%d\n',cand_5prime);
fprintf(fid,'Number of candidates for 5 prime sites:\t\t\t\t%d\n',cand_5prime);
fprintf(1,'Number of candidates for 3 prime sites:\t\t\t\t%d\n',cand_3prime);
fprintf(fid,'Number of candidates for 3 prime sites:\t\t\t\t%d\n',cand_3prime);



frac_exon = tot_skip/cand_exons;
frac_intron = num_intron/cand_introns;
frac_5prime = num_5prime/cand_5prime;
frac_3prime = num_3prime/cand_3prime;

nums.exon = tot_skip;
totals.exon = cand_exons;

nums.intron = num_intron ; 
totals.intron = cand_introns;

nums.prime5 = num_5prime ;
totals.prime5 = cand_5prime;

nums.prime3 = num_3prime ; 
totals.prime3 = cand_3prime;

