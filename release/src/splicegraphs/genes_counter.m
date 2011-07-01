function [frac_exon, frac_intron, frac_5prime, frac_3prime, nums, totals] = genes_counter(genes,fid)

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
uniq_idx = unique(idx_exon_skips);
num_exon = length(uniq_idx);
tot_skip = tot_skip + num_exon;
fprintf(1,'Number of unique single exon skips:\t\t\t\t%d\n',num_exon);
fprintf(fid,'Number of unique single exon skips:\t\t\t\t%d\n',num_exon);



[idx_xor_exons, exon_xor_exons] = detect_xorexons(genes, idx_alt) ;
uniq_idx = unique([uniq_idx,idx_xor_exons]);
num_xor = length(unique(idx_xor_exons));
tot_skip = tot_skip + num_xor;
fprintf(1,'Number of unique XOR exons :\t\t\t\t\t%d\n',num_xor);
fprintf(fid,'Number of unique XOR exons :\t\t\t\t\t%d\n',num_xor);





%[idx_multiple_skips, exon_multiple_skips] = detect_multipleskips(genes, idx_alt) ;
%uniq_idx = unique([uniq_idx,idx_multiple_skips]);
%num_unique = size(unique([idx_multiple_skips;exon_multiple_skips]','rows'),1);
%tot_skip = tot_skip + num_unique;
%fprintf(1,'Number of unique multiple exon skips:\t\t\t\t%d\n',num_unique);


fprintf(1,'Number of unique exon skips:\t\t\t\t\t%d\n',tot_skip);
fprintf(fid,'Number of unique exon skips:\t\t\t\t\t%d\n',tot_skip);








[idx_intron_reten,intron_intron_reten] = detect_intronreten(genes, idx_alt) ;
num_intron = length(unique(idx_intron_reten));
fprintf(1,'Number of unique intron retentions:\t\t\t\t%d\n',num_intron);
fprintf(fid,'Number of unique intron retentions:\t\t\t\t%d\n',num_intron);




[idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] = detect_altprime(genes, idx_alt);
num_5prime = length(unique(idx_alt_5prime));
fprintf(1,'Number of unique alternative 5 prime sites:\t\t\t%d\n',num_5prime);
fprintf(fid,'Number of unique alternative 5 prime sites:\t\t\t%d\n',num_5prime);
num_3prime = length(unique(idx_alt_3prime));
fprintf(1,'Number of unique alternative 5 prime sites:\t\t\t%d\n',num_3prime);
fprintf(fid,'Number of unique alternative 5 prime sites:\t\t\t%d\n',num_3prime);






frac_exon = tot_skip/length(genes);
frac_intron = num_intron/length(genes);
frac_5prime = num_5prime/length(genes);
frac_3prime = num_3prime/length(genes);

nums.exon = tot_skip;
totals.exon = length(genes);

nums.intron = num_intron ; 
totals.intron = length(genes);

nums.prime5 = num_5prime ;
totals.prime5 = length(genes);

nums.prime3 = num_3prime ; 
totals.prime3 = length(genes);

